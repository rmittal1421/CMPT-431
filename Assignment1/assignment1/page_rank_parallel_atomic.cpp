#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "utils.h"
#include "graph.h"

#include <thread>
#include <vector>
#include <mutex>
#include <atomic>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
typedef std::atomic<PageRankType> AtomicPageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
typedef std::atomic<float> AtomicPageRankType;
#endif

void threadFunction() {
    std::cout<<"I am inside the thread function"<<std::endl;
}

void pageRankSerial(Graph &g, int max_iters)
{
    uintV n = g.n_;
    std::cout<<"Number of nodes: " << n << std::endl;

    PageRankType *pr_curr = new PageRankType[n];
    // PageRankType *pr_next = new PageRankType[n];
    AtomicPageRankType *pr_next = new AtomicPageRankType[n];

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    //Lets take the number of threads to be 4 to start.
    std::thread threadList[4];
    t1.start();

    // std::mutex mtx;
    CustomBarrier myBarrier(4);
    std::vector<std::mutex> mutex_vector(n);

    uintV start = 0, eachWorkload = n/4, leftOver = n % 4;
    // std::cout<<start<<" "<<eachWorkload<<std::endl;
    for(int i = 0; i < 4; i++) {
        uintV carryOver = (leftOver > 0) ? 1 : 0;
        threadList[i] = std::thread ([&](uintV s, uintV workload) {
            for (int iter = 0; iter < max_iters; iter++) {
                // for each vertex 'u', process all its outNeighbors 'v'
                for (uintV u = s; u < (s + workload); u++)
                {
                    uintE out_degree = g.vertices_[u].getOutDegree();
                    for (uintE i = 0; i < out_degree; i++)
                    {
                        uintV v = g.vertices_[u].getOutNeighbor(i);

                        // std::lock_guard<std::mutex> lock(mutex_vector[v]);
                        // pr_next[v] = pr_next[v] + (pr_curr[u] / out_degree);
                        PageRankType current_val = pr_next[v];
                        // do {
                        //     current_val = pr_next[v];   
                        // }
                        while(!pr_next[v].compare_exchange_strong(current_val, current_val + (pr_curr[u]/out_degree), std::memory_order_release,
                                        std::memory_order_relaxed));
                    }
                }

                myBarrier.wait();

                for (uintV v = s; v < (s + workload); v++)
                {
                    pr_next[v] = PAGE_RANK(pr_next[v]);

                    // reset pr_curr for the next iteration
                    pr_curr[v] = pr_next[v];
                    pr_next[v] = 0.0;
                }

                myBarrier.wait();
            }
        }, start, eachWorkload + carryOver);
        start += eachWorkload + carryOver;
        if(leftOver > 0) leftOver--;
    }

    for(auto &t : threadList) {
        t.join();
        std::cout<<"Joined a thread"<<std::endl;
    }

    // for (int iter = 0; iter < max_iters; iter++)
    // {
    //     // for each vertex 'u', process all its outNeighbors 'v'
    //     for (uintV u = 0; u < n; u++)
    //     {
    //         uintE out_degree = g.vertices_[u].getOutDegree();
    //         for (uintE i = 0; i < out_degree; i++)
    //         {
    //             uintV v = g.vertices_[u].getOutNeighbor(i);
    //             pr_next[v] += (pr_curr[u] / out_degree);
    //         }
    //     }
    //     for (uintV v = 0; v < n; v++)
    //     {
    //         pr_next[v] = PAGE_RANK(pr_next[v]);

    //         // reset pr_curr for the next iteration
    //         pr_curr[v] = pr_next[v];
    //         pr_next[v] = 0.0;
    //     }
    // }
    time_taken = t1.stop();
    // -------------------------------------------------------------------
    // std::cout << "thread_id, time_taken\n";
    // Print the above statistics for each thread
    // Example output for 2 threads:
    // thread_id, time_taken
    // 0, 0.12
    // 1, 0.12

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++){
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";

    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    // delete[] pr_curr;
    // delete[] pr_next;
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("page_rank_push", "Calculate page_rank using serial and parallel execution");
    options.add_options("", {
                                {"nWorkers", "Number of workers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                                {"nIterations", "Maximum number of iterations", cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                                {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    std::cout << input_file_path << std::endl;
    g.read_graph_from_binary<int>(input_file_path);
    std::cout<<"Read the graph"<<std::endl;
    std::cout << "Created graph\n";
    pageRankSerial(g, max_iterations);

    return 0;
}
