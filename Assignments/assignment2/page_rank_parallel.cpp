#include <iostream>
#include <iomanip>
#include <thread>
#include <stdlib.h>
#include "core/utils.h"
#include "core/graph.h"

#include <vector>
#include <mutex>
#include <atomic>
#include <string.h>

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
typedef std::atomic<PageRankType> AtomicPageRankType;
#endif

void pageRankSerial(Graph &g, int max_iters)
{
    uintV n = g.n_;

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

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
    t1.start();
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = 0; u < n; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }
        for (uintV v = 0; v < n; v++)
        {
            pr_next[v] = PAGE_RANK(pr_next[v]);

            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++)
    {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

void pageRankParallelVertex(Graph &g, int max_iters, uint &n_workers)
{
    uintV n = g.n_;

    PageRankType *pr_curr = new PageRankType[n];
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
    std::thread threadList[n_workers];
    CustomBarrier myBarrier(n_workers);
    std::vector<std::mutex> mutex_vector(n);
    uintV start = 0, eachWorkload = n/n_workers, leftOver = n % n_workers;

    std::cout << "thread_id, time_taken\n";
    t1.start();

    for(int i = 0; i < n_workers; i++) {
        uintV carryOver = (leftOver > 0) ? 1 : 0;
        threadList[i] = std::thread ([&](uintV s, uintV workload, int iteration) {
            timer threadTimer;
            threadTimer.start();

            for (int iter = 0; iter < max_iters; iter++) {
                // for each vertex 'u', process all its outNeighbors 'v'
                for (uintV u = s; u < (s + workload); u++)
                {
                    uintE out_degree = g.vertices_[u].getOutDegree();
                    for (uintE i = 0; i < out_degree; i++)
                    {
                        uintV v = g.vertices_[u].getOutNeighbor(i);

                        PageRankType current_val = pr_next[v];
                        while(!pr_next[v].compare_exchange_weak(current_val, current_val + (pr_curr[u]/out_degree)));
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
            std::string threadInfo = std::to_string(iteration) + ", " + std::to_string(threadTimer.stop()) + "\n";
            std::cout << threadInfo;
        }, start, eachWorkload + carryOver, i);
        start += eachWorkload + carryOver;
        if(leftOver > 0) leftOver--;
    }

    for(auto &t : threadList) {
        t.join();
    }

    time_taken = t1.stop();
    // -------------------------------------------------------------------
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
    delete[] pr_curr;
    delete[] pr_next;
}

void pageRankParallelEdge(Graph &g, int max_iters, uint &n_workers)
{
    uintV n = g.n_;
    uintE m = g.m_;

    timer t1;
    t1.start();
    std::vector<intPair> range_of_vertices(n_workers);
    uint max_workload = m/n_workers, current_sum = 0, current_thread = 0;
    int start = 0;
    for(uintV u = 0; u < n; u++) {
        if(current_thread == n_workers - 1) {
            range_of_vertices[current_thread] = {start, n-1};
            break;
        }
        uintE out_degree = g.vertices_[u].getOutDegree();
        if(current_sum == 0) current_sum += out_degree;
        else if(current_sum + out_degree <= max_workload) current_sum += out_degree;
        else {
            current_sum = 0;
            range_of_vertices[current_thread] = {start, u - 1};
            start = u;
            u--;
            current_thread++;
        }
    }

    PageRankType *pr_curr = new PageRankType[n];
    AtomicPageRankType *pr_next = new AtomicPageRankType[n];

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    // timer t1;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    std::thread threadList[n_workers];
    CustomBarrier myBarrier(n_workers);

    std::cout << "thread_id, time_taken\n";
    // t1.start();

    for(int i = 0; i < n_workers; i++) {
        threadList[i] = std::thread ([&](int iteration) {
            timer threadTimer;
            threadTimer.start();

            for (int iter = 0; iter < max_iters; iter++) {
                // for each vertex 'u', process all its outNeighbors 'v'
                auto range_values = range_of_vertices[iteration];
                for (uintV u = range_values.first; u <= range_values.second; u++)
                {
                    uintE out_degree = g.vertices_[u].getOutDegree();
                    for (uintE i = 0; i < out_degree; i++)
                    {
                        uintV v = g.vertices_[u].getOutNeighbor(i);

                        PageRankType current_val = pr_next[v];
                        while(!pr_next[v].compare_exchange_weak(current_val, current_val + (pr_curr[u]/out_degree)));
                    }
                }

                myBarrier.wait();
                for (uintV v = range_values.first; v <= range_values.second; v++)
                {
                    pr_next[v] = PAGE_RANK(pr_next[v]);

                    // reset pr_curr for the next iteration
                    pr_curr[v] = pr_next[v];
                    pr_next[v] = 0.0;
                }

                myBarrier.wait();
            }
            std::string threadInfo = std::to_string(iteration) + ", " + std::to_string(threadTimer.stop()) + "\n";
            std::cout << threadInfo;
        }, i);
    }

    for(auto &t : threadList) {
        t.join();
    }

    time_taken = t1.stop();
    // -------------------------------------------------------------------
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
    delete[] pr_curr;
    delete[] pr_next;
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("page_rank_push", "Calculate page_rank using serial and parallel execution");
    options.add_options("", {
                                {"nWorkers", "Number of workers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                                {"nIterations", "Maximum number of iterations", cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                                {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                {"granularity", "Granularity to be used", cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
                                {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    uint granularity = cl_options["granularity"].as<uint>();    
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";
    std::cout << "Granularity : " << granularity << "\n";
    std::cout << "Iterations : " << max_iterations << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";
    switch (strategy)
    {
    case 0:
        std::cout << "\nSerial\n";
        pageRankSerial(g, max_iterations);
        break;
    case 1:
        std::cout << "\nVertex-based work partitioning\n";
        pageRankParallelVertex(g, max_iterations, n_workers);
        break;
    case 2:
        std::cout << "\nEdge-based work partitioning\n";
        pageRankParallelEdge(g, max_iterations, n_workers);
        break;
    case 3:
        std::cout << "\nDynamic work partitioning\n";
        break;
    default:
        break;
    }

    return 0;
}
