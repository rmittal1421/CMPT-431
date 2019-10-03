#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <thread>
#include "utils.h"
#include "graph.h"

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif
void thread_helper_pageRank(Graph &g,uintV subset_start,uintV subset_end,PageRankType *pr_curr,PageRankType *pr_next,int max_iters,CustomBarrier &my_barrier,std::mutex &lock){
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = subset_start; u < subset_end; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                lock.lock();
                pr_next[v] += (pr_curr[u] / out_degree);
                lock.unlock();
            }
        }
        my_barrier.wait();
        for (uintV v = subset_start; v < subset_end; v++)
        {
            pr_next[v] = PAGE_RANK(pr_next[v]);
            
            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
        my_barrier.wait();
    }
    
}
void pageRankSerial(Graph &g, int max_iters, uint no_of_threads)
{
    uintV n = g.n_;

    std::mutex lock;
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
    std::thread t[no_of_threads];
    
    CustomBarrier my_barrier(no_of_threads);
    
    uintV no_of_vertices_for_last_thread = n/no_of_threads + (n%no_of_threads) ;
    
    int i = 0 ;uintV subset_end = n/no_of_threads;uintV subset_start = 0;
    
    for (; i < no_of_threads - 1 ; i++) {
        t[i] = std::thread(thread_helper_pageRank,std::ref(g),subset_start,subset_end,std::ref(pr_curr),std::ref(pr_next),max_iters,std::ref(my_barrier),std::ref(lock));
        subset_start = subset_end;
        subset_end += n/no_of_threads;
    }
    
    subset_end = no_of_vertices_for_last_thread + subset_start;
    
    t[i] = std::thread(thread_helper_pageRank,std::ref(g),subset_start,subset_end,std::ref(pr_curr),std::ref(pr_next),max_iters,std::ref(my_barrier),std::ref(lock));
    for (int j =0 ; j < no_of_threads; j++) {
        
            // std::cout<<j;
            t[j].join();
            // std::cout<<", "<<t1.total()<<"\n";
        
    }

    // for(auto &thread: t) {
    //     thread.join();
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
    delete[] pr_curr;
    delete[] pr_next;
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
    g.read_graph_from_binary<int>(input_file_path);
    std::cout << "Created graph\n";
    pageRankSerial(g, max_iterations, n_workers);

    return 0;
}
