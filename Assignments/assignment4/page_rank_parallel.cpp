#include <iostream>
#include "core/utils.h"
#include "core/graph.h"

#include <mpi.h>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
#define PAGERANK_MPI_TYPE MPI_LONG
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
#define PAGERANK_MPI_TYPE MPI_FLOAT
typedef float PageRankType;
#endif

const int ROOT_PROCESSOR = 0;

void pageRankSerial(Graph &g, int max_iters)
{
    uintV n = g.n_;
    double time_taken;
    timer t1;
    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    t1.start();
    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    // -------------------------------------------------------------------
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
    // -------------------------------------------------------------------

    // For every thread, print the following statistics:
    // rank, num_edges, communication_time
    // 0, 344968860, 1.297778
    // 1, 344968860, 1.247763
    // 2, 344968860, 0.956243
    // 3, 344968880, 0.467028

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++)
    {
        sum_of_page_ranks += pr_curr[u];
    }
    time_taken = t1.stop();
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

void pageRankParallelS1(Graph &g, int max_iters, int world_rank, int P) {
    uintV n = g.n_;
    double time_taken = 0.0, communication_time = 0.0;
    timer overall_timer, comm_timer;

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    // if(world_rank == ROOT_PROCESSOR) {
    //     overall_timer.start();
    // }

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    uintV* start_indices = NULL;
    int start_vertex = 0, end_vertex = 0;
    for(int i = 0; i < P; i++) {

        if(i == ROOT_PROCESSOR) {
            start_indices = (uintV*) malloc(sizeof(uintV) * P);
        }

        start_vertex = end_vertex;
        long count = 0;

        while(end_vertex < g.n_) {
            count += g.vertices_[end_vertex].getOutDegree();
            end_vertex += 1;
            if(count >= g.m_/P) {
                break;
            }
        }

        if(i == ROOT_PROCESSOR) {
            start_indices[i] = start_vertex;
        }

        if(i == world_rank && i != ROOT_PROCESSOR) {
            break;
        }
    }

    if(world_rank == ROOT_PROCESSOR) {
        for(int i = 0; i < P; i++) {
            std::cout << start_indices[i] << " ";
        }
        std::cout<<std::endl;
    } 

    // Push based pagerank
    // -------------------------------------------------------------------

    // ----- Initializing required buffers ----

    PageRankType* buffer = NULL;
    PageRankType* starter = NULL;

    if(world_rank == ROOT_PROCESSOR) {
        // buffer = new PageRankType[n];
        buffer = (PageRankType*) malloc(sizeof(PageRankType) * n);
        starter = &pr_next[end_vertex];
    } else {
        buffer = (PageRankType*) malloc(sizeof(PageRankType) * (end_vertex - start_vertex));
        // buffer = new PageRankType[end_vertex - start_vertex];
    }

    // ---- Buffers initialization finish ----

    uint num_edges = 0;
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = start_vertex; u < end_vertex; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            num_edges += out_degree;
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }

        // ---- Synchronization phase 1 starts ----

        // If process != root, send values of all nodes to root
        // Root receives and adds them together
        // Sends it to all processes

        if(world_rank == ROOT_PROCESSOR) {

            // First receives from everyone and then send to everyone

            for(int i = 1; i < P; i++) {

                MPI_Recv(
                /* data         = */ buffer, 
                /* count        = */ n, 
                /* datatype     = */ PAGERANK_MPI_TYPE, 
                /* source       = */ i, 
                /* tag          = */ 0, 
                /* communicator = */ MPI_COMM_WORLD, 
                /* status       = */ MPI_STATUS_IGNORE);

                for(int j = 0; j < n; j++) {
                    pr_next[j] += buffer[j];
                }
            }

            printf("Received information from all other processors\n");

            for(int i = 1; i < P; i++) {
                int count;
                if(i != P - 1) {
                    count = start_indices[i+1] - start_indices[i];
                } else {
                    count = n - start_indices[i];
                }

                MPI_Send(
                /* data         = */ starter, 
                /* count        = */ count, 
                /* datatype     = */ PAGERANK_MPI_TYPE, 
                /* destination  = */ i, 
                /* tag          = */ 0, 
                /* communicator = */ MPI_COMM_WORLD);
            }

            printf("Sent information to all other processors from root\n");

        } else {

            // First send to root and then receive from root

            MPI_Send(
            /* data         = */ pr_next, 
            /* count        = */ n, 
            /* datatype     = */ PAGERANK_MPI_TYPE, 
            /* destination  = */ ROOT_PROCESSOR, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD);

            printf("Sent current values to root from processor : %d\n", world_rank);

            MPI_Recv(
            /* data         = */ buffer, 
            /* count        = */ (end_vertex - start_vertex), 
            /* datatype     = */ PAGERANK_MPI_TYPE, 
            /* source       = */ ROOT_PROCESSOR, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD, 
            /* status       = */ MPI_STATUS_IGNORE);

            printf("Recieved my updated values to next iteration : %d\n", world_rank);

        }

        // ---- Synchronization phase 1 ends ----

        for (uintV v = 0; v < n; v++)
        {
            if(v >= start_vertex && v < end_vertex) {
                pr_next[v] = PAGE_RANK(buffer[v - start_vertex]);

                // reset pr_curr for the next iteration
                pr_curr[v] = pr_next[v];
            }

            pr_next[v] = 0.0;
        }
    }

    // -------------------------------------------------------------------    

    PageRankType local_count = 0, sum_of_page_ranks = 0;
    for (uintV u = start_vertex; u < end_vertex; u++)
    {
        local_count += pr_curr[u];
    }

    // ---- Synchronization phase 2 starts ----

    MPI_Reduce(
    /*send_data      = */ &local_count,
    /*recv_data      = */ &sum_of_page_ranks,
    /*count          = */ 1,
    /*datatype       = */ PAGERANK_MPI_TYPE,
    /*op             = */ MPI_SUM,
    /*root           = */ ROOT_PROCESSOR,
    /*communicator   = */ MPI_COMM_WORLD);

    // ---- Synchronization phase 2 ends ----

    if(world_rank == ROOT_PROCESSOR) {
        std::cout << "Done from processor: " << world_rank << std::endl;
        std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
        // std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
    } else {
        std::cout << "Done from processor: " << world_rank << std::endl;
    }

    delete[] pr_curr;
    delete[] pr_next;
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("page_rank_push", "Calculate page_rank using serial and parallel execution");
    options.add_options("", {
                                {"nIterations", "Maximum number of iterations", cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                                {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                            });

    auto cl_options = options.parse(argc, argv);
    uint strategy = cl_options["strategy"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;

    // MPI initializing
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the processor
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(world_rank == ROOT_PROCESSOR) {
        // Get the world size and print it out here
        // std::cout << "World size : " << world_size << "\n"
        std::cout << "Communication strategy : " << strategy << "\n";
        std::cout << "Iterations : " << max_iterations << "\n";
    }

    Graph g;
    g.readGraphFromBinary<int>(input_file_path);

    switch (strategy)
    {
    case 0:
        pageRankSerial(g, max_iterations);
        break;
    case 1:
        pageRankParallelS1(g, max_iterations, world_rank, world_size);
        break;
    case 2:
        break;
    case 3:
        break;
    default:
        break;
    }

    MPI_Finalize();

    return 0;
}
