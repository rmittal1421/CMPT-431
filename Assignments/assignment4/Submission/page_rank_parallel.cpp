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

    if(world_rank == ROOT_PROCESSOR) {
        overall_timer.start();
    }

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    uintV* start_indices = NULL;
    int start_vertex = 0, end_vertex = 0;

    if(world_rank == ROOT_PROCESSOR) {
        start_indices = (uintV*) malloc(sizeof(uintV) * P);
    }

    for(int i = 0; i < P; i++) {

        start_vertex = end_vertex;
        long count = 0;

        while(end_vertex < g.n_) {
            count += g.vertices_[end_vertex].getOutDegree();
            end_vertex += 1;
            if(count >= g.m_/P) {
                break;
            }
        }

        if(world_rank == ROOT_PROCESSOR) {
            start_indices[i] = start_vertex;
        }

        if(i == world_rank && i != ROOT_PROCESSOR) {
            break;
        }
    }

    if(world_rank == ROOT_PROCESSOR) {
        // Assuming root processor is not the last one!

        // Reset the start index for root processor
        start_vertex = start_indices[ROOT_PROCESSOR];

        // End vertex should be the start vertex of next processor
        end_vertex = (world_rank != P - 1) ? start_indices[world_rank + 1] : n;
    }

    // Push based pagerank
    // -------------------------------------------------------------------

    // ----- Initializing required buffers ----

    PageRankType* buffer = NULL;

    if(world_rank == ROOT_PROCESSOR) {
        buffer = (PageRankType*) malloc(sizeof(PageRankType) * n);
    } else {
        buffer = (PageRankType*) malloc(sizeof(PageRankType) * (end_vertex - start_vertex));
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

        comm_timer.start();

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

            for(int i = 1; i < P; i++) {
                int count;
                count = ((i != P - 1) ? start_indices[i+1] : n) - start_indices[i];
                // if(i != P - 1) {
                //     count = start_indices[i+1] - start_indices[i];
                // } else {
                //     count = n - start_indices[i];
                // }

                MPI_Send(
                /* data         = */ &pr_next[start_indices[i]], 
                /* count        = */ count, 
                /* datatype     = */ PAGERANK_MPI_TYPE, 
                /* destination  = */ i, 
                /* tag          = */ 0, 
                /* communicator = */ MPI_COMM_WORLD);
            }

        } else {

            // First send to root and then receive from root

            MPI_Send(
            /* data         = */ pr_next, 
            /* count        = */ n, 
            /* datatype     = */ PAGERANK_MPI_TYPE, 
            /* destination  = */ ROOT_PROCESSOR, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD);

            MPI_Recv(
            /* data         = */ buffer, 
            /* count        = */ (end_vertex - start_vertex), 
            /* datatype     = */ PAGERANK_MPI_TYPE, 
            /* source       = */ ROOT_PROCESSOR, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD, 
            /* status       = */ MPI_STATUS_IGNORE);

        }

        communication_time += comm_timer.stop();

        // ---- Synchronization phase 1 ends ----

        for (uintV v = 0; v < n; v++)
        {
            if(v >= start_vertex && v < end_vertex) {
                if(world_rank == ROOT_PROCESSOR) pr_next[v] = PAGE_RANK(pr_next[v - start_vertex]);
                else pr_next[v] = PAGE_RANK(buffer[v - start_vertex]);

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

    std::string process_info = std::to_string(world_rank) + ", "
                                + std::to_string(num_edges) + ", " 
                                + std::to_string(communication_time) 
                                + "\n";
    std::cout << process_info;

    if(world_rank == ROOT_PROCESSOR) {
        time_taken += overall_timer.stop();

        // Print out overall statistics
        std::string overall_sum_information = "Sum of page rank : " + std::to_string(sum_of_page_ranks) + "\n";
        std::cout << overall_sum_information;
        
        std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
        free(start_indices);
    }

    delete[] pr_curr;
    delete[] pr_next;

    // Free the initialized buffer as well
    free(buffer);
}

void pageRankParallelS2(Graph &g, int max_iters, int world_rank, int P) {
    uintV n = g.n_;
    double time_taken = 0.0, communication_time = 0.0;
    timer overall_timer, comm_timer;

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    if(world_rank == ROOT_PROCESSOR) {
        overall_timer.start();
    }

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    uintV* start_indices = NULL;
    int* counts = NULL;
    int start_vertex = 0, end_vertex = 0;

    if(world_rank == ROOT_PROCESSOR) {
        start_indices = (uintV*) malloc(sizeof(uintV) * P);
        counts = (int*) malloc(sizeof(int) * P);
    }

    for(int i = 0; i < P; i++) {

        start_vertex = end_vertex;
        long count = 0;

        while(end_vertex < g.n_) {
            count += g.vertices_[end_vertex].getOutDegree();
            end_vertex += 1;
            if(count >= g.m_/P) {
                break;
            }
        }

        if(world_rank == ROOT_PROCESSOR) {
            start_indices[i] = start_vertex;
        }

        if(i == world_rank && i != ROOT_PROCESSOR) {
            break;
        }
    }

    if(world_rank == ROOT_PROCESSOR) {
        // Reset the start index for root processor
        start_vertex = start_indices[ROOT_PROCESSOR];

        // End vertex should be the start vertex of next processor
        end_vertex = (world_rank != P - 1) ? start_indices[world_rank + 1] : n;

        for(int i = 0; i < P; i++) {
            counts[i] = ((i != P - 1) ? start_indices[i+1] : n) - start_indices[i];
            // if(i != P - 1) {
            //     counts[i] = start_indices[i + 1] - start_indices[i];
            // } else {
            //     counts[i] = n - start_indices[i];
            // }
        }
    }

    // Push based pagerank
    // -------------------------------------------------------------------

    // ----- Initializing required buffers ----

    PageRankType* buffer = NULL;
    PageRankType* overall_buffer = NULL;

    if(world_rank == ROOT_PROCESSOR) {
        overall_buffer = (PageRankType*) malloc(sizeof(PageRankType) * n);
    } 
    buffer = (PageRankType*) malloc(sizeof(PageRankType) * (end_vertex - start_vertex));

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

       // Here, via reduce, everyone sends out their values of pr_next to root
       // Root gets the addition of them inside a buffer
       // Then scatters to everyone back about their respective nodes

       comm_timer.start();

        MPI_Reduce(
        /*send_data      = */ pr_next,
        /*recv_data      = */ overall_buffer,
        /*count          = */ n,
        /*datatype       = */ PAGERANK_MPI_TYPE,
        /*op             = */ MPI_SUM,
        /*root           = */ ROOT_PROCESSOR,
        /*communicator   = */ MPI_COMM_WORLD);

        MPI_Scatterv(
        /*sendbuf        = */ overall_buffer,
        /*sendcounts     = */ counts,
        /*displs         = */ start_indices,
        /*sendtype       = */ PAGERANK_MPI_TYPE,
        /*recvbuf        = */ buffer,
        /*recvcount      = */ (end_vertex - start_vertex),
        /*recvtype       = */ PAGERANK_MPI_TYPE,
        /*root           = */ ROOT_PROCESSOR,
        /*comm           = */ MPI_COMM_WORLD);

        communication_time += comm_timer.stop();

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

    std::string process_info = std::to_string(world_rank) + ", "
                                + std::to_string(num_edges) + ", " 
                                + std::to_string(communication_time) 
                                + "\n";
    std::cout << process_info;

    if(world_rank == ROOT_PROCESSOR) {
        time_taken += overall_timer.stop();

        // Print out overall statistics
        std::string overall_sum_information = "Sum of page rank : " + std::to_string(sum_of_page_ranks) + "\n";
        std::cout << overall_sum_information;

        std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";

        // Free buffers initialized just for root processor
        free(start_indices);
        free(counts);
        free(overall_buffer);
    }

    delete[] pr_curr;
    delete[] pr_next;

    // Free the initialized buffer for every processor
    free(buffer);
}

void pageRankParallelS3(Graph &g, int max_iters, int world_rank, int P) {
    uintV n = g.n_;
    double time_taken = 0.0, communication_time = 0.0;
    timer overall_timer, comm_timer;

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    if(world_rank == ROOT_PROCESSOR) {
        overall_timer.start();
    }

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    uintV* start_indices = NULL;
    int* counts = NULL;
    int start_vertex = 0, end_vertex = 0;

    start_indices = (uintV*) malloc(sizeof(uintV) * P);
    counts = (int*) malloc(sizeof(int) * P);

    for(int i = 0; i < P; i++) {

        start_vertex = end_vertex;
        long count = 0;

        while(end_vertex < g.n_) {
            count += g.vertices_[end_vertex].getOutDegree();
            end_vertex += 1;
            if(count >= g.m_/P) {
                break;
            }
        }

        start_indices[i] = start_vertex;
    }

    start_vertex = start_indices[world_rank];
    end_vertex = (world_rank != P - 1) ? start_indices[world_rank + 1] : n; 

    for(int i = 0; i < P; i++) {
        counts[i] = ((i != P - 1) ? start_indices[i+1] : n) - start_indices[i];
        // if(i != P - 1) {
        //     counts[i] = start_indices[i + 1] - start_indices[i];
        // } else {
        //     counts[i] = n - start_indices[i];
        // }
    }

    // Push based pagerank
    // -------------------------------------------------------------------

    // ----- Initializing required buffers ----

    PageRankType* buffer = NULL;
    buffer = (PageRankType*) malloc(sizeof(PageRankType) * (end_vertex - start_vertex));

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

       // Here all the processes get the updates values from other processes
       // And also contribute to other process's nodes in reduce

       comm_timer.start();

       for(int i = 0; i < P; i++) {
           MPI_Reduce(
            /*send_data      = */ &pr_next[start_indices[i]],
            /*recv_data      = */ buffer,
            /*count          = */ counts[i],
            /*datatype       = */ PAGERANK_MPI_TYPE,
            /*op             = */ MPI_SUM,
            /*root           = */ i,
            /*communicator   = */ MPI_COMM_WORLD);
       }

       communication_time += comm_timer.stop();

        // ---- Synchronization phase 1 ends ----

        for (uintV v = 0; v < n; v++)
        {
            if(v >= start_vertex && v < end_vertex) {
                pr_next[v] = PAGE_RANK(buffer[v - start_vertex]);
                // if(world_rank == ROOT_PROCESSOR) pr_next[v] = PAGE_RANK(pr_next[v - start_vertex]);
                // else pr_next[v] = PAGE_RANK(buffer[v - start_vertex]);

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

    std::string process_info = std::to_string(world_rank) + ", "
                                + std::to_string(num_edges) + ", " 
                                + std::to_string(communication_time) 
                                + "\n";
    std::cout << process_info;

    if(world_rank == ROOT_PROCESSOR) {
        time_taken += overall_timer.stop();

        // Print out overall statistics
        std::string overall_sum_information = "Sum of page rank : " + std::to_string(sum_of_page_ranks) + "\n";
        std::cout << overall_sum_information;

        std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
    }

    delete[] pr_curr;
    delete[] pr_next;

    // Free all the buffers initialized
    free(start_indices);
    free(counts);
    free(buffer);
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

    // MPI initializing
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the processor
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    Graph g;
    g.readGraphFromBinary<int>(input_file_path);

    if(world_rank == ROOT_PROCESSOR) {
        #ifdef USE_INT
            std::cout << "Using INT\n";
        #else
            std::cout << "Using FLOAT\n";
        #endif
            std::cout << std::fixed;

        // Get the world size and print it out here
        printf("World size : %d\n", world_size);
        printf("Communication strategy : %u\n", strategy);
        printf("Iterations : %u\n", max_iterations);
        printf("rank, num_edges, communication_time\n");
    }


    switch (strategy)
    {
    case 0:
        if(world_rank == ROOT_PROCESSOR) {
            pageRankSerial(g, max_iterations);
        }
        break;
    case 1:
        pageRankParallelS1(g, max_iterations, world_rank, world_size);
        break;
    case 2:
        pageRankParallelS2(g, max_iterations, world_rank, world_size);
        break;
    case 3:
        pageRankParallelS3(g, max_iterations, world_rank, world_size);
        break;
    default:
        break;
    }

    MPI_Finalize();

    return 0;
}
