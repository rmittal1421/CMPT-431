#include <iostream>
#include "core/utils.h"
#include "core/graph.h"

#include <mpi.h>

const int ROOT_PROCESSOR = 0;

uintV countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2, uintV u, uintV v)
{

    uintE i = 0, j = 0; // indexes for array1 and array2
    uintV count = 0;

    if (u == v)
        return count;

    while ((i < len1) && (j < len2))
    {
        if (array1[i] == array2[j])
        {
            if ((array1[i] != u) && (array1[i] != v))
            {
                count++;
            }
            else
            {
                // triangle with self-referential edge -> ignore
            }
            i++;
            j++;
        }
        else if (array1[i] < array2[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    return count;
}

void triangleCountSerial(Graph &g)
{
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken;
    timer t1;
    t1.start();
    for (uintV u = 0; u < n; u++)
    {
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++)
        {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(),
                                             u,
                                             v);
        }
    }

    // For every thread, print out the following statistics:
    // rank, edges, triangle_count, communication_time
    // 0, 17248443, 144441858, 0.000074
    // 1, 17248443, 152103585, 0.000020
    // 2, 17248443, 225182666, 0.000034
    // 3, 17248444, 185596640, 0.000022

    time_taken = t1.stop();

    // Print out overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

void triangleCountParallel(Graph &g, int& world_rank, int& P, int strategy) {
    double overall_time = 0.0, communication_time = 0.0;
    timer overall_timer, comm_timer;

    if(world_rank == ROOT_PROCESSOR) {
        overall_timer.start();
    }

    int start_vertex = 0, end_vertex = 0;
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

        if(i == world_rank) {
            break;
        }
    }

    long local_count = 0, global_count = 0;
    uint num_edges = 0;
    for(int u = start_vertex; u < end_vertex; u++) {
        uintE out_degree = g.vertices_[u].getOutDegree();
        num_edges += out_degree;
        for(uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            local_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                          g.vertices_[u].getInDegree(),
                                          g.vertices_[v].getOutNeighbors(),
                                          g.vertices_[v].getOutDegree(),
                                          u,
                                          v);
        }
    }

    // ---- Synchronization phase starts ----

    if(strategy == 1) {
        if(world_rank == ROOT_PROCESSOR) {
            comm_timer.start();

            global_count += local_count;
            int i;
            for(i = 1; i < P; i++) {
                long number;

                // Should I place it here or outside of the loop -- ASK!

                MPI_Recv(
                /* data         = */ &number, 
                /* count        = */ 1, 
                /* datatype     = */ MPI_LONG, 
                /* source       = */ i, 
                /* tag          = */ 0, 
                /* communicator = */ MPI_COMM_WORLD, 
                /* status       = */ MPI_STATUS_IGNORE);
                global_count += number;

            }

            communication_time += comm_timer.stop();
        } else {
            comm_timer.start();

            MPI_Send(
            /* data         = */ &local_count, 
            /* count        = */ 1, 
            /* datatype     = */ MPI_LONG, 
            /* destination  = */ ROOT_PROCESSOR, 
            /* tag          = */ 0, 
            /* communicator = */ MPI_COMM_WORLD);

            communication_time += comm_timer.stop();
        }
    } else if(strategy == 2) {
        long *sub_counts = NULL;
        if(world_rank == ROOT_PROCESSOR) {
            sub_counts = (long *) malloc(sizeof(long) * P);
        }

        comm_timer.start();

        // API call to gather data in sub_counts
        MPI_Gather(
            /*send_data       = */ &local_count,
            /*send_count      = */ 1,
            /*send_datatype   = */ MPI_LONG,
            /*recv_data       = */ sub_counts,
            /*recv_count      = */ 1,
            /*recv_datatype   = */ MPI_LONG,
            /*root            = */ ROOT_PROCESSOR,
            /*communicator    = */ MPI_COMM_WORLD
        );

        if(world_rank == ROOT_PROCESSOR) {
            int i;
            for(i = 0; i < P; i++) {
                global_count += sub_counts[i];
            }
        } 

        communication_time += comm_timer.stop();
    } else {
        comm_timer.start();

        // API call to gather data in sub_counts
        MPI_Reduce(
            /*send_data      = */ &local_count,
            /*recv_data      = */ &global_count,
            /*count          = */ 1,
            /*datatype       = */ MPI_LONG,
            /*op             = */ MPI_SUM,
            /*root           = */ ROOT_PROCESSOR,
            /*communicator   = */ MPI_COMM_WORLD
        );

        communication_time += comm_timer.stop();
    }

    // ---- Synchronization phase ends ----

    if(world_rank == ROOT_PROCESSOR) {
        overall_time = overall_timer.stop();
        
        std::cout << "rank, edges, triangle_count, communication_time" << std::endl;
        std::string process_info = std::to_string(world_rank) + ", "
                                    + std::to_string(num_edges) + ", " 
                                    + std::to_string(local_count) + ", "
                                    + std::to_string(communication_time) 
                                    + "\n";
        std::cout << process_info;

        // Print out overall statistics
        std::cout << "Number of triangles : " << global_count << "\n";
        std::cout << "Number of unique triangles : " << global_count / 3 << "\n";
        std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << overall_time << "\n";
    } else {
        std::string process_info = std::to_string(world_rank) + ", "
                                    + std::to_string(num_edges) + ", " 
                                    + std::to_string(local_count) + ", "
                                    + std::to_string(communication_time) 
                                    + "\n";
        std::cout << process_info;
    }
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("triangle_counting_serial", "Count the number of triangles using serial and parallel execution");
    options.add_options("custom", {
                                      {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                      {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                                  });

    auto cl_options = options.parse(argc, argv);
    uint strategy = cl_options["strategy"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

    // MPI initializing
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the processor
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::cout << std::fixed;

    if(world_rank == ROOT_PROCESSOR) {
        // Get the world size and print it out here
        std::cout << "World size : " << world_size << "\n";
        std::cout << "Communication strategy : " << strategy << "\n";
    }

    Graph g;
    g.readGraphFromBinary<int>(input_file_path);

    if(strategy == 0) {
        triangleCountSerial(g);
    } else {
        triangleCountParallel(g, world_rank, world_size, strategy);
    }

    MPI_Finalize();

    return 0;
}

