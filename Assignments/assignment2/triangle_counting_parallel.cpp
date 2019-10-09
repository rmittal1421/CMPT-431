#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <thread>
#include <future>
#include "core/utils.h"
#include "core/graph.h"
#include <string.h>
#include <atomic>

typedef std::atomic<uintV> AtomicInteger;

uintV countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2, uintV u, uintV v)
{

    uintE i = 0, j = 0; // indexes for array1 and array2
    uintV count = 0;
    while ((i < len1) && (j < len2))
    {
        if (array1[i] == array2[j])
        {
            if ((array1[i] != u) && (array1[j] != v))
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
    double time_taken = 0.0;
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
    time_taken = t1.stop();
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

void triangleCountParallelVertex(Graph &g, uint &n_workers)
{
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;

    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    std::future<uintV> futureList[n_workers];
    uintV start = 0, eachWorkload = n/n_workers, leftOver = n%n_workers;

    std::cout << "thread_id, triangle_count, time_taken\n";
    t1.start();

    // Process each edge <u,v>
    for(int i = 0; i < n_workers; i++) {
        uintV carryOver = (leftOver > 0) ? 1 : 0;
        futureList[i] = std::async (std::launch::async, [&](uintV s, uintV workload, int iteration){
            timer threadTimer;
            threadTimer.start();
            uintV local_count = 0;
            for (uintV u = s; u < (s + workload); u++)
            {
                // For each outNeighbor v, find the intersection of inNeighbor(u) and outNeighbor(v)
                uintE out_degree = g.vertices_[u].getOutDegree();
                for (uintE i = 0; i < out_degree; i++)
                {
                    uintV v = g.vertices_[u].getOutNeighbor(i);
                    local_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                                    g.vertices_[u].getInDegree(),
                                                    g.vertices_[v].getOutNeighbors(),
                                                    g.vertices_[v].getOutDegree(),
                                                    u,
                                                    v);
                }
            }
            std::string threadInfo = std::to_string(iteration) 
                                     + ", "
                                     + std::to_string(local_count) 
                                     + ", "
                                     + std::to_string(threadTimer.stop())
                                     + "\n";
            std::cout << threadInfo;
            return local_count;
        }, start, eachWorkload + carryOver, i);
        start += eachWorkload + carryOver;
        if(leftOver > 0) leftOver--;
    }

    for(auto &f: futureList) {
        triangle_count += f.get();
    }
    time_taken = t1.stop();
    
    // Print the overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

void triangleCountParallelEdge(Graph &g, uint &n_workers)
{
    uintV n = g.n_;
    uintE m = g.m_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;

    //First find all the edges and fill them in a vector.
    std::vector<intPair> edges;
    for(uintV u = 0; u < n; u++) {
        uintE out_degree = g.vertices_[u].getOutDegree();
        for(uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            edges.push_back({u, v});
        }
    }

    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    std::future<uintV> futureList[n_workers];
    uintV start = 0, eachWorkload = m/n_workers, leftOver = m%n_workers;

    std::cout << "thread_id, triangle_count, time_taken\n";
    t1.start();

    // Process each edge <u,v>
    for(int i = 0; i < n_workers; i++) {
        uintV carryOver = (leftOver > 0) ? 1 : 0;
        futureList[i] = std::async (std::launch::async, [&](uintV s, uintV workload, int iteration){
            timer threadTimer;
            threadTimer.start();
            uintV local_count = 0;
            for (uintV p = s; p < (s + workload); p++)
            {
                intPair pair = edges[p];
                uintV u = pair.first;
                uintV v = pair.second;
                local_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                                    g.vertices_[u].getInDegree(),
                                                    g.vertices_[v].getOutNeighbors(),
                                                    g.vertices_[v].getOutDegree(),
                                                    u,
                                                    v);
            }
            std::string threadInfo = std::to_string(iteration) 
                                     + ", "
                                     + std::to_string(local_count) 
                                     + ", "
                                     + std::to_string(threadTimer.stop())
                                     + "\n";
            std::cout << threadInfo;
            return local_count;
        }, start, eachWorkload + carryOver, i);
        start += eachWorkload + carryOver;
        if(leftOver > 0) leftOver--;
    }

    for(auto &f: futureList) {
        triangle_count += f.get();
    }
    time_taken = t1.stop();
    
    // Print the overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

void triangleCountParallelDynamic(Graph &g, uint &n_workers)
{
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;
    AtomicInteger index(0);

    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    std::future<uintV> futureList[n_workers];

    std::cout << "thread_id, triangle_count, time_taken\n";
    t1.start();

    // Process each edge <u,v>
    for(int i = 0; i < n_workers; i++) {
        futureList[i] = std::async (std::launch::async, [&](int iteration){
            timer threadTimer;
            threadTimer.start();
            uintV local_count = 0;
            while(true)
            {
                //Get the vertex to be processed next.
                uintV u = index.fetch_add(1);
                // while(!index.compare_exchange_weak(u, u + 1));
                if(u >= n) break;

                uintE out_degree = g.vertices_[u].getOutDegree();
                for (uintE i = 0; i < out_degree; i++)
                {
                    uintV v = g.vertices_[u].getOutNeighbor(i);
                    local_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                                    g.vertices_[u].getInDegree(),
                                                    g.vertices_[v].getOutNeighbors(),
                                                    g.vertices_[v].getOutDegree(),
                                                    u,
                                                    v);
                }
            }
            std::string threadInfo = std::to_string(iteration) 
                                     + ", "
                                     + std::to_string(local_count) 
                                     + ", "
                                     + std::to_string(threadTimer.stop())
                                     + "\n";
            std::cout << threadInfo;
            return local_count;
        }, i);
    }

    for(auto &f: futureList) {
        triangle_count += f.get();
    }
    time_taken = t1.stop();
    
    // Print the overall statistics
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("triangle_counting_serial", "Count the number of triangles using serial and parallel execution");
    options.add_options("custom", {
                                      {"nWorkers", "Number of workers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                                      {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                      {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/assignment1/input_graphs/roadNet-CA")},
                                  });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    switch (strategy)
    {
    case 0:
        std::cout << "\nSerial\n";
        triangleCountSerial(g);
        break;
    case 1:
        std::cout << "\nVertex-based work partitioning\n";
        triangleCountParallelVertex(g, n_workers);
        break;
    case 2:
        std::cout << "\nEdge-based work partitioning\n";
        triangleCountParallelEdge(g, n_workers);
        break;
    case 3:
        std::cout << "\nDynamic work partitioning\n";
        triangleCountParallelDynamic(g, n_workers);
        break;
    default:
        break;
    }

    return 0;
}
