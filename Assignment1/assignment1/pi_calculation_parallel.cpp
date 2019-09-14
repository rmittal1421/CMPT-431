#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "utils.h"

#include <vector>
#include <thread>
#include <future>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "12345678"

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed)
{
    return ((double)rand_r(random_seed)) / c_const;
}

uint get_points_in_circle(uint n, uint random_seed)
{
    std::cout<<"Thread id of async function=> "<<std::this_thread::get_id()<<std::endl;
    uint circle_count = 0;
    double x_coord, y_coord;
    for (uint i = 0; i < n; i++)
    {
        x_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
        y_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
        if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
            circle_count++;
    }
    return circle_count;
}

void piCalculation(uint n)
{
    timer serial_timer;
    double time_taken = 0.0;
    uint random_seed = 1;

    serial_timer.start();
    // Create threads and distribute the work across T threads

    //Taking 4 futures as of now.
    std::cout<<"Main thread => "<<std::this_thread::get_id()<<std::endl; 
    // std::vector<std::future<uint>> futureList(4);
    std::future<uint> futureList[4];
    // std::future<uint> fut = async(std::launch::async, get_points_in_circle, n, random_seed);
    // std::cout<<fut.get()<<std::endl;

    uint circle_points = 0;
    for(int i = 0; i < 4; i++) {
        futureList[i] = std::async(std::launch::async, get_points_in_circle, n/4, random_seed);
    }

    for(auto &f : futureList) 
        circle_points += f.get();
    std::cout<<"Total circle points allocated are: "<<circle_points<<std::endl;
    // -------------------------------------------------------------------
    // uint circle_points = get_points_in_circle(n, random_seed);
    double pi_value = 4.0 * (double)circle_points / (double)n;
    // -------------------------------------------------------------------
    time_taken = serial_timer.stop();
    
    // std::cout << "thread_id, points_generated, circle_points, time_taken\n"; 
    // // Print the above statistics for each thread
    // // Example output for 2 threads:
    // // thread_id, points_generated, circle_points, time_taken
    // // 1, 100, 90, 0.12
    // // 0, 100, 89, 0.12
    
    // // Print the overall statistics
    // std::cout << "Total points generated : " << n << "\n";
    // std::cout << "Total points in circle : " << circle_points << "\n";
    std::cout << "Result : " << std::setprecision(VAL_PRECISION) << pi_value << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
}

int main(int argc, char *argv[])
{
    // Initialize command line arguments
    cxxopts::Options options("pi_calculation", "Calculate pi using serial and parallel execution");
    options.add_options("custom", {
                                {"nPoints", "Number of points", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
                                {"nWorkers", "Number of workers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                            });

    auto cl_options = options.parse(argc, argv);
    uint n_points = cl_options["nPoints"].as<uint>();
    uint n_workers = cl_options["nWorkers"].as<uint>();
    std::cout << std::fixed;    
    std::cout << "Number of points : " << n_points << "\n";
    std::cout << "Number of workers : " << n_workers << "\n";

    piCalculation(n_points);
    
    return 0;
}
