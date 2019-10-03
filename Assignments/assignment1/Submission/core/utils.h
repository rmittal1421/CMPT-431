#ifndef UTILS_H
#define UTILS_H

#include "cxxopts.h"
#include "get_time.h"
#include <iostream>
#include <mutex>
#include <atomic>
#include <condition_variable>

#define intV int32_t
#define uintV int32_t

#define intE int32_t
#define uintE int32_t

#define DEFAULT_NUMBER_OF_WORKERS "1"
#define DEFAULT_MAX_ITER "10"
#define TIME_PRECISION 5
#define VAL_PRECISION 14
// #define ADDITIONAL_TIMER_LOGS 0

struct CustomBarrier
{
    int num_of_workers_;
    int current_waiting_;
    std::mutex my_mutex_;
    std::condition_variable my_cv_;

    CustomBarrier(int t_num_of_workers) : num_of_workers_(t_num_of_workers), current_waiting_(0) {}

    void wait()
    {
        std::unique_lock<std::mutex> u_lock(my_mutex_);
        current_waiting_++;
        if (current_waiting_ == num_of_workers_)
        {
            current_waiting_ = 0;
            // unlock and send signal to wake up
            u_lock.unlock();
            my_cv_.notify_all();
            return;
        }
        // unlock and continue. It will wait on the condition variable
        my_cv_.wait(u_lock);
        //  Condition has been reached. return
    }
};

#endif