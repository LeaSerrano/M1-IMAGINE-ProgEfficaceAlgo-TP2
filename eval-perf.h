#ifndef __EVAL_PERF_H
#define __EVAL_PERF_H

#include <chrono>
#include <x86intrin.h>
#include <vector>

uint64_t rdtsc() {
    unsigned int lo , hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

class EvalPerf {
    private:
    std::chrono::time_point<std::chrono::high_resolution_clock> t_start, t_end;
    uint64_t c_start, c_end;

    public:
    void start() {
        t_start = std::chrono::high_resolution_clock::now();
        c_start = rdtsc();
    }

    void stop() {
        t_end = std::chrono::high_resolution_clock::now();
        c_end = rdtsc();
    }

    void clear() {
        t_start = t_end;
    }

    auto millisecond() const {
        return std::chrono::duration<double, std::milli>(t_end-t_start).count();
    }

    auto seconds() const { 
        return std::chrono::duration<double>(t_end-t_start).count();
    }

    uint64_t cycles() const {
        return c_end - c_start;
    }

    double CPI(size_t N) const { 
        return double(c_end-c_start)/N;
    }

    double IPC(size_t N) const { 
        return N/double(c_end-c_start);
    }
  
    double Gflops(size_t N) const { 
        return double(N)/seconds();
    }
    
};

#endif