// Parallel execution with a bounded number of threads (cadcore uses at most 6).
//
// Tasks are collected first and run by run_all(); results are written by the tasks into
// caller-owned slots, so output order never depends on scheduling. Builds without threads
// (WebAssembly without pthreads) run the tasks in order on the calling thread.
#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

namespace cadcore {

constexpr unsigned kMaxThreads = 6;

// Per-piece output of a parallel loop, one cache line each: neighbouring pieces written by
// different threads would otherwise share cache lines (false sharing).
template <class T> struct alignas(64) Slot {
    T value;
};

// Threads to use for a request of `requested` (0 = all cores), capped at kMaxThreads.
unsigned thread_count(unsigned requested);

// Runs every task on up to `threads` threads and waits. The first exception thrown by a task
// (in task order) is rethrown after all tasks have finished.
void run_all(std::vector<std::function<void()>>& tasks, unsigned threads);

// Splits [0, n) into contiguous ranges of at least `grain` items and calls fn(begin, end) for each,
// on up to `threads` threads.
void parallel_for(std::size_t n, std::size_t grain, unsigned threads, const std::function<void(std::size_t, std::size_t)>& fn);

// Appends the vectors piece(0) .. piece(n-1) to dst: one resize, then the copies in parallel.
template <class T, class Piece> void append_all(std::vector<T>& dst, std::size_t n, Piece&& piece, unsigned threads) {
    std::vector<std::size_t> at(n + 1, dst.size());
    for (std::size_t k = 0; k < n; ++k) at[k + 1] = at[k] + piece(k).size();
    if (at[n] == dst.size()) return;
    dst.resize(at[n]);
    parallel_for(n, 1, threads, [&](std::size_t b, std::size_t e) {
        for (std::size_t k = b; k < e; ++k) std::copy(piece(k).begin(), piece(k).end(), dst.begin() + std::ptrdiff_t(at[k]));
    });
}

}  // namespace cadcore
