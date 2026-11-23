#include "core/parallel.hpp"

#include <algorithm>
#include <atomic>
#include <exception>

#if !defined(__EMSCRIPTEN__) || defined(__EMSCRIPTEN_PTHREADS__)
#define CADCORE_THREADS 1
#include <condition_variable>
#include <deque>
#include <mutex>
#include <set>
#include <thread>
#endif

namespace cadcore {

unsigned thread_count(unsigned requested) {
#ifdef CADCORE_THREADS
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    return std::min(kMaxThreads, requested ? requested : hw);
#else
    (void)requested;
    return 1;
#endif
}

#ifdef CADCORE_THREADS
namespace {

// One batch of tasks; whoever is free takes the next index.
struct Batch {
    std::vector<std::function<void()>>* tasks;
    std::vector<std::exception_ptr>* errors;
    std::atomic<std::size_t> next{0}, done{0};
    unsigned helpers;  // pool workers allowed to join (besides the caller)
    bool run_one() {
        const std::size_t i = next.fetch_add(1);
        if (i >= tasks->size()) return false;
        try { (*tasks)[i](); } catch (...) { (*errors)[i] = std::current_exception(); }
        done.fetch_add(1);
        return true;
    }
};

// kMaxThreads - 1 workers, started on first use and kept for the life of the process.
class Pool {
public:
    // Never destroyed: the detached workers wait on its condition variable until the process ends
    // (destroying a condition variable with waiters blocks).
    static Pool& get() { static Pool* p = new Pool; return *p; }
    void run(Batch& b) {
        {
            std::lock_guard<std::mutex> l(m_);
            for (unsigned k = 0; k < b.helpers; ++k) queue_.push_back(&b);
        }
        cv_.notify_all();
        while (b.run_one()) {}
        {  // withdraw helper slots no worker picked up
            std::lock_guard<std::mutex> l(m_);
            queue_.erase(std::remove(queue_.begin(), queue_.end(), &b), queue_.end());
        }
        std::unique_lock<std::mutex> l(m_);
        done_cv_.wait(l, [&] { return b.done.load() == b.tasks->size() && active_.count(&b) == 0; });
    }

private:
    Pool() {
        for (unsigned t = 1; t < kMaxThreads; ++t) std::thread([this] { work(); }).detach();
    }
    void work() {
        std::unique_lock<std::mutex> l(m_);
        while (true) {
            cv_.wait(l, [&] { return !queue_.empty(); });
            Batch* b = queue_.front();
            queue_.pop_front();
            active_.insert(b);
            l.unlock();
            while (b->run_one()) {}
            l.lock();
            active_.erase(active_.find(b));
            done_cv_.notify_all();
        }
    }
    std::mutex m_;
    std::condition_variable cv_, done_cv_;
    std::deque<Batch*> queue_;
    std::multiset<Batch*> active_;
};

}  // namespace
#endif

void run_all(std::vector<std::function<void()>>& tasks, unsigned threads) {
    std::vector<std::exception_ptr> errors(tasks.size());
    const unsigned n = unsigned(std::min<std::size_t>(thread_count(threads), tasks.size()));
#ifdef CADCORE_THREADS
    if (n > 1) {
        Batch b;
        b.tasks = &tasks;
        b.errors = &errors;
        b.helpers = n - 1;
        Pool::get().run(b);
    } else
#endif
    {
        (void)n;
        for (std::size_t i = 0; i < tasks.size(); ++i) {
            try { tasks[i](); } catch (...) { errors[i] = std::current_exception(); }
        }
    }
    for (auto& e : errors)
        if (e) std::rethrow_exception(e);
}

void parallel_for(std::size_t n, std::size_t grain, unsigned threads, const std::function<void(std::size_t, std::size_t)>& fn) {
    const unsigned t = thread_count(threads);
    const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>(t, n / std::max<std::size_t>(1, grain)));
    if (parts == 1) {
        if (n) fn(0, n);
        return;
    }
    std::vector<std::function<void()>> tasks;
    for (std::size_t k = 0; k < parts; ++k) {
        const std::size_t b = n * k / parts, e = n * (k + 1) / parts;
        tasks.push_back([&fn, b, e] { fn(b, e); });
    }
    run_all(tasks, t);
}

}  // namespace cadcore
