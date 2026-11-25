// Minimal self-contained test harness (no third-party code).
//
//   TEST_CASE("name") { ... }        registers a test
//   CHECK(expr)                       records a failure and continues
//   REQUIRE(expr)                     records a failure and ends the test case
//   CHECK_NOTHROW(expr), REQUIRE_NOTHROW(expr)
//   CAPTURE(value)                    adds "name := value" to failure reports in the current scope
#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <cstdio>
#include <exception>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace harness {

struct Case {
    const char* name;
    void (*fn)();
};

inline std::vector<Case>& registry() {
    static std::vector<Case> r;
    return r;
}
// Capture context of the running case (cases run on several threads: thread-local).
inline std::vector<std::function<std::string()>>& context() {
    thread_local std::vector<std::function<std::string()>> c;
    return c;
}
struct Stats {
    std::atomic<long> assertions{0}, failures{0};
};
inline long& case_failures() {  // failures of the case running on this thread
    thread_local long n = 0;
    return n;
}
inline std::mutex& output() {
    static std::mutex m;
    return m;
}
inline Stats& stats() {
    static Stats s;
    return s;
}
struct RequireFailed {};

inline void report(const char* file, int line, const char* expr, const char* what, const char* threw = nullptr) {
    ++stats().failures;
    ++case_failures();
    const std::lock_guard<std::mutex> lock(output());
    std::fprintf(stderr, "%s:%d: FAILED %s( %s )\n", file, line, what, expr);
    if (threw) std::fprintf(stderr, "    threw: %s\n", threw);
    for (auto& c : context()) std::fprintf(stderr, "    with %s\n", c().c_str());
}

struct Registrar {
    Registrar(const char* name, void (*fn)()) { registry().push_back({name, fn}); }
};

struct Capture {
    template <class T>
    Capture(const char* name, const T& value) {
        context().push_back([name, &value] {
            std::ostringstream o;
            o << name << " := " << value;
            return o.str();
        });
    }
    ~Capture() { context().pop_back(); }
};

// Runs every case whose name contains `filter` (all when null) on CADCORE_TEST_JOBS threads
// (default 4; cases must not share mutable state); lists cases slower than 0.5 s.
inline int run_all(const char* filter = nullptr) {
    std::vector<const Case*> cases;
    for (const Case& c : registry())
        if (!filter || std::string(c.name).find(filter) != std::string::npos) cases.push_back(&c);
    std::atomic<std::size_t> next{0};
    std::atomic<int> failed_cases{0};
    auto worker = [&] {
        for (std::size_t i; (i = next++) < cases.size();) {
            const Case& c = *cases[i];
            case_failures() = 0;
            const auto t0 = std::chrono::steady_clock::now();
            try {
                c.fn();
            } catch (const RequireFailed&) {
            } catch (const std::exception& e) {
                ++stats().failures;
                ++case_failures();
                const std::lock_guard<std::mutex> lock(output());
                std::fprintf(stderr, "test case '%s' threw: %s\n", c.name, e.what());
            }
            context().clear();
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            const std::lock_guard<std::mutex> lock(output());
            if (sec > 0.5) std::fprintf(stderr, "slow test case (%.1f s): %s\n", sec, c.name);
            if (case_failures()) {
                ++failed_cases;
                std::fprintf(stderr, "  in test case: %s\n", c.name);
            }
        }
    };
    const char* env = std::getenv("CADCORE_TEST_JOBS");
    const int jobs = std::max(1, env ? std::atoi(env) : 4);
    std::vector<std::thread> pool;
    for (int j = 1; j < jobs; ++j) pool.emplace_back(worker);
    worker();
    for (auto& t : pool) t.join();
    std::printf("test cases: %zu | %d failed\nassertions: %ld | %ld failed\nStatus: %s\n", cases.size(), failed_cases.load(), stats().assertions.load(),
                stats().failures.load(), stats().failures ? "FAILURE" : "SUCCESS");
    return stats().failures ? 1 : 0;
}

}  // namespace harness

#define HARNESS_CAT2(a, b) a##b
#define HARNESS_CAT(a, b) HARNESS_CAT2(a, b)
#define TEST_CASE(name)                                                                              \
    static void HARNESS_CAT(harness_case_, __LINE__)();                                              \
    static const harness::Registrar HARNESS_CAT(harness_reg_, __LINE__)(name, &HARNESS_CAT(harness_case_, __LINE__)); \
    static void HARNESS_CAT(harness_case_, __LINE__)()
#define CHECK(expr)                                                                      \
    do {                                                                                 \
        ++harness::stats().assertions;                                                   \
        if (!(expr)) harness::report(__FILE__, __LINE__, #expr, "CHECK");                \
    } while (0)
#define REQUIRE(expr)                                                                    \
    do {                                                                                 \
        ++harness::stats().assertions;                                                   \
        if (!(expr)) {                                                                   \
            harness::report(__FILE__, __LINE__, #expr, "REQUIRE");                       \
            throw harness::RequireFailed{};                                              \
        }                                                                                \
    } while (0)
#define CHECK_NOTHROW(expr)                                                              \
    do {                                                                                 \
        ++harness::stats().assertions;                                                   \
        try { expr; } catch (const std::exception& e_) {                                 \
            harness::report(__FILE__, __LINE__, #expr, "CHECK_NOTHROW", e_.what());     \
        }                                                                                \
    } while (0)
#define REQUIRE_NOTHROW(expr)                                                            \
    do {                                                                                 \
        ++harness::stats().assertions;                                                   \
        try { expr; } catch (const std::exception& e_) {                                 \
            harness::report(__FILE__, __LINE__, #expr, "REQUIRE_NOTHROW", e_.what());   \
            throw harness::RequireFailed{};                                              \
        }                                                                                \
    } while (0)
#define CAPTURE(x) const harness::Capture HARNESS_CAT(harness_capture_, __LINE__)(#x, x)
