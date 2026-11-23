// Minimal self-contained test harness (no third-party code).
//
//   TEST_CASE("name") { ... }        registers a test
//   CHECK(expr)                       records a failure and continues
//   REQUIRE(expr)                     records a failure and ends the test case
//   CHECK_NOTHROW(expr), REQUIRE_NOTHROW(expr)
//   CAPTURE(value)                    adds "name := value" to failure reports in the current scope
#pragma once

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
inline std::vector<std::function<std::string()>>& context() {
    static std::vector<std::function<std::string()>> c;
    return c;
}
struct Stats {
    long assertions = 0, failures = 0;
};
inline Stats& stats() {
    static Stats s;
    return s;
}
struct RequireFailed {};

inline void report(const char* file, int line, const char* expr, const char* what) {
    ++stats().failures;
    std::fprintf(stderr, "%s:%d: FAILED %s( %s )\n", file, line, what, expr);
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

inline int run_all() {
    int failed_cases = 0;
    for (const Case& c : registry()) {
        const long before = stats().failures;
        try {
            c.fn();
        } catch (const RequireFailed&) {
        } catch (const std::exception& e) {
            ++stats().failures;
            std::fprintf(stderr, "test case '%s' threw: %s\n", c.name, e.what());
        }
        context().clear();
        if (stats().failures != before) {
            ++failed_cases;
            std::fprintf(stderr, "  in test case: %s\n", c.name);
        }
    }
    std::printf("test cases: %zu | %d failed\nassertions: %ld | %ld failed\nStatus: %s\n", registry().size(), failed_cases, stats().assertions,
                stats().failures, stats().failures ? "FAILURE" : "SUCCESS");
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
            harness::report(__FILE__, __LINE__, #expr, "CHECK_NOTHROW");                 \
            std::fprintf(stderr, "    threw: %s\n", e_.what());                          \
        }                                                                                \
    } while (0)
#define REQUIRE_NOTHROW(expr)                                                            \
    do {                                                                                 \
        ++harness::stats().assertions;                                                   \
        try { expr; } catch (const std::exception& e_) {                                 \
            harness::report(__FILE__, __LINE__, #expr, "REQUIRE_NOTHROW");               \
            std::fprintf(stderr, "    threw: %s\n", e_.what());                          \
            throw harness::RequireFailed{};                                              \
        }                                                                                \
    } while (0)
#define CAPTURE(x) const harness::Capture HARNESS_CAT(harness_capture_, __LINE__)(#x, x)
