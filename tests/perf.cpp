// Speed and complexity tests (separate binary: cadcore_perf).
//
// Complexity: every format variant decodes generated models whose size doubles three times
// (long strips, one face per quad, many instances). The decode time is the minimum of
// repeated runs; the growth exponent of time against input size (least squares on log-log)
// must stay near 1 (linear). The exponent is a ratio, so it does not depend on machine speed.
//
// Speed: the throughput at the largest size must stay above a floor set at about a quarter
// of the speed measured when the floors were set (2026-11-22, release build). Floors catch
// gross regressions only; `cadcli benchsuite --check` compares against a per-machine
// baseline for finer tracking. Sanitizer builds skip the floors (not the exponents).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "cadcore/decode.hpp"
#include "cadcore/stats.hpp"
#include "harness.hpp"
#include "testgen/testgen.hpp"

using namespace cadcore;
using namespace cadcore::testgen;

#if defined(__SANITIZE_ADDRESS__)
constexpr bool kSanitized = true;
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
constexpr bool kSanitized = true;
#else
constexpr bool kSanitized = false;
#endif
#else
constexpr bool kSanitized = false;
#endif

namespace {

constexpr double kMaxExponent = 1.3;  // linear ~1.0, n log n ~1.05-1.1, quadratic ~2

struct Timing {
    double seconds = 0;
    std::size_t bytes = 0, triangles = 0;
};

Timing time_decode(const Files& files) {
    DecodeOptions opt;
    opt.load_resource = [&](std::string_view path) -> std::optional<std::vector<std::uint8_t>> {
        for (const auto& [name, bytes] : files)
            if (name == path) return bytes;
        return std::nullopt;
    };
    const auto& main = files.front();
    Timing t;
    for (const auto& f : files) t.bytes += f.second.size();
    const ByteView v(main.second.data(), main.second.size());
    t.triangles = measure(decode(v, main.first, opt)).triangles;  // warm-up + result check
    double best = 1e30, total = 0;
    for (int rep = 0; rep < 25 && (rep < 5 || total < 0.15); ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        const Scene s = decode(v, main.first, opt);
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (s.meshes.empty()) throw std::runtime_error("empty decode");
        best = std::min(best, dt);
        total += dt;
    }
    t.seconds = best;
    return t;
}

double exponent(const std::vector<Timing>& ts) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double n = double(ts.size());
    for (const auto& t : ts) {
        const double x = std::log(double(t.bytes)), y = std::log(t.seconds);
        sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// Minimum input MB/s at the largest size: a quarter of the slowest family measured on
// 2026-11-22 (release, -O2). 3MF is measured against its compressed size (~7x smaller).
double floor_mbps(const std::string& id) {
    static const std::map<std::string, double> floors = {
        {"stl-binary", 800}, {"stl-ascii", 90}, {"obj", 35}, {"ply-ascii", 30}, {"ply-le", 45}, {"ply-be", 70},
        {"off-ascii", 55}, {"off-binary", 240}, {"gltf", 6}, {"glb", 5}, {"3mf", 1.3}, {"catpart", 25},
    };
    for (const auto& [prefix, mbps] : floors)
        if (id.rfind(prefix, 0) == 0) return mbps;
    return 5;
}

struct Family {
    const char* name;
    std::vector<Model> (*models)();
    bool instanced;
};

std::vector<Model> strips() { return {grid_model(128, false), grid_model(181, false), grid_model(256, false), grid_model(362, false)}; }
std::vector<Model> quads() { return {grid_model(64, true), grid_model(90, true), grid_model(128, true), grid_model(181, true)}; }
std::vector<Model> nodes() { return {instances_model(1000), instances_model(2000), instances_model(4000), instances_model(8000)}; }

void run_family(const Family& fam) {
    const auto models = fam.models();
    for (const Variant& v : variants()) {
        if (fam.instanced && v.id.rfind("catpart", 0) == 0) continue;  // CATPart has no instances
        CAPTURE(fam.name);
        CAPTURE(v.id);
        std::vector<Files> files;
        for (const Model& m : models) files.push_back(v.write(m));
        // A load spike on the machine can distort one series: measure twice at most and keep
        // the better series (both checks still apply to it).
        std::vector<Timing> ts;
        double e = 0, mbps = 0;
        for (int attempt = 0; attempt < 2; ++attempt) {
            std::vector<Timing> run;
            for (std::size_t i = 0; i < models.size(); ++i) {
                Timing t;
                REQUIRE_NOTHROW(t = time_decode(files[i]));
                CHECK(t.triangles == models[i].expected.triangles);
                run.push_back(t);
            }
            const double re = exponent(run), rm = double(run.back().bytes) / 1e6 / run.back().seconds;
            if (attempt == 0 || re < e) { ts = run; e = re; mbps = rm; }
            if (e <= kMaxExponent && (kSanitized || mbps >= floor_mbps(v.id))) break;
        }
        const Timing& big = ts.back();
        CAPTURE(e);
        CAPTURE(mbps);
        std::fflush(stderr);
        std::printf("  %-8s %-20s exponent %.2f   %8.2f ms  %8.1f MB/s  %9zu triangles\n", fam.name, v.id.c_str(), e, big.seconds * 1e3, mbps,
                    big.triangles);
        std::fflush(stdout);
        CHECK(e <= kMaxExponent);
        if (!kSanitized) CHECK(mbps >= floor_mbps(v.id));
    }
}

}  // namespace

TEST_CASE("decode time grows linearly: long strips") { run_family({"strips", strips, false}); }
TEST_CASE("decode time grows linearly: one face per quad") { run_family({"quads", quads, false}); }
TEST_CASE("decode time grows linearly: many instances") { run_family({"nodes", nodes, true}); }
