#include "cadcore/compare.hpp"

#include <cmath>
#include <cstdint>
#include <unordered_map>

#include "cadcore/error.hpp"
#include "cadcore/stats.hpp"

namespace cadcore {
namespace {

template <class F>
void for_each_instance(const Scene& scene, F&& f) {
    for_each_node(scene, [&](const Node& n, const Matrix4& world) {
        for (auto mi : n.meshes) f(scene.meshes.at(mi), world);
    });
}

struct Grid {
    double cell;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells;
    static std::uint64_t key(std::int64_t x, std::int64_t y, std::int64_t z) {
        return (std::uint64_t(x) * 73856093u) ^ (std::uint64_t(y) * 19349663u) ^ (std::uint64_t(z) * 83492791u);
    }
    std::array<std::int64_t, 3> idx(const double* p) const {
        return {std::int64_t(std::floor(p[0] / cell)), std::int64_t(std::floor(p[1] / cell)), std::int64_t(std::floor(p[2] / cell))};
    }
};

double dist(const double* a, const double* b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template <std::size_t N, class Score>
MatchReport match(const std::vector<std::array<double, N>>& a, const std::vector<std::array<double, N>>& b, double tol, Score score) {
    MatchReport r;
    r.a = a.size();
    r.b = b.size();
    const double cell = tol > 0 ? tol * 2 : 1e-9;
    Grid g{cell, {}};
    auto centroid = [](const std::array<double, N>& e) {
        std::array<double, 3> c{0, 0, 0};
        for (std::size_t i = 0; i < N; ++i) c[i % 3] += e[i];
        for (auto& v : c) v /= double(N / 3);
        return c;
    };
    for (std::size_t i = 0; i < b.size(); ++i) {
        const auto c = centroid(b[i]);
        const auto k = g.idx(c.data());
        g.cells[Grid::key(k[0], k[1], k[2])].push_back(i);
    }
    std::vector<bool> used(b.size(), false);
    for (const auto& e : a) {
        const auto c = centroid(e);
        const auto k = g.idx(c.data());
        double best = -1;
        std::size_t best_j = 0;
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    auto it = g.cells.find(Grid::key(k[0] + dx, k[1] + dy, k[2] + dz));
                    if (it == g.cells.end()) continue;
                    for (const std::size_t j : it->second) {
                        if (used[j]) continue;
                        const double s = score(e, b[j]);
                        if (s <= tol && (best < 0 || s < best)) { best = s; best_j = j; }
                    }
                }
        if (best < 0) { ++r.unmatched_a; continue; }
        used[best_j] = true;
        r.max_deviation = std::max(r.max_deviation, best);
    }
    for (bool u : used) r.unmatched_b += u ? 0 : 1;
    return r;
}

}  // namespace

std::vector<WorldTriangle> world_triangles(const Scene& scene) {
    std::vector<WorldTriangle> out;
    for_each_instance(scene, [&](const Mesh& m, const Matrix4& t) {
        if (m.topology != Topology::Triangles) return;
        for (std::size_t k = 0; k + 2 < m.indices.size(); k += 3) {
            WorldTriangle tri;
            for (int c = 0; c < 3; ++c) {
                const auto p = apply(t, &m.positions[std::size_t(m.indices[k + std::size_t(c)]) * 3]);
                for (int i = 0; i < 3; ++i) tri[std::size_t(c * 3 + i)] = p[std::size_t(i)];
            }
            out.push_back(tri);
        }
    });
    return out;
}

std::vector<std::array<double, 3>> world_points(const Scene& scene) {
    std::vector<std::array<double, 3>> out;
    for_each_instance(scene, [&](const Mesh& m, const Matrix4& t) {
        if (m.topology != Topology::Points) return;
        for (std::size_t v = 0; v < m.vertex_count(); ++v) out.push_back(apply(t, &m.positions[v * 3]));
    });
    return out;
}

MatchReport match_triangles(const std::vector<WorldTriangle>& a, const std::vector<WorldTriangle>& b, double tolerance) {
    return match(a, b, tolerance, [](const WorldTriangle& x, const WorldTriangle& y) {
        double best = -1;
        for (int rot = 0; rot < 3; ++rot) {
            double worst = 0;
            for (int c = 0; c < 3; ++c) worst = std::max(worst, dist(&x[std::size_t(c * 3)], &y[std::size_t(((c + rot) % 3) * 3)]));
            if (best < 0 || worst < best) best = worst;
        }
        return best;
    });
}

MatchReport match_points(const std::vector<std::array<double, 3>>& a, const std::vector<std::array<double, 3>>& b, double tolerance) {
    return match(a, b, tolerance, [](const std::array<double, 3>& x, const std::array<double, 3>& y) { return dist(x.data(), y.data()); });
}

}  // namespace cadcore
