// Geometric comparison of two decoded scenes, independent of vertex order, face order and
// the starting corner of each triangle (orientation is preserved).
#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "cadcore/scene.hpp"

namespace cadcore {

using WorldTriangle = std::array<double, 9>;  // three xyz corners

// Every triangle of the scene in world space (instances expanded).
std::vector<WorldTriangle> world_triangles(const Scene& scene);
// Every point-cloud vertex of the scene in world space.
std::vector<std::array<double, 3>> world_points(const Scene& scene);

struct MatchReport {
    std::size_t a = 0, b = 0;          // element counts
    std::size_t unmatched_a = 0, unmatched_b = 0;
    double max_deviation = 0;          // largest corner distance among matched pairs
    bool identical() const { return unmatched_a == 0 && unmatched_b == 0; }
};

// One-to-one matching: a triangle of `a` matches a triangle of `b` when each corresponding
// corner (under some cyclic rotation) lies within `tolerance`.
MatchReport match_triangles(const std::vector<WorldTriangle>& a, const std::vector<WorldTriangle>& b, double tolerance);
MatchReport match_points(const std::vector<std::array<double, 3>>& a, const std::vector<std::array<double, 3>>& b, double tolerance);

}  // namespace cadcore
