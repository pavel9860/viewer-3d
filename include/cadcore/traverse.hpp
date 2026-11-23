// Walking a decoded scene in world space.
#pragma once

#include <array>
#include <functional>

#include "cadcore/scene.hpp"

namespace cadcore {

Matrix4 multiply(const Matrix4& a, const Matrix4& b);

// The point xyz under `m`.
std::array<double, 3> apply(const Matrix4& m, const float* xyz);

// Calls `visit(node, world)` for every node reached from the root, instances expanded
// (a node shared by several parents is visited once per path). A scene without nodes is
// treated as one root holding every mesh and polyline. Throws DecodeError(Malformed) when
// the tree is deeper than 512 (a cycle).
void for_each_node(const Scene& scene, const std::function<void(const Node&, const Matrix4&)>& visit);

}  // namespace cadcore
