// World-space measurements of a decoded scene (instances expanded through the node tree).
#pragma once

#include <array>
#include <cstddef>

#include "cadcore/scene.hpp"
#include "cadcore/traverse.hpp"

namespace cadcore {

struct SceneStats {
    std::size_t instances = 0;   // mesh instances reached from the root
    std::size_t triangles = 0;   // world triangles (instances expanded)
    double area = 0;             // sum of triangle areas
    double signed_volume = 0;    // divergence theorem; meaningful for closed, consistently wound meshes
    std::array<double, 3> min{0, 0, 0}, max{0, 0, 0};
    bool empty = true;
};

SceneStats measure(const Scene& scene);


}  // namespace cadcore
