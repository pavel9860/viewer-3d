#include "cadcore/stats.hpp"

#include <algorithm>
#include <cmath>

#include "cadcore/error.hpp"

namespace cadcore {

static void accumulate(const Mesh& m, const Matrix4& t, SceneStats& st) {
    const std::size_t nv = m.vertex_count();
    std::vector<double> w(nv * 3);
    for (std::size_t i = 0; i < nv; ++i) {
        const double x = m.positions[i * 3], y = m.positions[i * 3 + 1], z = m.positions[i * 3 + 2];
        for (int r = 0; r < 3; ++r) w[i * 3 + r] = t[r] * x + t[4 + r] * y + t[8 + r] * z + t[12 + r];
    }
    for (std::size_t i = 0; i < nv; ++i)
        for (int r = 0; r < 3; ++r) {
            if (st.empty) st.min[r] = st.max[r] = w[i * 3 + r];
            st.min[r] = std::min(st.min[r], w[i * 3 + r]);
            st.max[r] = std::max(st.max[r], w[i * 3 + r]);
            if (r == 2) st.empty = false;
        }
    for (std::size_t k = 0; k + 2 < m.indices.size(); k += 3) {
        const double* a = &w[m.indices[k] * 3];
        const double* b = &w[m.indices[k + 1] * 3];
        const double* c = &w[m.indices[k + 2] * 3];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        st.area += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        st.signed_volume += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0;
    }
    st.triangles += m.triangle_count();
    st.instances += 1;
}

SceneStats measure(const Scene& scene) {
    SceneStats st;
    for_each_node(scene, [&](const Node& n, const Matrix4& world) {
        for (auto mi : n.meshes) accumulate(scene.meshes.at(mi), world, st);
    });
    return st;
}

}  // namespace cadcore
