// Canonical models with analytic expectations.
#include <algorithm>
#include <cmath>
#include <limits>

#include "testgen/testgen.hpp"

namespace cadcore::testgen {

std::vector<std::array<std::uint32_t, 3>> triangulate(const Face& f) {
    std::vector<std::array<std::uint32_t, 3>> out;
    for (const Prim& p : f.prims) {
        const auto& r = p.refs;
        if (p.kind == PrimKind::Triangles) {
            for (std::size_t k = 0; k + 2 < r.size(); k += 3) out.push_back({r[k], r[k + 1], r[k + 2]});
        } else if (p.kind == PrimKind::Strip) {
            for (std::size_t t = 0; t + 2 < r.size(); ++t)
                out.push_back(t % 2 == 0 ? std::array<std::uint32_t, 3>{r[t], r[t + 1], r[t + 2]} : std::array<std::uint32_t, 3>{r[t + 1], r[t], r[t + 2]});
        } else {
            for (std::size_t t = 1; t + 1 < r.size(); ++t) out.push_back({r[0], r[t], r[t + 1]});
        }
    }
    return out;
}

namespace {

constexpr double kPi = 3.14159265358979323846;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

Body cube_body(float s, bool colored) {
    Body b;
    b.name = "cube";
    for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x) b.points.push_back({x * s, y * s, z * s});
    auto id = [](int x, int y, int z) { return std::uint32_t(x + 2 * y + 4 * z); };
    // Each face is one 4-ref strip p0 p1 p3 p2 of a counter-clockwise (outward) quad p0..p3.
    struct Q { std::uint32_t p0, p1, p2, p3; std::array<float, 3> n; Color c; };
    const Q quads[6] = {
        {id(0, 0, 0), id(0, 1, 0), id(1, 1, 0), id(1, 0, 0), {0, 0, -1}, {220, 40, 40, 255}},
        {id(0, 0, 1), id(1, 0, 1), id(1, 1, 1), id(0, 1, 1), {0, 0, 1}, {40, 200, 40, 255}},
        {id(0, 0, 0), id(1, 0, 0), id(1, 0, 1), id(0, 0, 1), {0, -1, 0}, {40, 40, 220, 255}},
        {id(0, 1, 0), id(0, 1, 1), id(1, 1, 1), id(1, 1, 0), {0, 1, 0}, {220, 220, 40, 255}},
        {id(0, 0, 0), id(0, 0, 1), id(0, 1, 1), id(0, 1, 0), {-1, 0, 0}, {220, 40, 220, 128}},
        {id(1, 0, 0), id(1, 1, 0), id(1, 1, 1), id(1, 0, 1), {1, 0, 0}, {40, 220, 220, 255}},
    };
    for (const Q& q : quads) {
        Face f;
        f.prims.push_back({PrimKind::Strip, {q.p0, q.p1, q.p3, q.p2}});
        f.color = colored ? q.c : Color{210, 210, 255, 255};
        f.has_normal = true;
        f.normal = q.n;
        b.faces.push_back(f);
    }
    return b;
}

// Cylinder around (cx, cy), z in [0, h]: side strip, top centre fan, bottom ring fan.
// mixed_top: the first 3 top triangles are written as an explicit triangle list (tri + fan record).
Body cylinder_body(int n, float r, float h, float cx, float cy, bool caps, bool mixed_top) {
    Body b;
    b.name = "cylinder" + std::to_string(n);
    b.solid = caps;
    for (int i = 0; i < n; ++i) {
        const double a = 2 * kPi * i / n;
        b.points.push_back({cx + float(r * std::cos(a)), cy + float(r * std::sin(a)), 0});
    }
    for (int i = 0; i < n; ++i) {
        const double a = 2 * kPi * i / n;
        b.points.push_back({cx + float(r * std::cos(a)), cy + float(r * std::sin(a)), h});
    }
    auto bot = [&](int i) { return std::uint32_t(i % n); };
    auto top = [&](int i) { return std::uint32_t(n + i % n); };
    Face side;
    Prim strip{PrimKind::Strip, {}};
    for (int i = 0; i <= n; ++i) {
        strip.refs.push_back(top(i));
        strip.refs.push_back(bot(i));
    }
    side.prims.push_back(strip);
    side.color = caps ? Color{180, 180, 190, 255} : Color{243, 254, 177, 255};
    b.faces.push_back(side);
    if (!caps) return b;
    b.points.push_back({cx, cy, h});
    const auto centre = std::uint32_t(b.points.size() - 1);
    Face tf;
    tf.has_normal = true;
    tf.normal = {0, 0, 1};
    tf.color = {200, 120, 60, 255};
    int start = 0;
    if (mixed_top) {
        Prim tris{PrimKind::Triangles, {}};
        for (int i = 0; i < 3; ++i) tris.refs.insert(tris.refs.end(), {centre, top(i), top(i + 1)});
        tf.prims.push_back(tris);
        start = 3;
    }
    Prim fan{PrimKind::Fan, {centre}};
    for (int i = start; i <= n; ++i) fan.refs.push_back(top(i));
    tf.prims.push_back(fan);
    b.faces.push_back(tf);
    Face bf;
    bf.has_normal = true;
    bf.normal = {0, 0, -1};
    bf.color = {60, 120, 200, 255};
    Prim bfan{PrimKind::Fan, {bot(0)}};
    for (int i = n - 1; i >= 1; --i) bfan.refs.push_back(bot(i));
    bf.prims.push_back(bfan);
    b.faces.push_back(bf);
    return b;
}

double polygon_area(int n, double r) { return 0.5 * n * r * r * std::sin(2 * kPi / n); }
double side_area(int n, double r, double h) { return n * 2 * r * std::sin(kPi / n) * h; }

}  // namespace

std::vector<Model> canonical_models() {
    std::vector<Model> out;
    out.push_back({"cube", {cube_body(10, true)}, {}, {12, 600, 1000, {0, 0, 0}, {10, 10, 10}}});
    {
        const int n = 32;
        const double r = 5, h = 10;
        out.push_back({"cylinder32", {cylinder_body(n, float(r), float(h), 0, 0, true, false)}, {},
                       {std::size_t(4 * n - 2), 2 * polygon_area(n, r) + side_area(n, r, h), polygon_area(n, r) * h, {-r, -r, 0}, {r, r, h}}});
    }
    {
        const int n = 300;  // > 255 strip/fan lengths (escaped counts) and > 256 elements (u16 refs)
        const double r = 50, h = 20;
        out.push_back({"cylinder300", {cylinder_body(n, float(r), float(h), 0, 0, true, true)}, {},
                       {std::size_t(4 * n - 2), 2 * polygon_area(n, r) + side_area(n, r, h), polygon_area(n, r) * h, {-r, -r, 0}, {r, r, h}}});
    }
    {
        const int n = 35000;  // 70001 elements: needs 24-bit references in CATPart
        const double r = 100, h = 5;
        out.push_back({"dense_cylinder35k", {cylinder_body(n, float(r), float(h), 0, 0, true, false)}, {},
                       {std::size_t(4 * n - 2), 2 * polygon_area(n, r) + side_area(n, r, h), polygon_area(n, r) * h, {-r, -r, 0}, {r, r, h}}});
    }
    {
        const int n = 32;
        Model m{"solid_and_surface", {cube_body(10, false), cylinder_body(n, 5, 10, 20, 5, false, false)}, {}, {}};
        m.expected = {std::size_t(12 + 2 * n), 600 + side_area(n, 5, 10), kNaN, {0, 0, 0}, {25, 10, 10}};
        out.push_back(m);
    }
    {
        Model m{"assembly3", {cube_body(10, true)}, {}, {36, 1800, 3000, {-10, 0, 0}, {30, 40, 10}}};
        Matrix4 shift = kIdentity;
        shift[12] = 20;
        Matrix4 rot = {0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 30, 0, 1};  // 90 deg about z, then +30 in y
        m.instances = {{0, kIdentity}, {0, shift}, {0, rot}};
        out.push_back(m);
    }
    return out;
}

Model grid_model(int k, bool face_per_quad) {
    Body b;
    b.name = "grid" + std::to_string(k);
    b.solid = false;
    const auto id = [k](int x, int y) { return std::uint32_t(y * (k + 1) + x); };
    double zmin = 0, zmax = 0;
    for (int y = 0; y <= k; ++y)
        for (int x = 0; x <= k; ++x) {
            const float z = float(((x * 7 + y * 13) % 11)) * 0.01f;  // varied heights: no repeated rows
            zmin = std::min<double>(zmin, z);
            zmax = std::max<double>(zmax, z);
            b.points.push_back({float(x), float(y), z});
        }
    const Color colors[2] = {{200, 80, 60, 255}, {60, 120, 200, 255}};
    for (int y = 0; y < k; ++y) {
        if (face_per_quad) {
            for (int x = 0; x < k; ++x) {
                Face f;
                f.prims.push_back({PrimKind::Strip, {id(x, y), id(x + 1, y), id(x, y + 1), id(x + 1, y + 1)}});
                f.color = colors[(x + y) & 1];
                b.faces.push_back(std::move(f));
            }
        } else {
            Face f;
            Prim strip{PrimKind::Strip, {}};
            for (int x = 0; x <= k; ++x) {
                strip.refs.push_back(id(x, y));
                strip.refs.push_back(id(x, y + 1));
            }
            f.prims.push_back(std::move(strip));
            f.color = colors[y & 1];
            b.faces.push_back(std::move(f));
        }
    }
    Model m{(face_per_quad ? "gridfaces" : "grid") + std::to_string(k), {b}, {}, {}};
    m.expected = {std::size_t(2) * std::size_t(k) * std::size_t(k), kNaN, kNaN, {0, 0, zmin}, {double(k), double(k), zmax}};
    return m;
}

Model instances_model(int count) {
    Model m{"instances" + std::to_string(count), {cube_body(1, true)}, {}, {}};
    for (int i = 0; i < count; ++i) {
        Matrix4 t = kIdentity;
        t[12] = 2.0 * (i % 100);
        t[13] = 2.0 * (i / 100);
        m.instances.push_back({0, t});
    }
    m.expected = {std::size_t(12) * std::size_t(count), 6.0 * count, 1.0 * count, {0, 0, 0},
                  {2.0 * std::min(count - 1, 99) + 1, 2.0 * ((count - 1) / 100) + 1, 1}};
    return m;
}

}  // namespace cadcore::testgen
