#include "formats/dxf/clip.hpp"

#include <algorithm>
#include <cmath>

namespace cadcore::formats::autocad {
namespace {

bool inside(const Polyline3& l, double x, double y) {
    bool in = false;
    for (std::size_t i = 0, j = l.size() - 1; i < l.size(); j = i++)
        if ((l[i].y > y) != (l[j].y > y) && x < (l[j].x - l[i].x) * (y - l[i].y) / (l[j].y - l[i].y) + l[i].x) in = !in;
    return in;
}

// Parameters t in (0, 1) where segment a-b crosses the boundary's edges.
void crossings(const Polyline3& poly, const Vec3& a, const Vec3& b, std::vector<double>& ts) {
    const double dx = b.x - a.x, dy = b.y - a.y;
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const Vec3 &p = poly[j], &q = poly[i];
        const double ex = q.x - p.x, ey = q.y - p.y;
        const double den = dx * ey - dy * ex;
        if (den == 0) continue;
        const double t = ((p.x - a.x) * ey - (p.y - a.y) * ex) / den;
        const double u = ((p.x - a.x) * dy - (p.y - a.y) * dx) / den;
        if (t > 0 && t < 1 && u >= 0 && u <= 1) ts.push_back(t);
    }
}

// Convex polygon (xy) clipped by a convex region given as a counter-clockwise triangle.
std::vector<Vec3> clip_convex(std::vector<Vec3> poly, const Vec3 r[3]) {
    for (int e = 0; e < 3 && !poly.empty(); ++e) {
        const Vec3 &p = r[e], &q = r[(e + 1) % 3];
        auto side = [&](const Vec3& v) { return (q.x - p.x) * (v.y - p.y) - (q.y - p.y) * (v.x - p.x); };
        std::vector<Vec3> out;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const Vec3 &a = poly[i], &b = poly[(i + 1) % poly.size()];
            const double sa = side(a), sb = side(b);
            if (sa >= 0) out.push_back(a);
            if ((sa >= 0) != (sb >= 0)) out.push_back(a + (b - a) * (sa / (sa - sb)));
        }
        poly = std::move(out);
    }
    return poly;
}

}  // namespace

Sink clip_to_boundary(const Sink& in, const Polyline3& boundary) {
    Sink out;
    Polyline3 poly = boundary;
    if (poly.size() >= 2 && poly.front() == poly.back()) poly.pop_back();
    if (poly.size() < 3) return out;
    // lines: split at every crossing, keep the pieces whose midpoint is inside
    std::vector<double> ts;
    for (const auto& l : in.lines) {
        Polyline3 cur;
        auto flush = [&] { if (cur.size() >= 2) out.line(cur, l.color); cur.clear(); };
        for (std::size_t i = 0; i + 1 < l.pts.size(); ++i) {
            const Vec3 a = l.pts[i], b = l.pts[i + 1];
            ts.assign({0.0, 1.0});
            crossings(poly, a, b, ts);
            std::sort(ts.begin(), ts.end());
            for (std::size_t k = 0; k + 1 < ts.size(); ++k) {
                if (ts[k + 1] - ts[k] <= 0) continue;
                const Vec3 m = a + (b - a) * ((ts[k] + ts[k + 1]) / 2);
                const Vec3 p0 = a + (b - a) * ts[k], p1 = a + (b - a) * ts[k + 1];
                if (inside(poly, m.x, m.y)) {
                    if (cur.empty() || !(cur.back() == p0)) { flush(); cur.push_back(p0); }
                    cur.push_back(p1);
                } else {
                    flush();
                }
            }
        }
        flush();
    }
    // triangles: the region as convex pieces (even-odd fill), each triangle clipped against them
    if (!in.tri_pts.empty()) {
        std::vector<Vec3> region;
        fill_even_odd({poly}, region);
        for (std::size_t t = 0; t < in.tri_colors.size(); ++t) {
            const Vec3 &a = in.tri_pts[3 * t], &b = in.tri_pts[3 * t + 1], &c = in.tri_pts[3 * t + 2];
            const double area = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if (std::fabs(area) < 1e-15) {  // seen edge-on from above: kept when its centre is inside
                const Vec3 m = (a + b + c) * (1.0 / 3);
                if (inside(poly, m.x, m.y)) out.tri(a, b, c, in.tri_colors[t]);
                continue;
            }
            // z on the triangle's plane, from barycentric coordinates in xy
            auto z_at = [&](const Vec3& p) {
                const double l1 = ((b.x - p.x) * (c.y - p.y) - (b.y - p.y) * (c.x - p.x)) / area;
                const double l2 = ((c.x - p.x) * (a.y - p.y) - (c.y - p.y) * (a.x - p.x)) / area;
                return l1 * a.z + l2 * b.z + (1 - l1 - l2) * c.z;
            };
            const double lo_x = std::min({a.x, b.x, c.x}), hi_x = std::max({a.x, b.x, c.x});
            const double lo_y = std::min({a.y, b.y, c.y}), hi_y = std::max({a.y, b.y, c.y});
            for (std::size_t r = 0; r + 2 < region.size(); r += 3) {
                const Vec3 rt[3] = {region[r], region[r + 1], region[r + 2]};
                if (std::max({rt[0].x, rt[1].x, rt[2].x}) < lo_x || std::min({rt[0].x, rt[1].x, rt[2].x}) > hi_x ||
                    std::max({rt[0].y, rt[1].y, rt[2].y}) < lo_y || std::min({rt[0].y, rt[1].y, rt[2].y}) > hi_y)
                    continue;
                std::vector<Vec3> piece = area > 0 ? std::vector<Vec3>{a, b, c} : std::vector<Vec3>{a, c, b};
                piece = clip_convex(std::move(piece), rt);
                if (piece.size() < 3) continue;
                for (auto& p : piece) p.z = z_at(p);
                if (area < 0) std::reverse(piece.begin(), piece.end());  // keep the original orientation
                for (std::size_t k = 1; k + 1 < piece.size(); ++k) out.tri(piece[0], piece[k], piece[k + 1], in.tri_colors[t]);
            }
        }
    }
    for (std::size_t i = 0; i < in.points.size(); ++i)
        if (inside(poly, in.points[i].x, in.points[i].y)) out.point(in.points[i], in.point_colors[i]);
    return out;
}

}  // namespace cadcore::formats::autocad
