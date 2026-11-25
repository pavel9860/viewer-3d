// Entity geometry. Coordinates: LINE, POINT, 3DFACE, ELLIPSE, SPLINE, MESH, 3D polylines and
// RAY/XLINE are in world (block) coordinates; CIRCLE, ARC, LWPOLYLINE, 2D POLYLINE, SOLID, TRACE
// and HATCH are in the object coordinate system of their extrusion direction (210).
#include "formats/dxf/draw.hpp"

#include "formats/dxf/colors.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cadcore::formats::autocad {
namespace {

Vec3 extrusion(const Entity& e) { return e.point(210, {0, 0, 1}); }

// Sweep from start to end angle (radians) counter-clockwise, in (0, 2 pi].
double ccw_sweep(double a0, double a1) {
    double s = std::fmod(a1 - a0, 2 * kPi);
    if (s <= 1e-12) s += 2 * kPi;
    return s;
}

Polyline3 to_wcs(const Ocs& o, Polyline3 p) {
    if (!o.identity())
        for (auto& v : p) v = o.to_wcs(v);
    return p;
}

// A curve with optional thickness: the curve, and when thick, its copy along `t` and the
// ruled surface between them; `cap` fills both ends (closed planar outlines: circles).
void emit_curve(Sink& out, Polyline3 c, Vec3 t, bool cap, Color col) {
    if (length(t) == 0) { out.line(std::move(c), col); return; }
    if (cap) {  // orient the outline counter-clockwise about t so walls and caps face outwards
        Vec3 n{};
        for (std::size_t i = 0; i < c.size(); ++i) n = n + cross(c[i], c[(i + 1) % c.size()]);
        if (dot(n, t) < 0) std::reverse(c.begin(), c.end());
    }
    Polyline3 top = c;
    for (auto& v : top) v = v + t;
    for (std::size_t i = 0; i + 1 < c.size(); ++i) out.quad(c[i], c[i + 1], top[i + 1], top[i], col);
    const std::size_t m = c.size() >= 2 && c.front() == c.back() ? c.size() - 1 : c.size();  // closed: last = first
    if (cap && m >= 3) {
        for (std::size_t i = 1; i + 1 < m; ++i) {
            out.tri(c[0], c[i + 1], c[i], col);
            out.tri(top[0], top[i], top[i + 1], col);
        }
    }
    out.line(std::move(c), col);
    out.line(std::move(top), col);
}

struct Vertex2 {
    Vec3 p;
    double w0 = 0, w1 = 0, bulge = 0;
};

// A 2D polyline (OCS) with bulges and widths: a centreline, or a filled band where widths are set.
void emit_polyline2d(Sink& out, const Ocs& o, const std::vector<Vertex2>& v, bool closed, double thickness, Color col) {
    if (v.empty()) return;
    const std::size_t segs = closed ? v.size() : v.size() - 1;
    bool wide = false;
    for (std::size_t i = 0; i < segs; ++i) wide |= v[i].w0 != 0 || v[i].w1 != 0;
    if (!wide) {
        Polyline3 c{v[0].p};
        for (std::size_t i = 0; i < segs; ++i) append_bulge(c, v[i].p, v[(i + 1) % v.size()].p, v[i].bulge);
        if (v.size() == 1) c.push_back(v[0].p);  // a one-vertex polyline is a point
        emit_curve(out, to_wcs(o, std::move(c)), o.az * thickness, false, col);
        return;
    }
    // band: per segment, sample the centreline and offset each sample by half the interpolated width
    for (std::size_t i = 0; i < segs; ++i) {
        const Vec3 a = v[i].p, b = v[(i + 1) % v.size()].p;
        Polyline3 c{a};
        append_bulge(c, a, b, v[i].bulge);
        if (c.size() < 2) continue;
        std::vector<Vec3> left, right;
        for (std::size_t k = 0; k < c.size(); ++k) {
            const Vec3 d = k + 1 < c.size() ? c[k + 1] - c[k] : c[k] - c[k - 1];
            const double l = std::sqrt(d.x * d.x + d.y * d.y);
            const Vec3 n = l > 0 ? Vec3{-d.y / l, d.x / l, 0} : Vec3{0, 0, 0};
            const double w = (v[i].w0 + (v[i].w1 - v[i].w0) * double(k) / double(c.size() - 1)) / 2;
            left.push_back(c[k] + n * w);
            right.push_back(c[k] - n * w);
        }
        for (std::size_t k = 0; k + 1 < c.size(); ++k)
            out.quad(o.to_wcs(right[k]), o.to_wcs(right[k + 1]), o.to_wcs(left[k + 1]), o.to_wcs(left[k]), col);
        if (thickness != 0) {  // a thick wide polyline: AutoCAD draws the band's outline walls
            Polyline3 outline = right;
            outline.insert(outline.end(), left.rbegin(), left.rend());
            outline.push_back(right.front());
            emit_curve(out, to_wcs(o, outline), o.az * thickness, false, col);
        }
    }
}

bool draw_lwpolyline(const Entity& e, Sink& out, Color col) {
    const Ocs o(extrusion(e));
    const double z = e.num(38);
    const double cw = e.num(43);
    std::vector<Vertex2> v;
    for (const Group* g = e.begin; g != e.end; ++g) {
        switch (g->code) {
            case 10: v.push_back({{g->num, 0, z}, cw, cw, 0}); break;
            case 20: if (!v.empty()) v.back().p.y = g->num; break;
            case 40: if (!v.empty()) v.back().w0 = g->num; break;
            case 41: if (!v.empty()) v.back().w1 = g->num; break;
            case 42: if (!v.empty()) v.back().bulge = g->num; break;
            default: break;
        }
    }
    emit_polyline2d(out, o, v, (e.integer(70) & 1) != 0, e.num(39), col);
    return true;
}

void draw_polyface(const Entity& e, Sink& out, Color col, std::vector<std::string>& losses) {
    std::vector<Vec3> pts;
    std::size_t bad = 0;
    for (const auto& v : e.followers) {
        const long long f = v.integer(70);
        if ((f & 128) && (f & 64)) { pts.push_back(v.point(10)); continue; }
        if (!(f & 128)) continue;
        long long idx[4] = {v.integer(71), v.integer(72), v.integer(73), v.integer(74)};
        std::vector<Vec3> face;
        for (long long i : idx) {
            if (i == 0) continue;
            const long long k = (i < 0 ? -i : i) - 1;
            if (k >= (long long)pts.size()) { ++bad; face.clear(); break; }
            face.push_back(pts[std::size_t(k)]);
        }
        if (face.size() == 3) out.tri(face[0], face[1], face[2], col);
        else if (face.size() == 4) out.quad(face[0], face[1], face[2], face[3], col);
        else if (face.size() == 2) out.line({face[0], face[1]}, col);
    }
    if (bad) losses.push_back("POLYLINE polyface mesh: face references a missing vertex");
}

void draw_polygon_mesh(const Entity& e, Sink& out, Color col, std::vector<std::string>& losses) {
    const long long flags = e.integer(70);
    std::vector<Vec3> pts;
    for (const auto& v : e.followers) pts.push_back(v.point(10));
    long long m = e.integer(71), n = e.integer(72);
    // a smoothed mesh (75 = 5, 6, 8) stores the fitted surface of 73 x 74 vertices after the control net
    const long long smooth = e.integer(75);
    if (smooth == 5 || smooth == 6 || smooth == 8) {
        const long long sm = e.integer(73), sn = e.integer(74);
        if (sm > 0 && sn > 0 && std::size_t(m * n + sm * sn) == pts.size()) {
            pts.erase(pts.begin(), pts.begin() + std::ptrdiff_t(m * n));
            m = sm;
            n = sn;
        }
    }
    if (m < 1 || n < 1 || std::size_t(m * n) != pts.size()) {
        losses.push_back("POLYLINE polygon mesh: vertex count is not M x N");
        return;
    }
    const bool cm = flags & 1, cn = flags & 32;
    auto at = [&](long long i, long long j) { return pts[std::size_t((i % m) * n + (j % n))]; };
    if (m == 1 || n == 1) {
        Polyline3 c(pts.begin(), pts.end());
        out.line(std::move(c), col);
        return;
    }
    for (long long i = 0; i < (cm ? m : m - 1); ++i)
        for (long long j = 0; j < (cn ? n : n - 1); ++j) out.quad(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), col);
}

void draw_polyline(const Entity& e, Sink& out, Color col, std::vector<std::string>& losses) {
    const long long flags = e.integer(70);
    if (flags & 64) { draw_polyface(e, out, col, losses); return; }
    if (flags & 16) { draw_polygon_mesh(e, out, col, losses); return; }
    const bool closed = flags & 1;
    if (flags & 8) {  // 3D polyline (WCS); spline frame control points (vertex flag 16) are not drawn
        Polyline3 c;
        for (const auto& v : e.followers)
            if (!(v.integer(70) & 16)) c.push_back(v.point(10));
        if (closed && c.size() > 2) c.push_back(c.front());
        out.line(std::move(c), col);
        return;
    }
    const Ocs o(extrusion(e));
    const double z = e.num(30, e.num(38));
    const double w0 = e.num(40), w1 = e.num(41);
    std::vector<Vertex2> v;
    for (const auto& f : e.followers) {
        if (f.integer(70) & 16) continue;
        Vec3 p = f.point(10);
        p.z = z;
        v.push_back({p, f.num(40, w0), f.num(41, w1), f.num(42)});
    }
    emit_polyline2d(out, o, v, closed, e.num(39), col);
}

void draw_solid(const Entity& e, Sink& out, Color col) {
    const Ocs o(extrusion(e));
    const double z = e.num(38);
    Vec3 p[4] = {e.point(10, {0, 0, z}), e.point(11, {0, 0, z}), e.point(12, {0, 0, z}), {}};
    p[3] = e.has(13) ? e.point(13, {0, 0, z}) : p[2];
    for (auto& q : p) q = o.to_wcs(q);
    // corner order 1-2-4-3 (a "Z")
    const bool tri = p[3] == p[2];
    const Vec3 t = o.az * e.num(39);
    if (length(t) == 0) {
        if (tri) out.tri(p[0], p[1], p[2], col);
        else out.quad(p[0], p[1], p[3], p[2], col);
        return;
    }
    // with thickness: a closed prism (outline closed by repeating its first corner)
    Polyline3 outline = tri ? Polyline3{p[0], p[1], p[2]} : Polyline3{p[0], p[1], p[3], p[2]};
    outline.push_back(outline.front());
    emit_curve(out, outline, t, true, col);
}

void draw_3dface(const Entity& e, Sink& out, Color col) {
    const Vec3 a = e.point(10), b = e.point(11), c = e.point(12);
    const Vec3 d = e.has(13) ? e.point(13) : c;
    if (d == c) out.tri(a, b, c, col);
    else out.quad(a, b, c, d, col);
}

// SPLINE (and HELIX, which carries spline data): control points and knots, or fit points only.
bool spline_points(const Entity& e, Polyline3& out, std::vector<std::string>& losses) {
    std::vector<double> knots, weights;
    std::vector<Vec3> ctrl, fit;
    for (const Group* g = e.begin; g != e.end; ++g) {
        if (g->code == 100 && g->str == "AcDbHelix") break;  // HELIX: its own 40/41 groups (radius, turns) follow
        switch (g->code) {
            case 40: knots.push_back(g->num); break;
            case 41: weights.push_back(g->num); break;
            case 10: ctrl.push_back({g->num, 0, 0}); break;
            case 20: if (!ctrl.empty()) ctrl.back().y = g->num; break;
            case 30: if (!ctrl.empty()) ctrl.back().z = g->num; break;
            case 11: fit.push_back({g->num, 0, 0}); break;
            case 21: if (!fit.empty()) fit.back().y = g->num; break;
            case 31: if (!fit.empty()) fit.back().z = g->num; break;
            default: break;
        }
    }
    const int degree = int(e.integer(71, 3));
    if (!weights.empty() && weights.size() != ctrl.size()) weights.clear();  // AutoCAD writes no weights for non-rational
    if (!ctrl.empty()) {
        if (sample_nurbs(out, degree, knots, ctrl, weights)) return true;
        losses.push_back("SPLINE: inconsistent control points / knots");
        return false;
    }
    if (fit.size() >= 2) {
        std::vector<double> U;
        std::vector<Vec3> P;
        if (interpolate_spline(fit, degree, e.point(12), e.point(13), U, P) && sample_nurbs(out, int(U.size() - P.size() - 1), U, P, {}))
            return true;
        losses.push_back("SPLINE: fit points could not be interpolated");
        return false;
    }
    losses.push_back("SPLINE without control or fit points");
    return false;
}

void draw_mesh(const Entity& e, Sink& out, Color col, std::vector<std::string>& losses) {
    std::vector<Vec3> pts;
    std::vector<long long> faces;
    enum { None, Verts, Faces, Done } state = None;
    long long expect = 0;
    for (const Group* g = e.begin; g != e.end; ++g) {
        if (g->code == 92 && state == None) { state = Verts; expect = g->ival; continue; }
        if (state == Verts) {
            if (g->code == 10) pts.push_back({g->num, 0, 0});
            else if (g->code == 20 && !pts.empty()) pts.back().y = g->num;
            else if (g->code == 30 && !pts.empty()) pts.back().z = g->num;
            else if (g->code == 93) { state = Faces; expect = g->ival; }
            continue;
        }
        if (state == Faces) {
            if (g->code == 90 && expect > 0) { faces.push_back(g->ival); --expect; }
            else if (g->code != 90) state = Done;
        }
    }
    std::size_t bad = 0;
    for (std::size_t i = 0; i < faces.size();) {
        const long long n = faces[i++];
        if (n < 0 || i + std::size_t(n) > faces.size()) { ++bad; break; }
        std::vector<Vec3> f;
        for (long long k = 0; k < n; ++k) {
            const long long j = faces[i + std::size_t(k)];
            if (j < 0 || j >= (long long)pts.size()) { ++bad; f.clear(); break; }
            f.push_back(pts[std::size_t(j)]);
        }
        i += std::size_t(n);
        for (std::size_t k = 1; k + 1 < f.size(); ++k) out.tri(f[0], f[k], f[k + 1], col);
    }
    if (bad) losses.push_back("MESH: face list out of range");
    if (e.integer(91) > 0) losses.push_back("MESH: subdivision smoothing not applied (level-0 cage shown)");
}

// --- HATCH -------------------------------------------------------------------------------------

struct HatchData {
    std::vector<Polyline3> loops;  // OCS, closed
    bool solid = false;
    int style = 0;                 // 75: 0 normal (odd parity), 1 outer, 2 ignore islands
    struct Family { double angle; Vec3 base, offset; std::vector<double> dashes; };
    std::vector<Family> families;
};

// One boundary edge: the candidate point sequences it may stand for. Clockwise arc and ellipse
// edges have two readings of their stored angles (as stored, or mirrored); the boundary's
// connectivity decides between them.
using EdgeOptions = std::vector<Polyline3>;

double gap(const Vec3& a, const Vec3& b) { return length(a - b); }

// Joins edges head to tail: each edge takes the reading and direction that continues from the
// previous edge's end (the first edge: the one that meets the second edge).
Polyline3 chain_edges(const std::vector<EdgeOptions>& edges) {
    Polyline3 loop;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const Polyline3* best = nullptr;
        bool reversed = false;
        double cost = std::numeric_limits<double>::infinity();
        for (const auto& opt : edges[i]) {
            if (opt.empty()) continue;
            for (int rev = 0; rev < 2; ++rev) {
                const Vec3 start = rev ? opt.back() : opt.front(), end = rev ? opt.front() : opt.back();
                double c;
                if (!loop.empty()) {
                    c = gap(loop.back(), start);
                } else if (edges.size() > 1) {  // first edge: its end should meet an end of the next edge
                    c = std::numeric_limits<double>::infinity();
                    for (const auto& nx : edges[1])
                        if (!nx.empty()) c = std::min({c, gap(end, nx.front()), gap(end, nx.back())});
                } else {
                    c = rev ? 1 : 0;
                }
                if (!best || c < cost) { cost = c; best = &opt; reversed = rev != 0; }
            }
        }
        if (!best) continue;
        Polyline3 seg = *best;
        if (reversed) std::reverse(seg.begin(), seg.end());
        if (!loop.empty() && gap(loop.back(), seg.front()) <= 1e-9 * (1 + length(seg.front()))) seg.erase(seg.begin());
        loop.insert(loop.end(), seg.begin(), seg.end());
    }
    return loop;
}

// Parses the boundary paths and pattern lines; the group order is fixed by the DXF reference.
bool parse_hatch(const Entity& e, HatchData& h, std::vector<std::string>& losses) {
    const Group* g = e.begin;
    const Group* end = e.end;
    auto next = [&](int code) -> const Group* {  // the next group if it has `code`
        if (g != end && g->code == code) return g++;
        return nullptr;
    };
    auto need = [&](int code) -> double {
        const Group* q = next(code);
        if (!q) throw 0;
        return q->num;
    };
    try {
        while (g != end && g->code != 91) {
            if (g->code == 70) h.solid = g->ival == 1;
            ++g;
        }
        if (g == end) return false;
        const long long paths = g++->ival;
        for (long long pi = 0; pi < paths; ++pi) {
            while (g != end && g->code != 92) ++g;
            if (g == end) throw 0;
            const long long type = g++->ival;
            Polyline3 loop;
            if (type & 2) {  // polyline path
                const Group* hb = next(72);
                const bool bulges = hb && hb->ival != 0;
                next(73);
                const long long n = (long long)need(93);
                std::vector<Vertex2> v;
                for (long long k = 0; k < n; ++k) {
                    Vertex2 x;
                    x.p.x = need(10);
                    x.p.y = need(20);
                    if (bulges) { const Group* b = next(42); x.bulge = b ? b->num : 0; }
                    v.push_back(x);
                }
                if (!v.empty()) {
                    loop.push_back(v[0].p);
                    for (std::size_t k = 0; k < v.size(); ++k) append_bulge(loop, v[k].p, v[(k + 1) % v.size()].p, v[k].bulge);
                }
            } else {
                const long long ne = (long long)need(93);
                std::vector<EdgeOptions> edges;
                for (long long k = 0; k < ne; ++k) {
                    const long long et = (long long)need(72);
                    EdgeOptions opts(1);
                    if (et == 1) {
                        Vec3 a{need(10), need(20), 0}, b{need(11), need(21), 0};
                        opts[0] = {a, b};
                    } else if (et == 2) {
                        Vec3 c{need(10), need(20), 0};
                        const double r = need(40), s = need(50) * kDeg, t = need(51) * kDeg;
                        const bool ccw = need(73) != 0;
                        if (ccw) {
                            append_arc(opts[0], c, r, s, ccw_sweep(s, t), false);
                        } else {
                            append_arc(opts[0], c, r, s, -ccw_sweep(t, s), false);  // as stored, clockwise
                            opts.emplace_back();
                            append_arc(opts[1], c, r, -s, -ccw_sweep(s, t), false);  // mirrored angles
                        }
                    } else if (et == 3) {
                        Vec3 c{need(10), need(20), 0}, u{need(11), need(21), 0};
                        const double ratio = need(40), s = need(50) * kDeg, t = need(51) * kDeg;
                        const bool ccw = need(73) != 0;
                        const Vec3 v{-u.y * ratio, u.x * ratio, 0};
                        if (ccw) {
                            append_ellipse(opts[0], c, u, v, s, s + ccw_sweep(s, t), false);
                        } else {
                            Polyline3 a, b;
                            append_ellipse(a, c, u, v, s - ccw_sweep(t, s), s, false);   // as stored, clockwise
                            append_ellipse(b, c, u, v, -s - ccw_sweep(s, t), -s, false);  // mirrored
                            opts[0].assign(a.rbegin(), a.rend());
                            opts.emplace_back(b.rbegin(), b.rend());
                        }
                    } else if (et == 4) {
                        const int degree = int(need(94));
                        const bool rational = need(73) != 0;
                        next(74);
                        const long long nk = (long long)need(95), nc = (long long)need(96);
                        std::vector<double> U, W;
                        std::vector<Vec3> P;
                        for (long long i = 0; i < nk; ++i) U.push_back(need(40));
                        for (long long i = 0; i < nc; ++i) {
                            Vec3 p{need(10), need(20), 0};
                            P.push_back(p);
                            if (rational) { const Group* w = next(42); W.push_back(w ? w->num : 1); }
                        }
                        std::vector<Vec3> fit;
                        Vec3 t0, t1;
                        if (const Group* nf = next(97)) {
                            for (long long i = 0; i < nf->ival; ++i) fit.push_back({need(11), need(21), 0});
                            if (next(12)) { t0.x = g[-1].num; t0.y = need(22); }
                            if (next(13)) { t1.x = g[-1].num; t1.y = need(23); }
                        }
                        bool ok = !P.empty() && sample_nurbs(opts[0], degree, U, P, W);
                        if (!ok && fit.size() >= 2) {  // fit points only: the interpolating spline
                            std::vector<double> FU;
                            std::vector<Vec3> FP;
                            opts[0].clear();
                            ok = interpolate_spline(fit, degree, t0, t1, FU, FP) && sample_nurbs(opts[0], int(FU.size() - FP.size() - 1), FU, FP, {});
                        }
                        if (!ok) losses.push_back("HATCH: spline edge with inconsistent knots");
                    } else {
                        throw 0;
                    }
                    edges.push_back(std::move(opts));
                }
                loop = chain_edges(edges);
            }
            if (const Group* nsrc = next(97))
                for (long long k = 0; k < nsrc->ival; ++k) next(330);
            if (loop.size() >= 3) h.loops.push_back(std::move(loop));
        }
        while (g != end && g->code != 75 && g->code != 78) ++g;
        if (const Group* st = next(75)) h.style = int(st->ival);
        // pattern definition lines
        while (g != end && g->code != 78) ++g;
        if (g != end) {
            const long long nl = g++->ival;
            for (long long k = 0; k < nl; ++k) {
                HatchData::Family f;
                f.angle = need(53) * kDeg;
                f.base = {need(43), need(44), 0};
                f.offset = {need(45), need(46), 0};
                const long long nd = (long long)need(79);
                for (long long d = 0; d < nd; ++d) f.dashes.push_back(need(49));
                h.families.push_back(std::move(f));
            }
        }
    } catch (int) {
        losses.push_back("HATCH: boundary data does not follow the group order");
        return false;
    }
    // hatch style: "outer" fills between the outermost loops and their first islands, "ignore"
    // fills the outermost loops only
    if (h.style == 1 || h.style == 2) {
        const std::vector<int> depth = nesting_depths(h.loops);
        std::vector<Polyline3> kept;
        for (std::size_t i = 0; i < h.loops.size(); ++i)
            if (depth[i] <= (h.style == 1 ? 1 : 0)) kept.push_back(std::move(h.loops[i]));
        h.loops = std::move(kept);
    }
    return true;
}

// Pattern lines clipped to the hatch region (even-odd), split into dashes.
void hatch_pattern(const HatchData& h, const Ocs& o, double z, Sink& out, Color col, std::vector<std::string>& losses) {
    double lo[2] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    for (const auto& l : h.loops)
        for (const auto& p : l) { lo[0] = std::min(lo[0], p.x); lo[1] = std::min(lo[1], p.y); hi[0] = std::max(hi[0], p.x); hi[1] = std::max(hi[1], p.y); }
    std::size_t total = 0;
    for (const auto& f : h.families) {
        const Vec3 d{std::cos(f.angle), std::sin(f.angle), 0};
        const Vec3 nrm{-d.y, d.x, 0};
        const double spacing = dot(f.offset, nrm);
        if (std::fabs(spacing) < 1e-12) continue;
        // range of line numbers k whose line base + k offset crosses the bounding box
        double kmin = std::numeric_limits<double>::infinity(), kmax = -kmin;
        for (double cx : {lo[0], hi[0]})
            for (double cy : {lo[1], hi[1]}) {
                const double k = dot(Vec3{cx, cy, 0} - f.base, nrm) / spacing;
                kmin = std::min(kmin, k);
                kmax = std::max(kmax, k);
            }
        const long long k0 = (long long)std::floor(kmin), k1 = (long long)std::ceil(kmax);
        if (k1 - k0 > 200000) { losses.push_back("HATCH: pattern too dense to draw (more than 200,000 lines)"); return; }
        double period = 0;
        for (double x : f.dashes) period += std::fabs(x);
        for (long long k = k0; k <= k1; ++k) {
            const Vec3 b = f.base + f.offset * double(k);
            // crossings of the line b + t d with every loop edge
            std::vector<double> ts;
            for (const auto& l : h.loops)
                for (std::size_t i = 0; i < l.size(); ++i) {
                    const Vec3 p = l[i], q = l[(i + 1) % l.size()];
                    const double sp = dot(p - b, nrm), sq = dot(q - b, nrm);
                    if ((sp > 0) == (sq > 0)) continue;
                    const double a = sp / (sp - sq);
                    const Vec3 x = p + (q - p) * a;
                    ts.push_back(dot(x - b, d));
                }
            std::sort(ts.begin(), ts.end());
            for (std::size_t i = 0; i + 1 < ts.size(); i += 2) {
                const double t0 = ts[i], t1 = ts[i + 1];
                auto seg = [&](double a, double c) {
                    Polyline3 s{b + d * a, b + d * c};
                    for (auto& v : s) v.z = z;
                    out.line(to_wcs(o, std::move(s)), col);
                    ++total;
                };
                if (f.dashes.empty() || period <= 0) { seg(t0, t1); continue; }
                // dash pattern starts at the line's base point (t = 0)
                double t = std::floor(t0 / period) * period;
                std::size_t guard = 0;
                while (t < t1 && guard++ < 1000000) {
                    for (double dl : f.dashes) {
                        const double len = std::fabs(dl);
                        if (dl > 0) { const double a = std::max(t, t0), c = std::min(t + len, t1); if (c > a) seg(a, c); }
                        else if (dl == 0 && t >= t0 && t <= t1) { Vec3 p = b + d * t; p.z = z; out.point(o.to_wcs(p), col); }
                        t += len;
                    }
                }
            }
        }
    }
    (void)total;
}

void draw_hatch(const Entity& e, Sink& out, Color col, std::vector<std::string>& losses) {
    HatchData h;
    if (!parse_hatch(e, h, losses)) return;
    const Ocs o(extrusion(e));
    const double z = e.num(30);
    if (h.loops.empty()) { losses.push_back("HATCH without a boundary"); return; }
    if (h.solid) {
        std::vector<Vec3> tris;
        fill_even_odd(h.loops, tris);
        for (std::size_t i = 0; i + 2 < tris.size(); i += 3) {
            Vec3 p[3] = {tris[i], tris[i + 1], tris[i + 2]};
            for (auto& q : p) { q.z = z; q = o.to_wcs(q); }
            out.tri(p[0], p[1], p[2], col);
        }
        if (e.integer(450) == 1) losses.push_back("HATCH: gradient fill shown in one colour");
        return;
    }
    if (h.families.empty()) { losses.push_back("HATCH: pattern fill without pattern lines"); return; }
    hatch_pattern(h, o, z, out, col, losses);
}

// --- MLINE ------------------------------------------------------------------------------------
// Vertices (11) with segment (12) and miter (13) directions; per vertex and element, parameters
// (74 count, 41 values): the offset along the miter, then the distance to the first dash and the
// dash / gap lengths. Element k at vertex v lies at vertex + miter * offset.
void draw_mline(const Entity& e, const DrawContext& ctx, Sink& out, std::vector<std::string>& losses) {
    struct V { Vec3 p, miter; std::vector<std::vector<double>> params; };
    std::vector<V> verts;
    int axis = 0;  // 0 none, 11 vertex, 13 miter
    for (const Group* g = e.begin; g != e.end; ++g) {
        switch (g->code) {
            case 11: verts.push_back({{g->num, 0, 0}, {}, {}}); axis = 11; break;
            case 21: if (axis == 11 && !verts.empty()) verts.back().p.y = g->num; break;
            case 31: if (axis == 11 && !verts.empty()) verts.back().p.z = g->num; break;
            case 12: axis = 12; break;
            case 13: if (!verts.empty()) verts.back().miter.x = g->num; axis = 13; break;
            case 23: if (axis == 13 && !verts.empty()) verts.back().miter.y = g->num; break;
            case 33: if (axis == 13 && !verts.empty()) verts.back().miter.z = g->num; break;
            case 74: if (!verts.empty()) verts.back().params.emplace_back(); axis = 74; break;
            case 41: if (axis == 74 && !verts.empty() && !verts.back().params.empty()) verts.back().params.back().push_back(g->num); break;
            case 75: axis = 75; break;
            default: break;
        }
    }
    if (verts.size() < 2) { losses.push_back("MLINE with fewer than 2 vertices"); return; }
    const std::size_t ne = std::size_t(e.integer(73, 0));
    const bool closed = (e.integer(71) & 2) != 0;
    const MlineStyle* st = nullptr;
    if (ctx.mline_styles) {
        std::string key(e.str(2, "STANDARD"));
        for (auto& c : key) if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
        if (auto it = ctx.mline_styles->find(key); it != ctx.mline_styles->end()) st = &it->second;
    }
    auto color_of = [&](int aci) { return aci == 256 ? ctx.layer_color : aci == 0 || aci == 257 ? ctx.color : aci_color(aci); };
    std::vector<Polyline3> lines(ne);
    for (const V& v : verts) {
        if (v.params.size() != ne) { losses.push_back("MLINE: element parameters do not match the element count"); return; }
        for (std::size_t k = 0; k < ne; ++k) {
            const auto& pr = v.params[k];
            if (pr.empty()) { losses.push_back("MLINE: element without an offset"); return; }
            if (pr.size() > 2 || (pr.size() == 2 && pr[1] != 0)) { losses.push_back("MLINE: element breaks (dashes) not applied"); return; }
            lines[k].push_back(v.p + v.miter * pr[0]);
        }
    }
    for (std::size_t k = 0; k < ne; ++k) {
        if (closed) lines[k].push_back(lines[k].front());
        const int aci = st && k < st->elements.size() ? st->elements[k].aci : 256;
        out.line(lines[k], color_of(aci));
    }
    if (!st) return;
    // outermost elements: largest and smallest offset
    std::size_t hi = 0, lo = 0;
    for (std::size_t k = 0; k < ne && k < st->elements.size(); ++k) {
        if (st->elements[k].offset > st->elements[hi].offset) hi = k;
        if (st->elements[k].offset < st->elements[lo].offset) lo = k;
    }
    if (ne >= 2 && (st->flags & 1)) {  // fill between the outermost elements
        for (std::size_t i = 0; i + 1 < lines[hi].size(); ++i) out.quad(lines[lo][i], lines[lo][i + 1], lines[hi][i + 1], lines[hi][i], color_of(st->fill_aci));
    }
    if (!closed && ne >= 2) {
        const long long flags = e.integer(71);
        const Color cap = color_of(256);
        if ((st->flags & 16) && !(flags & 4)) out.line({lines[lo].front(), lines[hi].front()}, cap);
        if ((st->flags & 256) && !(flags & 8)) out.line({lines[lo].back(), lines[hi].back()}, cap);
        if (st->flags & (32 | 64 | 512 | 1024)) losses.push_back("MLINE: round caps / inner arcs not drawn");
    }
}

// --- infinite lines ------------------------------------------------------------------------------

void draw_infinite(const Entity& e, const DrawContext& ctx, Sink& out, Color col, bool ray) {
    // shown as a segment reaching past the drawing extents on both sides (AutoCAD draws it to the
    // screen edge); a drawing without extents gets a segment of 1000 units around the base point
    const Vec3 p = e.point(10), d = normalized(e.point(11, {1, 0, 0}));
    double reach = 1000;
    if (ctx.has_extents) {
        const Vec3 mid = (ctx.extmin + ctx.extmax) * 0.5;
        reach = std::max(1.0, length(ctx.extmax - ctx.extmin) + length(p - mid));
    }
    out.line({ray ? p : p - d * reach, p + d * reach}, col);
    if (ctx.notes) ctx.notes->push_back("RAY/XLINE: infinite construction lines shown as finite segments");
}

}  // namespace

void Sink::tri(const Vec3& a, const Vec3& b, const Vec3& c, Color col) {
    tri_pts.insert(tri_pts.end(), {a, b, c});
    tri_colors.push_back(col);
}

bool draw_entity(const Entity& e, const DrawContext& ctx, Sink& out, std::vector<std::string>& losses) {
    const Color col = ctx.color;
    const std::string_view t = e.type;
    if (t == "LINE" || t == "3DLINE") {
        emit_curve(out, {e.point(10), e.point(11)}, normalized(extrusion(e)) * e.num(39), false, col);
    } else if (t == "POINT") {
        const Vec3 p = e.point(10);
        const double th = e.num(39);
        if (th != 0) out.line({p, p + normalized(extrusion(e)) * th}, col);
        else out.point(p, col);
    } else if (t == "CIRCLE" || t == "ARC") {
        const Ocs o(extrusion(e));
        const Vec3 c = e.point(10, {0, 0, e.num(38)});
        const double r = e.num(40);
        Polyline3 pts;
        const bool circle = t == "CIRCLE";
        const double a0 = circle ? 0 : e.num(50) * kDeg;
        append_arc(pts, c, r, a0, circle ? 2 * kPi : ccw_sweep(a0, e.num(51) * kDeg), false);
        emit_curve(out, to_wcs(o, std::move(pts)), o.az * e.num(39), circle, col);
    } else if (t == "ELLIPSE") {
        const Vec3 c = e.point(10), u = e.point(11), n = normalized(extrusion(e));
        const Vec3 v = cross(n, u) * e.num(40, 1);
        double t0 = e.num(41, 0), t1 = e.num(42, 2 * kPi);
        if (t1 <= t0) t1 += 2 * kPi;
        Polyline3 pts;
        append_ellipse(pts, c, u, v, t0, t1, false);
        out.line(std::move(pts), col);
    } else if (t == "LWPOLYLINE") {
        draw_lwpolyline(e, out, col);
    } else if (t == "POLYLINE") {
        draw_polyline(e, out, col, losses);
    } else if (t == "3DFACE") {
        draw_3dface(e, out, col);
    } else if (t == "SOLID" || t == "TRACE") {
        draw_solid(e, out, col);
    } else if (t == "SPLINE" || t == "HELIX") {
        Polyline3 pts;
        if (spline_points(e, pts, losses)) out.line(std::move(pts), col);
    } else if (t == "MESH") {
        draw_mesh(e, out, col, losses);
    } else if (t == "HATCH") {
        draw_hatch(e, out, col, losses);
    } else if (t == "RAY" || t == "XLINE") {
        draw_infinite(e, ctx, out, col, t == "RAY");
    } else if (t == "MLINE") {
        draw_mline(e, ctx, out, losses);
    } else if (t == "LEADER") {
        Polyline3 pts;
        for (const Group* g = e.begin; g != e.end; ++g) {
            if (g->code == 10) pts.push_back({g->num, 0, 0});
            else if (g->code == 20 && !pts.empty()) pts.back().y = g->num;
            else if (g->code == 30 && !pts.empty()) pts.back().z = g->num;
        }
        if (e.integer(72) == 1 && pts.size() >= 2) {  // spline path through the vertices
            std::vector<double> U;
            std::vector<Vec3> P;
            Polyline3 s;
            if (interpolate_spline(pts, 3, {}, {}, U, P) && sample_nurbs(s, int(U.size() - P.size() - 1), U, P, {})) pts = std::move(s);
        }
        out.line(std::move(pts), col);  // the arrowhead is added by the decoder (it may be a block)
    } else {
        return false;
    }
    return true;
}

}  // namespace cadcore::formats::autocad
