#include "formats/dxf/geometry.hpp"

#include <algorithm>
#include <limits>

namespace cadcore::formats::autocad {

Ocs::Ocs(Vec3 n) {
    az = normalized(n);
    const double limit = 1.0 / 64.0;
    ax = std::fabs(az.x) < limit && std::fabs(az.y) < limit ? normalized(cross({0, 1, 0}, az)) : normalized(cross({0, 0, 1}, az));
    ay = normalized(cross(az, ax));
}

int arc_segments(double sweep) {
    static const double step = 2 * std::acos(1 - kChordTolerance);  // ~0.0894 rad
    const double a = std::fabs(sweep);
    int n = int(std::ceil(a / step - 1e-9));
    // a full circle: 72 segments (5 degrees), within tolerance for display and a round number
    n = std::max(n, int(std::ceil(a / (2 * kPi) * 72 - 1e-9)));
    return std::clamp(n, 1, 100000);
}

void append_arc(Polyline3& out, Vec3 c, double r, double a0, double sweep, bool skip_first) {
    const int n = arc_segments(sweep);
    for (int i = skip_first ? 1 : 0; i <= n; ++i) {
        const double a = a0 + sweep * i / n;
        out.push_back({c.x + r * std::cos(a), c.y + r * std::sin(a), c.z});
    }
}

void append_bulge(Polyline3& out, Vec3 p0, Vec3 p1, double b) {
    const double dx = p1.x - p0.x, dy = p1.y - p0.y;
    const double chord = std::sqrt(dx * dx + dy * dy);
    if (b == 0 || chord == 0 || !std::isfinite(b)) { out.push_back(p1); return; }
    const double sweep = 4 * std::atan(b);           // signed included angle
    const double r = chord / (2 * std::sin(std::fabs(sweep) / 2));
    // centre: from the chord midpoint along the left normal of p0->p1 by r cos(sweep / 2), which
    // is negative (centre on the right) for clockwise arcs and for counter-clockwise arcs > 180 degrees
    const double h = r * std::cos(sweep / 2) * (b > 0 ? 1 : -1);
    const double mx = (p0.x + p1.x) / 2, my = (p0.y + p1.y) / 2;
    const double nx = -dy / chord, ny = dx / chord;
    const Vec3 c{mx + nx * h, my + ny * h, p0.z};
    const double a0 = std::atan2(p0.y - c.y, p0.x - c.x);
    const int n = arc_segments(sweep);
    for (int i = 1; i < n; ++i) {
        const double a = a0 + sweep * i / n;
        const double t = double(i) / n;
        out.push_back({c.x + r * std::cos(a), c.y + r * std::sin(a), p0.z + (p1.z - p0.z) * t});
    }
    out.push_back(p1);
}

void append_ellipse(Polyline3& out, Vec3 c, Vec3 u, Vec3 v, double t0, double t1, bool skip_first) {
    const double ru = length(u), rv = length(v);
    const double ratio = ru > 0 ? std::min(ru, rv) / std::max(ru, rv) : 1;
    // chord deviation of an ellipse is largest at the flat ends; scale the count by 1/sqrt(ratio)
    int n = arc_segments(t1 - t0);
    if (ratio > 0 && ratio < 1) n = std::min(100000, int(std::ceil(n / std::sqrt(ratio))));
    for (int i = skip_first ? 1 : 0; i <= n; ++i) {
        const double t = t0 + (t1 - t0) * i / n;
        out.push_back(c + u * std::cos(t) + v * std::sin(t));
    }
}

namespace {

// Rational point at parameter t (de Boor in homogeneous space); span s with knots[s] <= t < knots[s+1].
Vec3 de_boor(int p, const std::vector<double>& U, const std::vector<Vec3>& P, const std::vector<double>& W, int s, double t) {
    double d[64][4];
    for (int j = 0; j <= p; ++j) {
        const int i = s - p + j;
        const double w = W.empty() ? 1.0 : W[std::size_t(i)];
        d[j][0] = P[std::size_t(i)].x * w; d[j][1] = P[std::size_t(i)].y * w; d[j][2] = P[std::size_t(i)].z * w; d[j][3] = w;
    }
    for (int r = 1; r <= p; ++r)
        for (int j = p; j >= r; --j) {
            const int i = s - p + j;
            const double den = U[std::size_t(i + p - r + 1)] - U[std::size_t(i)];
            const double a = den == 0 ? 0 : (t - U[std::size_t(i)]) / den;
            for (int k = 0; k < 4; ++k) d[j][k] = (1 - a) * d[j - 1][k] + a * d[j][k];
        }
    const double w = d[p][3] == 0 ? 1 : d[p][3];
    return {d[p][0] / w, d[p][1] / w, d[p][2] / w};
}

}  // namespace

bool sample_nurbs(Polyline3& out, int p, const std::vector<double>& U, const std::vector<Vec3>& P, const std::vector<double>& W) {
    const int n = int(P.size());
    if (p < 1 || p > 63 || n < p + 1 || int(U.size()) != n + p + 1) return false;
    if (!W.empty() && int(W.size()) != n) return false;
    for (std::size_t i = 1; i < U.size(); ++i)
        if (U[i] < U[i - 1]) return false;
    const double lo = U[std::size_t(p)], hi = U[std::size_t(n)];
    if (!(hi > lo)) return false;
    const int per_span = std::max(8, 4 * p);
    bool first = true;
    for (int s = p; s < n; ++s) {
        const double a = U[std::size_t(s)], b = U[std::size_t(s + 1)];
        if (!(b > a)) continue;
        for (int i = first ? 0 : 1; i <= per_span; ++i) {
            const double t = a + (b - a) * i / per_span;
            out.push_back(de_boor(p, U, P, W, s, std::min(t, b)));
        }
        first = false;
    }
    return !first;
}

bool interpolate_spline(const std::vector<Vec3>& Q, int p, Vec3 t0, Vec3 t1, std::vector<double>& U, std::vector<Vec3>& P) {
    const int m = int(Q.size());
    if (m < 2) return false;
    p = std::min(p, m - 1 + (length(t0) > 0) + (length(t1) > 0));
    if (p < 1) return false;
    // chord-length parameters
    std::vector<double> u(std::size_t(m), 0.0);
    double total = 0;
    for (int i = 1; i < m; ++i) total += length(Q[std::size_t(i)] - Q[std::size_t(i - 1)]);
    if (total == 0) return false;
    for (int i = 1; i < m; ++i) u[std::size_t(i)] = u[std::size_t(i - 1)] + length(Q[std::size_t(i)] - Q[std::size_t(i - 1)]) / total;
    u[std::size_t(m - 1)] = 1;
    const bool d0 = length(t0) > 0, d1 = length(t1) > 0;
    const int n = m + d0 + d1;  // control points
    if (n < p + 1) return false;
    // parameters the rows interpolate (position rows plus derivative rows at the ends)
    std::vector<double> rows_u;
    std::vector<int> rows_kind;  // 0 position (index into Q), 1 start tangent, 2 end tangent
    std::vector<int> rows_q;
    rows_u.push_back(0); rows_kind.push_back(0); rows_q.push_back(0);
    if (d0) { rows_u.push_back(0); rows_kind.push_back(1); rows_q.push_back(-1); }
    for (int i = 1; i < m - 1; ++i) { rows_u.push_back(u[std::size_t(i)]); rows_kind.push_back(0); rows_q.push_back(i); }
    if (d1) { rows_u.push_back(1); rows_kind.push_back(2); rows_q.push_back(-1); }
    rows_u.push_back(1); rows_kind.push_back(0); rows_q.push_back(m - 1);
    // knots by averaging (NURBS book eq. 9.8 over the extended parameter list)
    U.assign(std::size_t(n + p + 1), 0.0);
    for (int i = 0; i <= p; ++i) { U[std::size_t(i)] = 0; U[std::size_t(n + i)] = 1; }
    for (int j = 1; j <= n - p - 1; ++j) {
        double s = 0;
        for (int i = j; i < j + p; ++i) s += rows_u[std::size_t(std::min(i, n - 1))];
        U[std::size_t(j + p)] = s / p;
    }
    // basis function values and first derivatives at t
    auto basis = [&](double t, std::vector<double>& N, std::vector<double>& dN) {
        N.assign(std::size_t(n), 0.0);
        dN.assign(std::size_t(n), 0.0);
        int s = p;
        while (s < n - 1 && !(t < U[std::size_t(s + 1)])) ++s;
        // Cox-de Boor triangle of degree p and p-1
        std::vector<std::vector<double>> T(std::size_t(p + 1), std::vector<double>(std::size_t(p + 1), 0.0));
        T[0][0] = 1;
        for (int k = 1; k <= p; ++k)
            for (int j = 0; j <= k; ++j) {
                const int i = s - k + j;
                double v = 0;
                if (j > 0) { double den = U[std::size_t(i + k)] - U[std::size_t(i)]; if (den != 0) v += (t - U[std::size_t(i)]) / den * T[std::size_t(k - 1)][std::size_t(j - 1)]; }
                if (j < k) { double den = U[std::size_t(i + k + 1)] - U[std::size_t(i + 1)]; if (den != 0) v += (U[std::size_t(i + k + 1)] - t) / den * T[std::size_t(k - 1)][std::size_t(j)]; }
                T[std::size_t(k)][std::size_t(j)] = v;
            }
        for (int j = 0; j <= p; ++j) N[std::size_t(s - p + j)] = T[std::size_t(p)][std::size_t(j)];
        for (int j = 0; j <= p; ++j) {
            const int i = s - p + j;
            double v = 0;
            if (j > 0) { double den = U[std::size_t(i + p)] - U[std::size_t(i)]; if (den != 0) v += p / den * T[std::size_t(p - 1)][std::size_t(j - 1)]; }
            if (j < p) { double den = U[std::size_t(i + p + 1)] - U[std::size_t(i + 1)]; if (den != 0) v -= p / den * T[std::size_t(p - 1)][std::size_t(j)]; }
            dN[std::size_t(i)] = v;
        }
    };
    // tangents are directions: scale to the curve's parametric speed (total chord length)
    std::vector<std::vector<double>> A(std::size_t(n), std::vector<double>(std::size_t(n), 0.0));
    std::vector<Vec3> B(static_cast<std::size_t>(n));
    std::vector<double> N, dN;
    for (int r = 0; r < n; ++r) {
        basis(rows_u[std::size_t(r)], N, dN);
        if (rows_kind[std::size_t(r)] == 0) { A[std::size_t(r)] = N; B[std::size_t(r)] = Q[std::size_t(rows_q[std::size_t(r)])]; }
        else { A[std::size_t(r)] = dN; B[std::size_t(r)] = normalized(rows_kind[std::size_t(r)] == 1 ? t0 : t1) * total; }
    }
    // Gaussian elimination with partial pivoting
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r)
            if (std::fabs(A[std::size_t(r)][std::size_t(c)]) > std::fabs(A[std::size_t(piv)][std::size_t(c)])) piv = r;
        if (std::fabs(A[std::size_t(piv)][std::size_t(c)]) < 1e-14) return false;
        std::swap(A[std::size_t(c)], A[std::size_t(piv)]);
        std::swap(B[std::size_t(c)], B[std::size_t(piv)]);
        for (int r = c + 1; r < n; ++r) {
            const double f = A[std::size_t(r)][std::size_t(c)] / A[std::size_t(c)][std::size_t(c)];
            if (f == 0) continue;
            for (int k = c; k < n; ++k) A[std::size_t(r)][std::size_t(k)] -= f * A[std::size_t(c)][std::size_t(k)];
            B[std::size_t(r)] = B[std::size_t(r)] - B[std::size_t(c)] * f;
        }
    }
    P.assign(std::size_t(n), Vec3{});
    for (int r = n - 1; r >= 0; --r) {
        Vec3 s = B[std::size_t(r)];
        for (int k = r + 1; k < n; ++k) s = s - P[std::size_t(k)] * A[std::size_t(r)][std::size_t(k)];
        P[std::size_t(r)] = s * (1.0 / A[std::size_t(r)][std::size_t(r)]);
    }
    return true;
}

// --- region fill --------------------------------------------------------------------------------
namespace {

struct Edge { double x0, y0, x1, y1; };  // y0 < y1

double x_at(const Edge& e, double y) { return e.x0 + (e.x1 - e.x0) * (y - e.y0) / (e.y1 - e.y0); }

bool inside(const Polyline3& l, double x, double y) {
    bool in = false;
    for (std::size_t i = 0, j = l.size() - 1; i < l.size(); j = i++)
        if ((l[i].y > y) != (l[j].y > y) && x < (l[j].x - l[i].x) * (y - l[i].y) / (l[j].y - l[i].y) + l[i].x) in = !in;
    return in;
}

}  // namespace

std::vector<int> nesting_depths(const std::vector<Polyline3>& L) {
    std::vector<int> depth(L.size(), 0);
    for (std::size_t i = 0; i < L.size(); ++i)
        for (std::size_t j = 0; j < L.size(); ++j)
            if (i != j && !L[i].empty() && L[j].size() >= 3 && inside(L[j], L[i][0].x, L[i][0].y)) ++depth[i];
    return depth;
}

void fill_even_odd(const std::vector<Polyline3>& loops, std::vector<Vec3>& tris) {
    std::vector<Edge> edges;
    std::vector<double> ys;
    for (const auto& l : loops) {
        if (l.size() < 3) continue;
        for (std::size_t i = 0; i < l.size(); ++i) {
            const Vec3 &a = l[i], &b = l[(i + 1) % l.size()];
            ys.push_back(a.y);
            if (a.y == b.y || !std::isfinite(a.y + b.y + a.x + b.x)) continue;  // horizontal edges bound no slab
            edges.push_back(a.y < b.y ? Edge{a.x, a.y, b.x, b.y} : Edge{b.x, b.y, a.x, a.y});
        }
    }
    if (edges.empty()) return;
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) { return a.y0 < b.y0; });
    std::vector<const Edge*> active;
    std::size_t next = 0;
    struct Cut { double xb, xt; const Edge* e; };
    std::vector<Cut> cuts;
    // Slab [ya, yb]: edges sorted by x at ya; when no neighbours swap order at yb, no two edges
    // cross inside the slab and consecutive pairs bound trapezoids. Otherwise split at the
    // crossings of swapped neighbours and recurse.
    auto slab = [&](auto&& self, double ya, double yb, int depth) -> void {
        cuts.clear();
        for (const Edge* e : active)
            if (e->y0 <= ya && e->y1 >= yb) cuts.push_back({x_at(*e, ya), x_at(*e, yb), e});
        std::sort(cuts.begin(), cuts.end(), [](const Cut& p, const Cut& q) { return p.xb != q.xb ? p.xb < q.xb : p.xt < q.xt; });
        std::vector<double> split;
        for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
            const double db = cuts[k].xb - cuts[k + 1].xb, dt = cuts[k].xt - cuts[k + 1].xt;
            if (dt > 0 && db < 0) {
                const double y = ya + (yb - ya) * db / (db - dt);
                if (y > ya && y < yb) split.push_back(y);
            }
        }
        if (!split.empty() && depth < 64) {
            std::sort(split.begin(), split.end());
            split.erase(std::unique(split.begin(), split.end()), split.end());
            double y = ya;
            for (double ys_ : split) { self(self, y, ys_, depth + 1); y = ys_; }
            self(self, y, yb, depth + 1);
            return;
        }
        for (std::size_t k = 0; k + 1 < cuts.size(); k += 2) {
            const Cut &l = cuts[k], &r = cuts[k + 1];
            const Vec3 p0{l.xb, ya, 0}, p1{r.xb, ya, 0}, p2{r.xt, yb, 0}, p3{l.xt, yb, 0};
            if (r.xb > l.xb) tris.insert(tris.end(), {p0, p1, p2});
            if (r.xt > l.xt) tris.insert(tris.end(), {p0, p2, p3});
        }
    };
    for (std::size_t s = 0; s + 1 < ys.size(); ++s) {
        const double ya = ys[s], yb = ys[s + 1];
        while (next < edges.size() && edges[next].y0 <= ya) active.push_back(&edges[next++]);
        active.erase(std::remove_if(active.begin(), active.end(), [&](const Edge* e) { return e->y1 <= ya; }), active.end());
        slab(slab, ya, yb, 0);
    }
}

}  // namespace cadcore::formats::autocad
