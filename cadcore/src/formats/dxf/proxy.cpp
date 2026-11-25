#include "formats/dxf/proxy.hpp"

#include <array>
#include <cmath>
#include <cstring>

#include "formats/dxf/colors.hpp"

namespace cadcore::formats::autocad {
namespace {

struct Bad {};  // a command whose data does not match its size

class Data {
public:
    Data(const std::uint8_t* p, std::size_t n) : p_(p), n_(n) {}
    std::size_t left() const { return n_ - at_; }
    std::int32_t i32() { need(4); std::int32_t v; std::memcpy(&v, p_ + at_, 4); at_ += 4; return v; }
    double f64() { need(8); double v; std::memcpy(&v, p_ + at_, 8); at_ += 8; return v; }
    Vec3 pt() { const double x = f64(), y = f64(), z = f64(); return {x, y, z}; }
    std::size_t count() {  // an element count that fits the remaining data
        const std::int32_t v = i32();
        if (v < 0 || std::size_t(v) > left()) throw Bad{};
        return std::size_t(v);
    }

private:
    void need(std::size_t k) const { if (k > left()) throw Bad{}; }
    const std::uint8_t* p_;
    std::size_t n_, at_ = 0;
};

using Mat = std::array<double, 12>;  // row-major 3 x 4
constexpr Mat kId{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
Vec3 xform(const Mat& m, const Vec3& p) {
    return {m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3], m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7], m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11]};
}
Mat multiply(const Mat& a, const Mat& b) {  // a * b
    Mat r{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            double s = j == 3 ? a[std::size_t(i * 4 + 3)] : 0;
            for (int k = 0; k < 3; ++k) s += a[std::size_t(i * 4 + k)] * b[std::size_t(k * 4 + j)];
            r[std::size_t(i * 4 + j)] = s;
        }
    }
    return r;
}

// Circle through three points: centre and radius; false when collinear.
bool circumcircle(const Vec3& a, const Vec3& b, const Vec3& c, Vec3& centre, double& r, Vec3& normal) {
    const Vec3 ab = b - a, ac = c - a;
    const Vec3 n = cross(ab, ac);
    const double nn = dot(n, n);
    if (nn == 0) return false;
    centre = a + (cross(n, ab) * dot(ac, ac) + cross(ac, n) * dot(ab, ab)) * (1.0 / (2 * nn));
    r = length(a - centre);
    normal = normalized(n);
    return true;
}

// Points of the arc around `centre` in the plane of `normal`, from direction `start` sweeping `sweep`.
Polyline3 arc_points(const Vec3& centre, double r, const Vec3& normal, const Vec3& start, double sweep) {
    const Vec3 u = normalized(start), v = cross(normalized(normal), u);
    Polyline3 pts;
    const int n = arc_segments(sweep);
    for (int i = 0; i <= n; ++i) {
        const double t = sweep * i / n;
        pts.push_back(centre + (u * std::cos(t) + v * std::sin(t)) * r);
    }
    return pts;
}

}  // namespace

std::vector<std::uint8_t> proxy_bytes(const Entity& e) {
    std::vector<std::uint8_t> out;
    for (const Group* g = e.begin; g != e.end; ++g) {
        if (g->code != 310) continue;
        if (g->hex) {
            auto nib = [](char c) { return c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1; };
            for (std::size_t i = 0; i + 1 < g->str.size(); i += 2) {
                const int hi = nib(g->str[i]), lo = nib(g->str[i + 1]);
                if (hi < 0 || lo < 0) return {};
                out.push_back(std::uint8_t(hi << 4 | lo));
            }
        } else {
            out.insert(out.end(), g->str.begin(), g->str.end());
        }
    }
    return out;
}

bool draw_proxy(const std::vector<std::uint8_t>& bytes, Color entity_color, Sink& out, std::vector<std::string>& losses) {
    Sink drawn;  // committed only when the whole list frames
    std::vector<std::string> lost;
    try {
        Data head(bytes.data(), bytes.size());
        const std::int32_t total = head.i32(), commands = head.i32();
        if (total < 8 || std::size_t(total) > bytes.size() || commands < 0) return false;
        std::size_t at = 8;
        Color col = entity_color;
        bool fill = false;
        std::vector<Mat> stack{kId};
        auto line = [&](Polyline3 p) {
            for (auto& v : p) v = xform(stack.back(), v);
            drawn.line(std::move(p), col);
        };
        for (std::int32_t c = 0; c < commands; ++c) {
            if (at + 8 > std::size_t(total)) return false;
            std::int32_t size, type;
            std::memcpy(&size, bytes.data() + at, 4);
            std::memcpy(&type, bytes.data() + at + 4, 4);
            if (size < 8 || at + std::size_t(size) > std::size_t(total)) return false;
            Data d(bytes.data() + at + 8, std::size_t(size) - 8);
            bool exact = true;  // the command must use all of its data
            switch (type) {
                case 1: d.pt(); d.pt(); break;  // extents
                case 2: {
                    const Vec3 c0 = d.pt();
                    const double r = d.f64();
                    const Vec3 n = d.pt();
                    const Ocs o(n);
                    line(arc_points(c0, r, n, o.ax, 2 * kPi));
                    break;
                }
                case 3: {
                    const Vec3 a = d.pt(), b = d.pt(), e = d.pt();
                    Vec3 c0, n;
                    double r;
                    if (circumcircle(a, b, e, c0, r, n)) line(arc_points(c0, r, n, a - c0, 2 * kPi));
                    break;
                }
                case 4: {
                    const Vec3 c0 = d.pt();
                    const double r = d.f64();
                    const Vec3 n = d.pt(), sv = d.pt();
                    const double sweep = d.f64();
                    const std::int32_t kind = d.i32();
                    Polyline3 p = arc_points(c0, r, n, sv, sweep);
                    if (kind == 1) { p.insert(p.begin(), c0); p.push_back(c0); }  // sector
                    else if (kind == 2) p.push_back(p.front());                    // chord
                    line(std::move(p));
                    break;
                }
                case 5: {
                    const Vec3 a = d.pt(), b = d.pt(), e = d.pt();
                    d.i32();
                    Vec3 c0, n;
                    double r;
                    if (!circumcircle(a, b, e, c0, r, n)) { line({a, e}); break; }
                    // from a through b to e, counter-clockwise about n
                    auto ang = [&](const Vec3& q) {
                        const Vec3 u = normalized(a - c0), v = cross(n, u), w = q - c0;
                        double t = std::atan2(dot(w, v), dot(w, u));
                        return t < 0 ? t + 2 * kPi : t;
                    };
                    line(arc_points(c0, r, n, a - c0, ang(e)));
                    break;
                }
                case 6:
                case 7:
                case 32: {
                    const std::size_t n = d.count();
                    Polyline3 p;
                    for (std::size_t i = 0; i < n; ++i) p.push_back(d.pt());
                    if (type == 32) d.pt();  // normal
                    if (type == 7 && fill && p.size() >= 3) {
                        for (auto& v : p) v = xform(stack.back(), v);
                        for (std::size_t i = 1; i + 1 < p.size(); ++i) drawn.tri(p[0], p[i], p[i + 1], col);
                        break;
                    }
                    if (type == 7 && p.size() >= 3) p.push_back(p.front());
                    line(std::move(p));
                    break;
                }
                case 9: {
                    const std::size_t nv = d.count();
                    std::vector<Vec3> v;
                    for (std::size_t i = 0; i < nv; ++i) v.push_back(xform(stack.back(), d.pt()));
                    const std::size_t nl = d.count();
                    std::vector<std::int32_t> list;
                    for (std::size_t i = 0; i < nl; ++i) list.push_back(d.i32());
                    for (std::size_t i = 0; i < list.size();) {
                        const std::int32_t k = list[i++];
                        if (k < 0) { lost.push_back("proxy graphics: shell face with a hole"); i += std::size_t(-k); continue; }
                        if (i + std::size_t(k) > list.size()) throw Bad{};
                        for (std::int32_t j = 0; j < k; ++j)
                            if (list[i + std::size_t(j)] < 0 || std::size_t(list[i + std::size_t(j)]) >= v.size()) throw Bad{};
                        for (std::int32_t j = 1; j + 1 < k; ++j)
                            drawn.tri(v[std::size_t(list[i])], v[std::size_t(list[i + std::size_t(j)])], v[std::size_t(list[i + std::size_t(j) + 1])], col);
                        i += std::size_t(k);
                    }
                    exact = false;  // edge / face / vertex attribute data follow, framed by the size
                    break;
                }
                case 12:
                case 13: {
                    const Vec3 a = d.pt(), b = d.pt();
                    const Vec3 dir = normalized(b - a);
                    line({type == 13 ? a : a - dir * 1000.0, a + dir * 1000.0});
                    lost.push_back("proxy graphics: infinite line shown as a 1000-unit segment");
                    break;
                }
                case 44: {
                    const Vec3 c0 = d.pt(), n = d.pt();
                    const double rmaj = d.f64(), rmin = d.f64(), t0 = d.f64(), t1 = d.f64(), tilt = d.f64();
                    d.i32();
                    const Ocs o(n);
                    const Vec3 u = o.ax * std::cos(tilt) + o.ay * std::sin(tilt), v = cross(o.az, u);
                    Polyline3 p;
                    append_ellipse(p, c0, u * rmaj, v * rmin, t0, t1 > t0 ? t1 : t1 + 2 * kPi, false);
                    line(std::move(p));
                    break;
                }
                case 14: {
                    const std::int32_t aci = d.i32();
                    col = aci == 0 || aci == 256 || aci == 257 ? entity_color : aci_color(aci < 0 ? -aci : aci);
                    break;
                }
                case 22: {
                    const std::uint32_t v = std::uint32_t(d.i32());
                    const std::uint32_t method = v >> 24;
                    if (method == 0xC2) col = {std::uint8_t(v >> 16), std::uint8_t(v >> 8), std::uint8_t(v), 255};
                    else if (method == 0xC3) col = aci_color(int(v & 0xFF));
                    else col = entity_color;
                    break;
                }
                case 20: fill = d.i32() != 0; break;
                case 29: {
                    Mat m;
                    double all[16];
                    for (double& x : all) x = d.f64();
                    for (int k = 0; k < 12; ++k) m[std::size_t(k)] = all[k];
                    stack.push_back(multiply(stack.back(), m));
                    break;
                }
                case 31:
                    if (stack.size() > 1) stack.pop_back();
                    break;
                case 16: case 18: case 19: case 23: case 51: d.i32(); break;
                case 24: case 25: d.f64(); break;
                case 34: d.i32(); d.i32(); break;  // 8-byte attribute (always 150 in the corpus), no geometry
                case 11: case 38: lost.push_back("proxy graphics: text is not rendered (no font engine)"); exact = false; break;
                case 33: lost.push_back("proxy graphics: command 33 (length-prefixed bit-packed data) not decoded"); exact = false; break;
                default: lost.push_back("proxy graphics: command " + std::to_string(type) + " not decoded"); exact = false; break;
            }
            if (exact && d.left() != 0) return false;
            at += std::size_t(size);
        }
    } catch (const Bad&) {
        return false;
    }
    for (auto& l : drawn.lines) out.line(std::move(l.pts), l.color);
    for (std::size_t i = 0; i < drawn.tri_colors.size(); ++i) out.tri(drawn.tri_pts[3 * i], drawn.tri_pts[3 * i + 1], drawn.tri_pts[3 * i + 2], drawn.tri_colors[i]);
    for (std::size_t i = 0; i < drawn.points.size(); ++i) out.point(drawn.points[i], drawn.point_colors[i]);
    losses.insert(losses.end(), lost.begin(), lost.end());
    return true;
}

}  // namespace cadcore::formats::autocad
