#include "formats/catia/cgm_model.hpp"

#include <algorithm>
#include <cstring>
#include <deque>
#include <functional>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include "cadcore/error.hpp"
#include "formats/catia/cgm_compact.hpp"

namespace cadcore::formats::catia::cgm {
namespace {

[[noreturn]] void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CGM record: " + what, at); }

struct Obj {
    std::uint16_t type = 0;
    std::vector<std::uint32_t> refs;
    ByteView p;  // payload after the reference list
    std::size_t at = 0;
    std::string error;  // compact record whose body did not decode (raised when the object is used)
};

// Sequential reader over a payload in the stream's byte order.
struct Rd {
    ByteView b;
    bool big;
    std::size_t k = 0, at = 0;
    void need(std::size_t n) const { if (k + n > b.size) fail(ErrorCode::Truncated, "CGM payload ends early", at + k); }
    double f64() {
        need(8);
        std::uint8_t t[8];
        for (int i = 0; i < 8; ++i) t[i] = b[k + (big ? 7 - i : i)];
        double d;
        std::memcpy(&d, t, 8);
        k += 8;
        return d;
    }
    std::int32_t i32() { need(4); auto v = std::int32_t(load_u32(b.data + k, big)); k += 4; return v; }
    std::int16_t i16() { need(2); auto v = std::int16_t(load_u16(b.data + k, big)); k += 2; return v; }
    std::uint8_t u8() { need(1); return b[k++]; }
    std::uint32_t var() { return varint(b, k); }
    // A count: a varint, or a raw u32 in streams of the oldest releases (low byte < 0x80, then 00 00 00).
    std::uint32_t count() {
        if (k + 4 <= b.size && b[k] < 0x80 && b[k + 1] == 0 && b[k + 2] == 0 && b[k + 3] == 0) return std::uint32_t(i32());
        return var();
    }
    std::size_t left() const { return b.size - k; }
};

struct KnotVector { int degree = 0; int first = 0; std::vector<double> knots; std::vector<int> mults; std::size_t poles() const { std::size_t n = 0; for (int m : mults) n += std::size_t(m); return n - std::size_t(degree) - 1; } };

KnotVector knot_vector(Rd& r) {
    KnotVector kv;
    kv.degree = r.i32();
    r.i32();
    kv.first = r.i32();
    const int n = r.i32();
    r.i32();
    if (n < 2 || n > 100000) bad("knot count out of range", r.at + r.k);
    for (int i = 0; i < n; ++i) kv.knots.push_back(r.f64());
    for (int i = 0; i < n; ++i) kv.mults.push_back(r.i32());
    return kv;
}

class Decoder {
public:
    Decoder(const Stream& s, std::vector<std::vector<std::uint32_t>> visible, const CompactLayout* layout, const DecodeOptions& opt) : s_(s), visible_(std::move(visible)), opt_(opt) {
        if (layout) {
            // compact bodies in canonical (varint-dialect) form
            auto canon = canonical_records(s, *layout);
            for (std::size_t i = 0; i < s.records.size(); ++i) {
                const Record& rec = s.records[i];
                Obj o;
                o.type = rec.type;
                o.at = rec.offset;
                if (canon[i].error.empty()) {
                    o.refs = std::move(canon[i].c.refs);
                    payloads_.push_back(std::move(canon[i].c.payload));
                    o.p = ByteView(payloads_.back().data(), payloads_.back().size());
                } else {
                    o.error = std::move(canon[i].error);
                }
                objs_.emplace(rec.id, std::move(o));
            }
            return;
        }
        for (const Record& rec : s.records) {
            Obj o;
            o.type = rec.type;
            o.at = rec.offset;
            std::size_t k = 0;
            const std::uint32_t n = varint(rec.body, k);
            for (std::uint32_t i = 0; i < n; ++i) o.refs.push_back(varint(rec.body, k));
            o.p = ByteView(rec.body.data + k, rec.body.size - k);
            objs_.emplace(rec.id, std::move(o));
        }
    }

    Ir run() {
        // shells of volumes (525 -> 513 -> 520) bound solids; other shells are sheets
        for (const Record& rec : s_.records)
            if (rec.type == 513 && usable(rec.id))
                for (std::uint32_t r : obj(rec.id).refs) solid_.insert(r);
        for (const Record& rec : s_.records)
            if (rec.type == 520) shell(rec.id);
        // bodies: 500 -> 520 (sheets), 500 -> 525 -> 513 -> 520 (solids); versions of a body share
        // their unchanged shells, so a shell can belong to several bodies
        for (const Record& rec : s_.records) {
            if (rec.type != 500 || !usable(rec.id)) continue;
            Ir::Body bd;
            bd.id = rec.id;
            std::vector<std::uint32_t> todo = obj(rec.id).refs;
            for (std::size_t i = 0; i < todo.size() && i < 100000; ++i) {
                if (!usable(todo[i])) continue;
                const Obj& o = obj(todo[i]);
                if (o.type == 525 || o.type == 513) todo.insert(todo.end(), o.refs.begin(), o.refs.end());
                if (o.type == 520)
                    if (auto it = shell_.find(todo[i]); it != shell_.end() && std::find(bd.shells.begin(), bd.shells.end(), it->second) == bd.shells.end()) bd.shells.push_back(it->second);
            }
            if (std::find(s_.roots.begin(), s_.roots.end(), rec.id) != s_.roots.end()) bd.flags |= 1;
            for (std::size_t j = 0; j < visible_.size() && j < 30; ++j)
                if (std::find(visible_[j].begin(), visible_[j].end(), rec.id) != visible_[j].end()) bd.flags |= 2u << j;
            out_.bodies.push_back(std::move(bd));
        }
        out_.vertices = std::uint32_t(vertex_.size());
        return std::move(out_);
    }

private:
    const Stream& s_;
    std::vector<std::vector<std::uint32_t>> visible_;
    DecodeOptions opt_;
    std::unordered_map<std::uint32_t, Obj> objs_;
    std::deque<std::vector<std::uint8_t>> payloads_;  // canonical payloads of compact records
    std::unordered_map<std::uint32_t, std::uint32_t> surface_, pcurve_, edge_, vertex_;
    std::unordered_map<std::uint32_t, std::optional<std::uint32_t>> face_;
    std::unordered_set<std::uint32_t> solid_;
    std::unordered_map<std::uint32_t, std::uint32_t> shell_;
    Ir out_;

    const Obj& obj(std::uint32_t id) const {
        auto it = objs_.find(id);
        if (it == objs_.end()) bad("reference to a missing object " + std::to_string(id), 0);
        if (!it->second.error.empty()) throw DecodeError(ErrorCode::Unsupported, it->second.error, it->second.at);
        return it->second;
    }
    Rd rd(const Obj& o) const { return Rd{o.p, s_.big, 0, o.at}; }
    bool usable(std::uint32_t id) const { auto it = objs_.find(id); return it != objs_.end() && it->second.error.empty(); }
    static std::uint32_t first(const Obj& o) {
        if (o.refs.empty()) bad("record without its first reference", o.at);
        return o.refs[0];
    }

    void pending(const Obj& o, const char* what) const {
        if (!opt_.pending_kinds) throw DecodeError(ErrorCode::Unsupported, "surface type " + std::to_string(o.type) + " (" + what + ") awaits kernel support", o.at);
    }

    // A revolution's profile (cgm_model.hpp SurfRevolution): edge curves stand for their first curve.
    void profile(std::uint32_t c, std::vector<double>& p) {
        for (int depth = 0; is_edge_curve(obj(c).type); ++depth) {
            if (depth > 8) bad("edge curves nest too deep", obj(c).at);
            c = first(obj(c));
        }
        const Obj& o = obj(c);
        Rd r = rd(o);
        switch (o.type) {
            case 150: case 151: case 160:
                p.push_back(1);
                p.push_back(double(pcurve(c)));
                return;
            case 100:
                p.push_back(2);
                for (int i = 0; i < 6; ++i) p.push_back(r.f64());
                return;
            case 130: {
                const KnotVector kv = knot_vector(r);
                const bool rat = r.i16() != 0;
                const std::size_t np = kv.poles();
                p.push_back(3);
                p.push_back(double(kv.degree));
                std::vector<double> knots;
                for (std::size_t i = 0; i < kv.knots.size(); ++i)
                    for (int m = 0; m < kv.mults[i]; ++m) knots.push_back(kv.knots[i]);
                p.push_back(double(knots.size()));
                p.insert(p.end(), knots.begin(), knots.end());
                p.push_back(double(np));
                p.push_back(rat ? 1.0 : 0.0);
                for (std::size_t i = 0; i < 3 * np; ++i) p.push_back(r.f64());
                if (rat) for (std::size_t i = 0; i < np; ++i) p.push_back(r.f64());
                return;
            }
            default: throw DecodeError(ErrorCode::Unsupported, "revolution profile of class " + std::to_string(o.type), o.at);
        }
    }

    std::uint32_t surface(std::uint32_t id) {
        if (auto it = surface_.find(id); it != surface_.end()) return it->second;
        const Obj& o = obj(id);
        if (o.type == 212 || o.type == 241 || o.type == 242 || o.type == 244) {
            const std::uint32_t r = surface(first(o));
            surface_[id] = r;
            return r;
        }
        Ir::Surface sf;
        sf.id = id;
        Rd r = rd(o);
        auto vec = [&](int n) { for (int i = 0; i < n; ++i) sf.p.push_back(r.f64()); };
        switch (o.type) {
            case 200: {
                sf.kind = SurfPlane;
                if (o.p.size == 120) {
                    vec(9);
                } else if (o.p.size % 8 == 1) {
                    vec(3);
                    const std::uint8_t mask = r.u8();
                    for (int j = 0; j < 6; ++j) sf.p.push_back(mask >> j & 1 ? 0.0 : r.f64());
                } else {
                    bad("plane payload length", o.at);
                }
                break;
            }
            case 201: {
                if (o.p.size / 8 < 15) bad("cylinder payload length", o.at);
                sf.kind = SurfCylinder;
                vec(9);
                const double rad = r.f64();
                r.k += 4 * 8;
                const double su = r.f64();
                sf.p.push_back(rad);
                sf.p.push_back(rad * su);
                break;
            }
            case 202: {
                if (o.p.size != 184) bad("cone payload length", o.at);
                sf.kind = SurfCone;
                vec(14);          // o X Y Z, half-angle, radius at v = 0
                r.k += 4 * 8;     // u0 u1 v0 v1
                sf.p.push_back(r.f64());  // su
                break;
            }
            case 210: {
                // profile: a p-curve, or an edge curve whose first curve is one (same parameter)
                std::uint32_t prof = first(o);
                if (is_edge_curve(obj(prof).type)) prof = first(obj(prof));
                const double pc = double(pcurve(prof));
                sf.kind = SurfExtrusion;
                sf.p.push_back(pc);
                vec(3);
                break;
            }
            case 214: {
                if (o.refs.size() != 2) bad("offset surface without carrier and base", o.at);
                const double base = double(surface(o.refs[1]));
                sf.kind = SurfOffset;
                sf.p = {base, r.f64()};
                break;
            }
            case 203: sf.kind = SurfSphere; vec(12); sf.p.push_back(r.f64()); break;
            case 204: sf.kind = SurfTorus; vec(14); break;
            case 211: {
                pending(o, "revolution");
                sf.kind = SurfRevolution;
                vec(12);
                r.k += 4 * 8 + 8;  // box, two ints
                sf.p.push_back(r.f64());  // k
                profile(first(o), sf.p);
                break;
            }
            case 216: {
                pending(o, "fillet");
                sf.kind = SurfFillet;
                const std::int32_t n = r.i32();
                const KnotVector kv = knot_vector(r);
                if (kv.degree != 5 || kv.knots.size() != std::size_t(n)) bad("fillet knot vector", o.at);
                sf.p.push_back(double(n));
                sf.p.insert(sf.p.end(), kv.knots.begin(), kv.knots.end());
                vec(30 * n + 4);
                break;
            }
            case 231: {
                sf.kind = SurfNurbs;
                const KnotVector u = knot_vector(r), v = knot_vector(r);
                const bool rat = r.i16() != 0;
                const std::size_t nu = u.poles(), nv = v.poles();
                sf.p = {double(u.degree), double(v.degree), double(nu), double(nv), rat ? 1.0 : 0.0};
                for (const KnotVector* kv : {&u, &v})
                    for (std::size_t i = 0; i < kv->knots.size(); ++i)
                        for (int m = 0; m < kv->mults[i]; ++m) sf.p.push_back(kv->knots[i]);
                vec(int(3 * nu * nv));
                if (rat) vec(int(nu * nv));
                break;
            }
            default: throw DecodeError(ErrorCode::Unsupported, "surface type " + std::to_string(o.type), o.at);
        }
        const auto i = std::uint32_t(out_.surfaces.size());
        out_.surfaces.push_back(std::move(sf));
        surface_[id] = i;
        return i;
    }

    std::uint32_t pcurve(std::uint32_t id) {
        if (auto it = pcurve_.find(id); it != pcurve_.end()) return it->second;
        const Obj& o = obj(id);
        Ir::PCurve pc;
        pc.id = id;
        pc.surface = surface(first(o));
        Rd r = rd(o);
        switch (o.type) {
            case 150:
                pc.kind = PcLine;
                if (o.p.size == 48) {
                    for (int i = 0; i < 4; ++i) pc.p.push_back(r.f64());
                } else if (o.p.size == 64) {
                    pc.p.push_back(r.f64());
                    pc.p.push_back(r.f64());
                    r.f64();
                    pc.p.push_back(r.f64());
                    pc.p.push_back(r.f64());
                } else {
                    bad("line p-curve payload length", o.at);
                }
                break;
            case 151: {
                pc.kind = PcCircle;
                pc.p.push_back(r.f64());
                pc.p.push_back(r.f64());
                r.i32();
                r.i32();
                const double rad = r.f64();
                r.f64();
                r.f64();
                pc.p.push_back(rad);
                pc.p.push_back(r.f64());
                break;
            }
            case 160: {
                pc.kind = PcHermite;
                const KnotVector kv = knot_vector(r);
                const auto n = std::size_t(r.i32());
                if (n != kv.knots.size()) bad("Hermite value count differs from the knot count", r.at + r.k);
                pc.p.push_back(double(n));
                pc.p.insert(pc.p.end(), kv.knots.begin(), kv.knots.end());
                for (std::size_t i = 0; i < 4 * n; ++i) pc.p.push_back(r.f64());
                const bool has2 = r.i16() != 0;
                pc.p.push_back(has2 ? 1.0 : 0.0);
                if (has2) for (std::size_t i = 0; i < 2 * n; ++i) pc.p.push_back(r.f64());
                break;
            }
            default: throw DecodeError(ErrorCode::Unsupported, "p-curve type " + std::to_string(o.type), o.at);
        }
        const auto i = std::uint32_t(out_.pcurves.size());
        out_.pcurves.push_back(std::move(pc));
        pcurve_[id] = i;
        return i;
    }

    // Parameter maps of an edge curve: f[i](q) is the parameter on refs[i] for q on refs[0].
    std::vector<std::function<double(double)>> edge_maps(const Obj& ec) {
        Rd r = rd(ec);
        std::vector<std::function<double(double)>> f{[](double q) { return q; }};
        const std::uint32_t n = r.count();
        bool tail = false;
        for (std::uint32_t j = 0; j < n; ++j) {
            const int map = r.i16();
            const double a0 = r.f64(), a1 = r.f64();
            r.f64();
            if (map == 0) continue;
            if (map == 1) { r.f64(); r.f64(); r.f64(); f.push_back([](double q) { return q; }); continue; }
            const double b0 = r.f64(), b1 = r.f64();
            if (map == 2) {
                f.push_back([=](double q) { return a1 != a0 ? b0 + (q - a0) * (b1 - b0) / (a1 - a0) : b0; });
            } else if (map == 3) {
                const double s = r.f64(), off = r.f64();
                f.push_back([=](double q) { return (q - off) / s; });
            } else if (map == 5) {
                const std::uint32_t m = r.count();
                std::vector<double> A(m), B(m);
                for (auto& x : A) x = r.f64();
                for (auto& x : B) x = r.f64();
                f.push_back([A, B](double q) {
                    std::size_t i = std::size_t(std::upper_bound(A.begin(), A.end(), q) - A.begin());
                    i = std::min(std::max<std::size_t>(i, 1), A.size() - 1) - 1;
                    const double h = A[i + 1] - A[i];
                    return h != 0 ? B[i] + (q - A[i]) / h * (B[i + 1] - B[i]) : B[i];
                });
            } else {
                bad("edge-curve map " + std::to_string(map), r.at + r.k);
            }
            tail = true;
        }
        if (tail) r.k += 24;
        if (r.left() != 4) bad("edge-curve payload length", ec.at);
        // some compact edge curves list their curves twice (A B A B) with the maps of the first copy
        const std::size_t h = ec.refs.size() / 2;
        if (f.size() == h && ec.refs.size() == 2 * h && std::equal(ec.refs.begin(), ec.refs.begin() + std::ptrdiff_t(h), ec.refs.begin() + std::ptrdiff_t(h)))
            for (std::size_t i = 0; i < h; ++i) f.push_back(f[i]);
        if (f.size() != ec.refs.size()) bad("edge-curve maps differ from its curves", ec.at);
        return f;
    }

    static bool is_edge_curve(std::uint16_t t) { return t == 191 || t == 192 || t == 193; }

    // Parameter on `target` of the point with parameter q on curve `c`, following edge curves
    // (which may nest: an edge curve can list another edge curve among its curves).
    std::optional<double> through(std::uint32_t c, double q, std::uint32_t target, int depth) {
        if (c == target) return q;
        const Obj& e = obj(c);
        if (!is_edge_curve(e.type) || depth > 8) return std::nullopt;
        const auto maps = edge_maps(e);
        for (std::size_t j = 0; j < e.refs.size(); ++j)
            if (auto r = through(e.refs[j], maps[j](q), target, depth + 1)) return r;
        return std::nullopt;
    }

    // Start parameter of segment `span` (1-based) of curve `c`: a point's parameter is an offset
    // within that segment. Hermite p-curves have one segment per knot interval; an edge curve's
    // parameter is its first curve's. Single-segment curves only use span 1.
    double span_start(std::uint32_t c, std::int32_t span, std::size_t at) {
        if (span == 1) return 0;
        const Obj& o = obj(c);
        if (is_edge_curve(o.type)) return span_start(first(o), span, at);
        if (o.type == 160) {
            Rd r = rd(o);
            const KnotVector kv = knot_vector(r);
            const std::int32_t j = span + kv.first;  // segment numbers are shifted by the knot vector's offset field (<= 0)
            if (j < 1 || std::size_t(j) >= kv.knots.size()) bad("point segment " + std::to_string(span) + " outside Hermite curve " + std::to_string(c) + " (first " + std::to_string(kv.first) + ", " + std::to_string(kv.knots.size()) + " knots)", at);
            return kv.knots[std::size_t(j) - 1];
        }
        throw DecodeError(ErrorCode::Unsupported, "point segment " + std::to_string(span) + " on curve type " + std::to_string(o.type), at);
    }

    // Parameter of a point (type 14: parameters on several curves) on p-curve `pc`.
    double point_param(std::uint32_t point, std::uint32_t pc) {
        const Obj& pt = obj(point);
        Rd r = rd(pt);
        const std::uint32_t n = r.count();
        if (n != pt.refs.size()) bad("point parameter count differs from its curves", pt.at);
        std::vector<double> q(n);
        for (std::uint32_t i = 0; i < n; ++i) { q[i] = r.f64(); q[i] += span_start(pt.refs[i], r.i32(), pt.at); }
        if (r.left() != 0) bad("point payload length", pt.at);
        for (std::uint32_t i = 0; i < n; ++i)
            if (auto t = through(pt.refs[i], q[i], pc, 0)) return *t;
        throw DecodeError(ErrorCode::Unsupported, "p-curve " + std::to_string(pc) + " not reachable from point " + std::to_string(point), pt.at);
    }

    std::uint32_t vertex(std::uint32_t id) {
        auto [it, fresh] = vertex_.emplace(id, std::uint32_t(vertex_.size()));
        return it->second;
    }

    std::uint32_t edge(std::uint32_t id) {
        if (auto it = edge_.find(id); it != edge_.end()) return it->second;
        const Obj& e = obj(id);
        if (e.refs.size() != 5 && e.refs.size() != 4) bad("edge without its vertices and points", e.at);
        const std::size_t b = e.refs.size() - 4;  // [edge curve,] v0 v1 p0 p1
        Ir::Edge ed{vertex(e.refs[b]), vertex(e.refs[b + 1]), id};
        const auto i = std::uint32_t(out_.edges.size());
        out_.edges.push_back(ed);
        edge_[id] = i;
        return i;
    }

    // Appends the loop's coedges (throws when the loop cannot be decoded).
    void loop(std::uint32_t id) {
        const Obj& l = obj(id);
        // V5R8 and later end the references with the surface; the oldest releases do not.
        const bool with_surface = l.refs.size() % 2 == 1;
        const std::size_t n = l.refs.size() / 2;
        if (n == 0) bad("loop without coedges", l.at);
        Rd r = rd(l);
        if (r.var() != n) bad("loop count differs from its coedges", l.at);
        const std::size_t flags = 2 * (3 * n + (with_surface ? 4 : 3));
        // A non-zero closing i16 (3, or 1 in compact streams) announces a 73-byte block (not decoded;
        // it holds no topology).
        const bool block = with_surface && r.left() == flags + 73 && load_u16(r.b.data + r.k + flags - 2, s_.big) != 0;
        if (r.left() != flags && !block) throw DecodeError(ErrorCode::Unsupported, "loop flag layout of " + std::to_string(r.left()) + " bytes", l.at);
        r.i16();
        r.i16();
        std::vector<Ir::Coedge> ces;
        for (std::size_t j = 0; j < n; ++j) {
            r.i16();
            const int sense = r.i16();
            r.i16();
            const std::uint32_t pcid = l.refs[2 * j], eid = l.refs[2 * j + 1];
            const Obj& e = obj(eid);
            Ir::Coedge c;
            c.pcurve = pcurve(pcid);
            c.edge = edge(eid);
            c.sense = sense < 0 ? -1 : 1;
            if (e.refs.size() < 4) bad("edge without its points", e.at);
            c.t0 = point_param(e.refs[e.refs.size() - 2], pcid);
            c.t1 = point_param(e.refs[e.refs.size() - 1], pcid);
            ces.push_back(c);
        }
        out_.loops.push_back({std::uint32_t(out_.coedges.size()), std::uint32_t(n)});
        out_.coedges.insert(out_.coedges.end(), ces.begin(), ces.end());
    }

    std::optional<std::uint32_t> face(std::uint32_t id) {
        if (auto it = face_.find(id); it != face_.end()) return it->second;
        const auto r = decode_face(id);
        face_[id] = r;
        return r;
    }

    std::optional<std::uint32_t> decode_face(std::uint32_t id) {
        const Obj& f = obj(id);
        const auto nc = out_.coedges.size(), nl = out_.loops.size();
        try {
            if (f.refs.size() < 2) bad("face without surface and loop", f.at);
            Rd r = rd(f);
            Ir::Face fc;
            fc.id = id;
            fc.surface = surface(f.refs[0]);
            fc.sense = r.i16() < 0 ? -1 : 1;
            fc.first_loop = std::uint32_t(out_.loops.size());
            for (std::size_t i = 1; i < f.refs.size(); ++i) loop(f.refs[i]);
            fc.loops = std::uint32_t(out_.loops.size()) - fc.first_loop;
            out_.faces.push_back(fc);
            return std::uint32_t(out_.faces.size() - 1);
        } catch (const DecodeError& e) {
            out_.coedges.resize(nc);
            out_.loops.resize(nl);
            out_.skipped.push_back({id, 512, e.what()});
            return std::nullopt;
        }
    }

    void shell(std::uint32_t id) {
        Ir::Shell sh;
        sh.id = id;
        sh.solid = solid_.count(id) != 0;
        try {
            const Obj& s = obj(id);
            Rd r = rd(s);
            const std::uint32_t n = r.var();
            if (n != s.refs.size()) bad("shell count differs from its faces", s.at);
            r.i16();
            std::vector<int> senses(n);
            for (auto& x : senses) x = r.i16();
            r.i16();
            if (r.left() != 0) bad("shell payload length", s.at);
            for (std::uint32_t i = 0; i < n; ++i)
                if (auto f = face(s.refs[i])) { sh.faces.push_back(*f); sh.senses.push_back(senses[i] < 0 ? -1 : 1); }
        } catch (const DecodeError& e) {
            out_.skipped.push_back({id, 520, e.what()});
            return;
        }
        shell_[id] = std::uint32_t(out_.shells.size());
        out_.shells.push_back(std::move(sh));
    }
};

struct Out {
    std::vector<std::uint8_t> b;
    void u32(std::uint32_t v) { for (int i = 0; i < 4; ++i) b.push_back(std::uint8_t(v >> (8 * i))); }
    void i32(std::int32_t v) { u32(std::uint32_t(v)); }
    void f64(double d) { std::uint64_t v; std::memcpy(&v, &d, 8); for (int i = 0; i < 8; ++i) b.push_back(std::uint8_t(v >> (8 * i))); }
};

}  // namespace

Ir decode(const Stream& stream, const std::vector<std::vector<std::uint32_t>>& visible, const DecodeOptions& options) {
    if (stream.varint_bodies) return Decoder(stream, visible, nullptr, options).run();
    const auto layout = compact_layout(stream);
    if (!layout) throw DecodeError(ErrorCode::Unsupported, std::string("CGM ") + dialect_name(stream.dialect) + " dialect bodies of level " + std::to_string(stream.level) + " (" + stream.version + ") not decoded", 0);
    return Decoder(stream, visible, &*layout, options).run();
}

std::vector<std::uint8_t> serialize(const Ir& ir) {
    Out o;
    const char magic[8] = {'C', 'G', 'M', 'I', 'R', 0, 0, 2};
    o.b.assign(magic, magic + 8);
    for (std::size_t n : {ir.surfaces.size(), ir.pcurves.size(), std::size_t(ir.vertices), ir.edges.size(), ir.coedges.size(), ir.loops.size(), ir.faces.size(), ir.shells.size(), ir.bodies.size(), ir.skipped.size()}) o.u32(std::uint32_t(n));
    for (const auto& s : ir.surfaces) { o.u32(s.kind); o.u32(s.id); o.u32(std::uint32_t(s.p.size())); for (double d : s.p) o.f64(d); }
    for (const auto& c : ir.pcurves) { o.u32(c.kind); o.u32(c.surface); o.u32(c.id); o.u32(std::uint32_t(c.p.size())); for (double d : c.p) o.f64(d); }
    for (const auto& e : ir.edges) { o.u32(e.v0); o.u32(e.v1); o.u32(e.id); }
    for (const auto& c : ir.coedges) { o.u32(c.pcurve); o.u32(c.edge); o.i32(c.sense); o.u32(0); o.f64(c.t0); o.f64(c.t1); }
    for (const auto& l : ir.loops) { o.u32(l.first); o.u32(l.count); }
    for (const auto& f : ir.faces) { o.u32(f.surface); o.i32(f.sense); o.u32(f.first_loop); o.u32(f.loops); o.u32(f.id); }
    for (const auto& s : ir.shells) { o.u32(s.id); o.u32(s.solid ? 1u : 0u); o.u32(std::uint32_t(s.faces.size())); for (std::size_t i = 0; i < s.faces.size(); ++i) { o.u32(s.faces[i]); o.i32(s.senses[i]); } }
    for (const auto& b : ir.bodies) { o.u32(b.id); o.u32(b.flags); o.u32(std::uint32_t(b.shells.size())); for (auto i : b.shells) o.u32(i); }
    for (const auto& k : ir.skipped) { o.u32(k.id); o.u32(k.type); o.u32(std::uint32_t(k.why.size())); o.b.insert(o.b.end(), k.why.begin(), k.why.end()); }
    return o.b;
}

}  // namespace cadcore::formats::catia::cgm
