// DXF entities against analytic geometry: curves, object coordinate systems, block references,
// colour inheritance, hatches, thickness. Each case is a hand-written DXF (group code / value
// lines) whose expected geometry follows from the DXF reference.
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <set>
#include <string>

#include "cadcore/decode.hpp"
#include "cadcore/stats.hpp"
#include "harness.hpp"

using namespace cadcore;

namespace {

// "code value" pairs separated by '|', e.g. "0|LINE|10|1.0".
std::string dxf(std::initializer_list<std::string> sections) {
    std::string s;
    for (const auto& sec : sections) s += sec;
    s += "0|EOF";
    std::string out;
    for (char c : s) out += c == '|' ? '\n' : c;
    return out + "\n";
}
std::string entities(const std::string& body) { return "0|SECTION|2|ENTITIES|" + body + "0|ENDSEC|"; }
std::string blocks(const std::string& body) { return "0|SECTION|2|BLOCKS|" + body + "0|ENDSEC|"; }
std::string tables(const std::string& body) { return "0|SECTION|2|TABLES|0|TABLE|2|LAYER|" + body + "0|ENDTAB|0|ENDSEC|"; }

Scene load(const std::string& text) { return decode(ByteView(text.data(), text.size()), "t.dxf"); }

struct P { double x, y, z; };
using WorldLine = std::vector<P>;

std::vector<std::pair<WorldLine, Color>> world_lines(const Scene& s) {
    std::vector<std::pair<WorldLine, Color>> out;
    for_each_node(s, [&](const Node& n, const Matrix4& w) {
        for (auto pi : n.polylines) {
            const Polyline& pl = s.polylines[pi];
            WorldLine l;
            for (std::size_t i = 0; i + 2 < pl.points.size(); i += 3) {
                const auto p = apply(w, &pl.points[i]);
                l.push_back({p[0], p[1], p[2]});
            }
            out.push_back({l, pl.color});
        }
    });
    return out;
}

bool near(double a, double b, double tol = 1e-5) { return std::fabs(a - b) <= tol; }

std::set<std::uint32_t> triangle_colors(const Scene& s) {
    std::set<std::uint32_t> out;
    auto add = [&](Color c) { out.insert(std::uint32_t(c.r) << 24 | c.g << 16 | c.b << 8 | c.a); };
    for_each_node(s, [&](const Node& n, const Matrix4&) {
        for (auto mi : n.meshes) {
            const Mesh& m = s.meshes[mi];
            if (m.groups.empty()) add(m.color);
            for (const auto& g : m.groups) add(g.color);
        }
    });
    return out;
}

constexpr double kPi = 3.14159265358979323846;

}  // namespace

TEST_CASE("dxf: circles and arcs lie on their radius and span their angles") {
    const Scene s = load(dxf({entities("0|CIRCLE|8|0|10|1|20|1|30|0|40|2|0|ARC|8|0|10|0|20|0|40|3|50|0|51|90|0|ARC|8|0|10|0|20|0|40|1|50|350|51|10|")}));
    CHECK(s.warnings.empty());
    const auto lines = world_lines(s);
    REQUIRE(lines.size() == 3);
    const auto& c = lines[0].first;
    CHECK(c.size() == 73);  // 72 segments, closed
    for (const auto& p : c) CHECK(near(std::hypot(p.x - 1, p.y - 1), 2, 1e-5));
    CHECK(near(c.front().x, c.back().x) && near(c.front().y, c.back().y));
    const auto& a = lines[1].first;
    CHECK(near(a.front().x, 3) && near(a.front().y, 0));
    CHECK(near(a.back().x, 0) && near(a.back().y, 3));
    const auto& w = lines[2].first;  // 350 -> 10 wraps through 0: 20 degrees
    CHECK(near(std::atan2(w.front().y, w.front().x), -10 * kPi / 180));
    CHECK(near(std::atan2(w.back().y, w.back().x), 10 * kPi / 180));
    for (const auto& p : w) CHECK(p.x > 0.98);
}

TEST_CASE("dxf: polyline bulges are arcs on the correct side") {
    // bulge 1 = counter-clockwise half circle from (0,0) to (2,0): centre (1,0), passes (1,-1)
    const Scene s = load(dxf({entities("0|LWPOLYLINE|8|0|90|2|70|0|10|0|20|0|42|1|10|2|20|0|"
                                       "0|LWPOLYLINE|8|0|90|2|70|0|10|0|20|0|42|-1|10|2|20|0|"
                                       "0|LWPOLYLINE|8|0|90|2|70|0|10|0|20|0|42|0.41421356237|10|1|20|1|")}));
    const auto lines = world_lines(s);
    REQUIRE(lines.size() == 3);
    double ymin = 1e9, ymax = -1e9;
    for (const auto& p : lines[0].first) { ymin = std::min(ymin, p.y); CHECK(near(std::hypot(p.x - 1, p.y), 1)); }
    for (const auto& p : lines[1].first) ymax = std::max(ymax, p.y);
    CHECK(near(ymin, -1));
    CHECK(near(ymax, 1));
    // bulge tan(22.5 deg): quarter circle (0,0) -> (1,1) counter-clockwise, centre (0,1)
    for (const auto& p : lines[2].first) CHECK(near(std::hypot(p.x, p.y - 1), 1));
}

TEST_CASE("dxf: object coordinate system of extrusion (0,0,-1) mirrors x") {
    const Scene s = load(dxf({entities("0|CIRCLE|8|0|10|1|20|2|30|3|40|0.5|210|0|220|0|230|-1|")}));
    const auto lines = world_lines(s);
    REQUIRE(lines.size() == 1);
    double cx = 0, cy = 0, cz = 0;
    for (std::size_t i = 0; i + 1 < lines[0].first.size(); ++i) { cx += lines[0].first[i].x; cy += lines[0].first[i].y; cz += lines[0].first[i].z; }
    const double n = double(lines[0].first.size() - 1);
    CHECK(near(cx / n, -1));  // Ax = Wy x N = (-1, 0, 0)
    CHECK(near(cy / n, 2));
    CHECK(near(cz / n, -3));
}

TEST_CASE("dxf: ellipses and splines") {
    const Scene s = load(dxf({entities("0|ELLIPSE|8|0|10|0|20|0|30|0|11|3|21|0|31|0|40|0.5|41|0|42|6.283185307179586|"
                                       "0|SPLINE|8|0|70|0|71|1|72|5|73|3|40|0|40|0|40|1|40|2|40|2|10|0|20|0|30|0|10|1|20|1|30|0|10|2|20|0|30|0|"
                                       "0|SPLINE|8|0|70|0|71|3|72|0|73|0|74|4|11|0|21|0|31|0|11|1|21|1|31|0|11|2|21|0|31|0|11|3|21|1|31|0|")}));
    CHECK(s.warnings.empty());
    const auto lines = world_lines(s);
    REQUIRE(lines.size() == 3);
    double xmax = 0, ymax = 0;
    for (const auto& p : lines[0].first) { xmax = std::max(xmax, std::fabs(p.x)); ymax = std::max(ymax, std::fabs(p.y)); CHECK(near(p.x * p.x / 9 + p.y * p.y / 2.25, 1)); }
    CHECK(near(xmax, 3));
    CHECK(near(ymax, 1.5, 1.5e-3));  // sampled within the chord tolerance (1e-3 of the radius)
    // degree 1: the control polygon itself
    bool apex = false;
    for (const auto& p : lines[1].first) apex |= near(p.x, 1) && near(p.y, 1);
    CHECK(apex);
    CHECK(near(lines[1].first.front().x, 0));
    CHECK(near(lines[1].first.back().x, 2));
    // fit points only: the interpolating curve passes through every fit point
    for (P q : {P{0, 0, 0}, P{1, 1, 0}, P{2, 0, 0}, P{3, 1, 0}}) {
        double best = 1e9;
        for (const auto& p : lines[2].first) best = std::min(best, std::hypot(p.x - q.x, p.y - q.y));
        CHECK(best < 2e-2);  // sampled: within the chord deviation
    }
    CHECK(near(lines[2].first.front().x, 0) && near(lines[2].first.back().x, 3));
}

TEST_CASE("dxf: INSERT places a block with base point, scale, rotation, array") {
    const std::string blk = blocks("0|BLOCK|8|0|2|B|70|0|10|1|20|0|30|0|0|3DFACE|8|0|10|1|20|0|30|0|11|2|21|0|31|0|12|1|22|1|32|0|13|1|23|1|33|0|0|ENDBLK|8|0|");
    // scale 2, rotation 90: block point (2,0) -> (1,0) relative to base -> scaled (2,0) -> rotated (0,2) -> + (10,0)
    const Scene s = load(dxf({blk, entities("0|INSERT|8|0|2|B|10|10|20|0|30|0|41|2|42|2|43|2|50|90|")}));
    CHECK(s.warnings.empty());
    const SceneStats st = measure(s);
    CHECK(st.triangles == 1);
    CHECK(near(st.area, 2.0));  // area 0.5 scaled by 4
    CHECK(near(st.min[0], 8) && near(st.max[0], 10));
    CHECK(near(st.min[1], 0) && near(st.max[1], 2));
    const Scene arr = load(dxf({blk, entities("0|INSERT|8|0|2|B|10|0|20|0|70|3|71|2|44|5|45|7|")}));
    const SceneStats sa = measure(arr);
    CHECK(sa.triangles == 6);
    CHECK(near(sa.min[0], 0) && near(sa.max[0], 11));
    CHECK(near(sa.min[1], 0) && near(sa.max[1], 8));
}

TEST_CASE("dxf: BYBLOCK colour and layer 0 take the inserting entity's") {
    const std::string lay = tables("0|LAYER|2|0|70|0|62|7|0|LAYER|2|GREEN|70|0|62|3|0|LAYER|2|HIDDEN|70|0|62|-2|");
    const std::string blk = blocks("0|BLOCK|8|0|2|B|70|0|10|0|20|0|30|0|"
                                   "0|3DFACE|8|0|62|0|10|0|20|0|30|0|11|1|21|0|31|0|12|0|22|1|32|0|13|0|23|1|33|0|"
                                   "0|3DFACE|8|0|10|5|20|0|30|0|11|6|21|0|31|0|12|5|22|1|32|0|13|5|23|1|33|0|0|ENDBLK|8|0|");
    const Scene s = load(dxf({lay, blk, entities("0|INSERT|8|GREEN|62|1|2|B|10|0|20|0|0|INSERT|8|0|62|5|2|B|10|0|20|5|0|INSERT|8|HIDDEN|2|B|10|0|20|9|")}));
    CHECK(s.warnings.empty());
    // first insert: BYBLOCK face red, layer-0 BYLAYER face green; second: blue and white
    const std::set<std::uint32_t> expect{0xFF0000FFu, 0x00FF00FFu, 0x0000FFFFu, 0xFFFFFFFFu};
    CHECK(triangle_colors(s) == expect);
    // each block holds two triangles; the insert on the off layer HIDDEN shows nothing, since both
    // of the block's entities are on layer 0 and so take the off layer
    CHECK(measure(s).triangles == 4);
    bool noted = false;
    for (const auto& n : s.notes) noted |= n.find("off or frozen") != std::string::npos;
    CHECK(noted);
}

TEST_CASE("dxf: hatch fills follow the even-odd rule and the hatch style") {
    auto square = [](double x0, double y0, double x1, double y1) {
        return "92|3|72|0|73|1|93|4|10|" + std::to_string(x0) + "|20|" + std::to_string(y0) + "|10|" + std::to_string(x1) + "|20|" + std::to_string(y0) +
               "|10|" + std::to_string(x1) + "|20|" + std::to_string(y1) + "|10|" + std::to_string(x0) + "|20|" + std::to_string(y1) + "|97|0|";
    };
    auto hatch = [&](const std::string& paths, int n, int style) {
        return "0|HATCH|8|0|10|0|20|0|30|0|210|0|220|0|230|1|2|SOLID|70|1|71|0|91|" + std::to_string(n) + "|" + paths + "75|" + std::to_string(style) + "|76|1|98|0|";
    };
    // square with a hole
    CHECK(near(measure(load(dxf({entities(hatch(square(0, 0, 10, 10) + square(2, 2, 4, 4), 2, 0))}))).area, 96));
    // style 2 ignores islands
    CHECK(near(measure(load(dxf({entities(hatch(square(0, 0, 10, 10) + square(2, 2, 4, 4), 2, 2))}))).area, 100));
    // two crossing squares: the overlap (25) is outside under even-odd
    CHECK(near(measure(load(dxf({entities(hatch(square(0, 0, 10, 10) + square(5, 5, 15, 15), 2, 0))}))).area, 150));
    // edge path: a half disc of radius 5 closed by its diameter, written with a clockwise arc
    const Scene hd = load(dxf({entities("0|HATCH|8|0|10|0|20|0|30|0|210|0|220|0|230|1|2|SOLID|70|1|71|0|91|1|92|1|93|2|"
                                        "72|1|10|-5|20|0|11|5|21|0|72|2|10|0|20|0|40|5|50|0|51|180|73|0|97|0|75|0|76|1|98|0|")}));
    CHECK(near(measure(hd).area, kPi * 25 / 2, 0.05));
    // pattern: horizontal lines 1 apart through a 10 x 10 square
    const Scene pat = load(dxf({entities("0|HATCH|8|0|10|0|20|0|30|0|210|0|220|0|230|1|2|USER|70|0|71|0|91|1|" + square(0, 0.5, 10, 10.5) +
                                         "75|0|76|0|52|0|41|1|77|0|78|1|53|0|43|0|44|0|45|0|46|1|79|0|98|0|")}));
    const auto lines = world_lines(pat);
    CHECK(lines.size() == 10);
    for (const auto& l : lines) CHECK(near(std::hypot(l.first.back().x - l.first.front().x, l.first.back().y - l.first.front().y), 10));
}

TEST_CASE("dxf: thickness extrudes curves and solids") {
    // line: a 3 x 2 wall
    CHECK(near(measure(load(dxf({entities("0|LINE|8|0|39|2|10|0|20|0|30|0|11|3|21|0|31|0|")}))).area, 6));
    // solid: a closed unit cube, outward facing for positive and negative thickness
    for (const char* th : {"1", "-1"}) {
        const SceneStats st = measure(load(dxf({entities(std::string("0|SOLID|8|0|39|") + th + "|10|0|20|0|11|1|21|0|12|0|22|1|13|1|23|1|")})));
        CHECK(near(st.area, 6));
        CHECK(near(st.signed_volume, 1));
    }
    // circle r = 1, thickness 2: a closed 72-gon prism
    const double poly = 72 * std::sin(kPi / 72) * std::cos(kPi / 72);  // area of the 72-gon
    const SceneStats cy = measure(load(dxf({entities("0|CIRCLE|8|0|39|2|10|0|20|0|40|1|")})));
    CHECK(near(cy.signed_volume, 2 * poly, 1e-4));
}

TEST_CASE("dxf: polyface mesh, REPEAT array, units, losses reported") {
    const std::string cube =
        "0|POLYLINE|8|0|66|1|70|64|71|8|72|6|"
        "0|VERTEX|8|0|10|0|20|0|30|0|70|192|0|VERTEX|8|0|10|1|20|0|30|0|70|192|0|VERTEX|8|0|10|1|20|1|30|0|70|192|0|VERTEX|8|0|10|0|20|1|30|0|70|192|"
        "0|VERTEX|8|0|10|0|20|0|30|1|70|192|0|VERTEX|8|0|10|1|20|0|30|1|70|192|0|VERTEX|8|0|10|1|20|1|30|1|70|192|0|VERTEX|8|0|10|0|20|1|30|1|70|192|"
        "0|VERTEX|8|0|10|0|20|0|30|0|70|128|71|1|72|4|73|3|74|2|0|VERTEX|8|0|10|0|20|0|30|0|70|128|71|5|72|6|73|7|74|8|"
        "0|VERTEX|8|0|10|0|20|0|30|0|70|128|71|1|72|2|73|6|74|5|0|VERTEX|8|0|10|0|20|0|30|0|70|128|71|2|72|3|73|7|74|6|"
        "0|VERTEX|8|0|10|0|20|0|30|0|70|128|71|3|72|4|73|8|74|7|0|VERTEX|8|0|10|0|20|0|30|0|70|128|71|4|72|1|73|5|74|8|0|SEQEND|8|0|";
    const Scene s = load(dxf({"0|SECTION|2|HEADER|9|$INSUNITS|70|6|0|ENDSEC|", entities(cube)}));
    const SceneStats st = measure(s);
    CHECK(st.triangles == 12);
    CHECK(near(st.signed_volume, 1));
    CHECK(s.unit_to_mm == 1000);
    const Scene r = load(dxf({entities("0|REPEAT|8|0|0|3DFACE|8|0|10|0|20|0|30|0|11|1|21|0|31|0|12|0|22|1|32|0|13|0|23|1|33|0|0|ENDREP|8|0|70|3|71|2|40|2|41|3|")}));
    const SceneStats sr = measure(r);
    CHECK(sr.triangles == 6);
    CHECK(near(sr.max[0], 5) && near(sr.max[1], 4));
    const Scene t = load(dxf({entities("0|TEXT|8|0|10|0|20|0|40|1|1|hello|0|3DSOLID|8|0|1|acis|")}));
    CHECK(t.warnings.size() == 2);
}

TEST_CASE("dxf: broken streams are rejected") {
    auto rejects = [](const std::string& text) {
        try { load(text); } catch (const DecodeError&) { return true; }
        return false;
    };
    CHECK(rejects("0\nSECTION\n2\nENTITIES\n0\nLINE\n10\nabc\n0\nENDSEC\n0\nEOF\n"));       // non-numeric coordinate
    CHECK(rejects("0\nSECTION\n2\nENTITIES\n0\nLINE\n10\n"));                              // value missing
    CHECK(rejects("0\nSECTION\n2\nBLOCKS\n0\nBLOCK\n2\nB\n0\nLINE\n0\nENDSEC\n0\nEOF\n"));  // block without ENDBLK
    CHECK(rejects(dxf({blocks("0|BLOCK|2|B|10|0|20|0|0|INSERT|2|B|10|0|20|0|0|ENDBLK|"), entities("0|INSERT|2|B|10|0|20|0|")})));  // self insert
}

namespace {
// Proxy graphics blob as the 310 hex text of an entity.
struct Blob {
    std::vector<std::uint8_t> b{0, 0, 0, 0, 0, 0, 0, 0};
    int commands = 0;
    void i32(std::int32_t v) { for (int k = 0; k < 4; ++k) b.push_back(std::uint8_t(std::uint32_t(v) >> (8 * k))); }
    void f64(double v) { std::uint64_t u; std::memcpy(&u, &v, 8); for (int k = 0; k < 8; ++k) b.push_back(std::uint8_t(u >> (8 * k))); }
    void cmd(std::int32_t type, const std::vector<double>& pts_or_values, std::vector<std::int32_t> ints_first = {}) {
        i32(std::int32_t(8 + 4 * ints_first.size() + 8 * pts_or_values.size()));
        i32(type);
        for (auto x : ints_first) i32(x);
        for (auto x : pts_or_values) f64(x);
        ++commands;
    }
    std::string groups() {
        std::uint32_t n = std::uint32_t(b.size());
        for (int k = 0; k < 4; ++k) b[std::size_t(k)] = std::uint8_t(n >> (8 * k));
        for (int k = 0; k < 4; ++k) b[std::size_t(4 + k)] = std::uint8_t(std::uint32_t(commands) >> (8 * k));
        static const char* hx = "0123456789ABCDEF";
        std::string hex;
        for (auto c : b) { hex += hx[c >> 4]; hex += hx[c & 15]; }
        return "92|" + std::to_string(n) + "|310|" + hex + "|";
    }
};
}  // namespace

TEST_CASE("dxf: proxy graphics draw unknown entities") {
    Blob p;
    p.cmd(22, {}, {std::int32_t(0xC2FF0000)});                       // true colour red
    p.cmd(6, {0, 0, 0, 1, 0, 0, 1, 1, 0}, {3});                        // polyline, 3 points
    p.cmd(20, {}, {1});                                                // fill on
    p.cmd(29, {1, 0, 0, 10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});     // push: translate x by 10 (row-major)
    p.cmd(7, {0, 0, 0, 2, 0, 0, 0, 2, 0}, {3});                        // filled triangle, area 2
    p.cmd(31, {});
    const Scene s = load(dxf({entities("0|SOMEPROXY|8|0|" + p.groups())}));
    CHECK(s.warnings.empty());
    const auto lines = world_lines(s);
    REQUIRE(lines.size() == 1);
    CHECK(lines[0].first.size() == 3);
    CHECK(lines[0].second == (Color{255, 0, 0, 255}));
    const SceneStats st = measure(s);
    CHECK(near(st.area, 2));
    CHECK(near(st.min[0], 10) && near(st.max[0], 12));
    // a blob that does not frame draws nothing and the entity is reported
    Blob bad;
    bad.cmd(6, {0, 0, 0, 1, 0, 0}, {5});  // says 5 points, holds 2
    const Scene b = load(dxf({entities("0|SOMEPROXY|8|0|" + bad.groups())}));
    CHECK(world_lines(b).empty());
    CHECK(b.warnings.size() == 1);
}

TEST_CASE("dxf: MLINE elements follow the miter offsets; style caps") {
    const std::string objects = "0|SECTION|2|OBJECTS|0|MLINESTYLE|2|TWO|70|272|62|256|71|2|49|0.5|62|1|6|BYLAYER|49|-0.5|62|5|6|BYLAYER|0|ENDSEC|";
    // straight two-vertex mline along x: miter (0,1); offsets +1 / -1
    const std::string ml = "0|MLINE|8|0|2|TWO|40|2|70|1|71|1|72|2|73|2|10|0|20|0|30|0|210|0|220|0|230|1|"
                           "11|0|21|0|31|0|12|1|22|0|32|0|13|0|23|1|33|0|74|2|41|1|41|0|75|0|74|2|41|-1|41|0|75|0|"
                           "11|5|21|0|31|0|12|1|22|0|32|0|13|0|23|1|33|0|74|2|41|1|41|0|75|0|74|2|41|-1|41|0|75|0|";
    const Scene s = load(dxf({entities(ml), objects}));
    CHECK(s.warnings.empty());
    const auto lines = world_lines(s);
    REQUIRE(lines.size() == 4);  // two elements, start and end caps (flags 16 + 256)
    CHECK(near(lines[0].first[0].y, 1) && near(lines[0].first[1].x, 5));
    CHECK(lines[0].second == (Color{255, 0, 0, 255}));
    CHECK(near(lines[1].first[0].y, -1));
    CHECK(lines[1].second == (Color{0, 0, 255, 255}));
    CHECK(near(lines[2].first[0].x, 0) && near(lines[3].first[0].x, 5));
}

TEST_CASE("dxf: leader arrowheads, default and block") {
    // default closed filled arrow, DIMASZ 0.5: length 0.5, base 0.5 / 3, area half their product
    const std::string hdr = "0|SECTION|2|HEADER|9|$DIMASZ|40|0.5|9|$DIMSCALE|40|1|0|ENDSEC|";
    const Scene a = load(dxf({hdr, entities("0|LEADER|8|0|3|STANDARD|76|2|10|0|20|0|30|0|10|4|20|0|30|0|")}));
    CHECK(a.warnings.empty());
    const SceneStats sa = measure(a);
    CHECK(near(sa.area, 0.5 * 0.5 * (0.5 / 3)));
    CHECK(near(sa.min[0], 0) && near(sa.max[0], 0.5));
    // DIMLDRBLK override by handle: a unit block placed at the tip, rotated toward -x, scaled 2
    const std::string tbl = "0|SECTION|2|TABLES|0|TABLE|2|BLOCK_RECORD|0|BLOCK_RECORD|5|AB|2|ARR|0|ENDTAB|0|ENDSEC|";
    const std::string blk = blocks("0|BLOCK|8|0|2|ARR|70|0|10|0|20|0|30|0|0|3DFACE|8|0|10|0|20|0|30|0|11|-1|21|0.5|31|0|12|-1|22|-0.5|32|0|13|-1|23|-0.5|33|0|0|ENDBLK|8|0|");
    const Scene b = load(dxf({tbl, blk, entities("0|LEADER|8|0|3|STANDARD|76|2|10|0|20|0|30|0|10|4|20|0|30|0|1001|ACAD|1000|DSTYLE|1002|{|1070|41|1040|2|1070|341|1005|AB|1002|}|")}));
    CHECK(b.warnings.empty());
    const SceneStats sb = measure(b);
    CHECK(near(sb.area, 4 * 0.5));      // block triangle area 0.5, scaled by 2
    CHECK(near(sb.min[0], 0) && near(sb.max[0], 2));  // the leader points to -x, so the arrow body lies on +x
}

TEST_CASE("dxf: XCLIP clips a block reference to its boundary in block coordinates") {
    // block: a 10 x 10 square outline and a filled 10 x 10 square; reference at (100, 0);
    // boundary in world coordinates (102,2)-(105,20), first matrix = translate(-100, 0)
    const std::string blk = blocks("0|BLOCK|8|0|2|B|70|0|10|0|20|0|30|0|"
                                   "0|LWPOLYLINE|8|0|90|4|70|1|10|0|20|0|10|10|20|0|10|10|20|10|10|0|20|10|"
                                   "0|SOLID|8|0|10|0|20|0|11|10|21|0|12|0|22|10|13|10|23|10|0|ENDBLK|8|0|");
    auto file = [&](int enabled) {
        const std::string m1 = "40|1|40|0|40|0|40|-100|40|0|40|1|40|0|40|0|40|0|40|0|40|1|40|0|";
        const std::string m2 = "40|1|40|0|40|0|40|0|40|0|40|1|40|0|40|0|40|0|40|0|40|1|40|0|";
        const std::string objects = "0|SECTION|2|OBJECTS|0|DICTIONARY|5|D1|3|ACAD_FILTER|360|D2|0|DICTIONARY|5|D2|3|SPATIAL|360|F1|"
                                    "0|SPATIAL_FILTER|5|F1|70|2|10|102|20|2|10|105|20|20|210|0|220|0|230|1|11|0|21|0|31|0|71|" + std::to_string(enabled) +
                                    "|72|0|73|0|" + m1 + m2 + "0|ENDSEC|";
        return dxf({blk, entities("0|INSERT|8|0|360|D1|2|B|10|100|20|0|30|0|"), objects});
    };
    const Scene s = load(file(1));
    CHECK(s.warnings.empty());
    const SceneStats st = measure(s);
    CHECK(near(st.area, 3 * 8));  // the filled square inside x 102..105, y 2..10
    CHECK(near(st.min[0], 102) && near(st.max[0], 105));
    CHECK(near(st.min[1], 2) && near(st.max[1], 10));
    // the outline keeps its two pieces inside the band: bottom and top edges are outside (y 0, y 10 is inside)
    double lo = 1e9, hi = -1e9;
    for (const auto& l : world_lines(s))
        for (const auto& p : l.first) { lo = std::min(lo, p.x); hi = std::max(hi, p.x); }
    CHECK(near(lo, 102) && near(hi, 105));
    // clipping off (71 = 0): the whole block
    const SceneStats off = measure(load(file(0)));
    CHECK(near(off.area, 100));
}
