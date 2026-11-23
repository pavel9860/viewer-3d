// AutoCAD DXF — ASCII and binary, R2.10 (AC2.10) to 2018 (AC1032).
//
// Model space (the ENTITIES section) becomes the root node. Each block is decoded once per
// display context and shared by every INSERT that places it in that context; an INSERT (and
// each cell of a MINSERT array) is a node whose transform maps block coordinates to the
// parent's: OCS(210) * T(10) * Rz(50) * T(array cell) * S(41, 42, 43) * T(-block base).
// The display context is what block entities inherit from the INSERT: the BYBLOCK colour and
// the layer that replaces layer "0".
//
// Dimensions and tables draw through their anonymous blocks (*D..., *T...), as AutoCAD does.
// Text, ACIS solids and raster images need engines this build does not have; every such
// entity is reported as a warning with its count, as is any content of a drawn entity that is
// not shown (e.g. a leader arrowhead).
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <initializer_list>
#include <iterator>

#include "formats/dxf/colors.hpp"
#include "formats/dxf/document.hpp"
#include "formats/dxf/draw.hpp"
#include "formats/dxf/proxy.hpp"
#include "formats/dxf/clip.hpp"
#include "cadcore/traverse.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

using namespace autocad;

Node named_node(std::string name) {
    Node n;
    n.name = std::move(name);
    return n;
}

int sniff(ByteView d, std::string_view ext) {
    if (is_binary_dxf(d) || looks_like_ascii_dxf(d)) return 2;
    return ext == "dxf" ? 1 : 0;
}

// Millimetres per drawing unit of $INSUNITS (0 = unitless and unknown codes: 1).
double unit_to_mm(int insunits) {
    static const double mm[] = {1, 25.4, 304.8, 1609344.0, 1, 10, 1000, 1e6, 25.4e-6, 0.0254, 914.4, 1e-7, 1e-6, 1e-3, 100, 1e4, 1e5, 1e12,
                                1.495978707e14, 9.4607304725808e18, 3.0856775814913673e19, 1200.0 / 3937.0 * 1000.0};
    return insunits >= 0 && insunits < int(std::size(mm)) ? mm[insunits] : 1.0;
}

// What block entities inherit from the INSERT that places them.
struct Inherit {
    Color byblock{255, 255, 255, 255};  // top level: BYBLOCK shows as colour 7
    std::string layer;                  // upper-case layer that replaces "0"; empty at top level
};

enum class Kind { Draw, Insert, BlockRef, Repeat, Text, Brep, NonGraphic, Lost };

Kind classify(std::string_view t) {
    static const std::map<std::string_view, Kind> kinds = [] {
        std::map<std::string_view, Kind> m;
        auto add = [&](Kind k, std::initializer_list<std::string_view> names) { for (auto n : names) m[n] = k; };
        add(Kind::Draw, {"MLINE", "LINE", "3DLINE", "POINT", "CIRCLE", "ARC", "ELLIPSE", "LWPOLYLINE", "POLYLINE", "3DFACE", "SOLID", "TRACE", "SPLINE",
                         "HELIX", "MESH", "HATCH", "RAY", "XLINE", "LEADER"});
        add(Kind::Insert, {"INSERT", "ACAD_TABLE"});
        add(Kind::BlockRef, {"DIMENSION", "ARC_DIMENSION", "LARGE_RADIAL_DIMENSION"});
        add(Kind::Repeat, {"REPEAT"});
        add(Kind::Text, {"TEXT", "MTEXT", "ATTRIB", "ATTDEF", "SHAPE", "ARCALIGNEDTEXT", "RTEXT", "TOLERANCE"});
        add(Kind::Brep, {"3DSOLID", "REGION", "BODY", "SURFACE", "PLANESURFACE", "EXTRUDEDSURFACE", "LOFTEDSURFACE", "REVOLVEDSURFACE",
                         "SWEPTSURFACE", "NURBSURFACE"});
        add(Kind::NonGraphic, {"VIEWPORT", "LIGHT", "SUN", "SECTIONOBJECT", "CAMERA"});
        return m;
    }();
    const auto it = kinds.find(t);
    return it == kinds.end() ? Kind::Lost : it->second;
}

class Builder {
public:
    Builder(const Document& doc, Scene& scene) : doc_(doc), scene_(scene) {
        dctx_.has_extents = doc.has_extents;
        dctx_.extmin = doc.extmin;
        dctx_.extmax = doc.extmax;
    }

    void build_model_space() {
        scene_.nodes.push_back(named_node("Model space"));
        Sink sink;
        for (const Entity& e : doc_.entities) {
            if (e.integer(67) == 1) { ++paper_; continue; }
            process(e, Inherit{}, sink, 0, false);
        }
        emit(sink, 0, "Model space");
        for (const auto& [key, b] : doc_.blocks)
            if (key.rfind("*PAPER_SPACE", 0) == 0) paper_ += b.entities.size();
    }

    void report() {
        for (const auto& [why, n] : losses_) scene_.warnings.push_back(std::to_string(n) + " x " + why);
        if (paper_) scene_.notes.push_back(std::to_string(paper_) + " paper-space (layout) entities not shown; the model space is shown");
        if (hidden_) scene_.notes.push_back(std::to_string(hidden_) + " entities on layers that are off or frozen, not shown (as in AutoCAD)");
        for (const auto& [t, n] : non_graphic_) scene_.notes.push_back(std::to_string(n) + " " + t + " entities (not graphical)");
        for (const auto& [n, k] : notes_) scene_.notes.push_back(std::to_string(k) + " x " + n);
        if (attdefs_) scene_.notes.push_back(std::to_string(attdefs_) + " attribute definitions inside blocks (not displayed in block references)");
    }

private:
    void lose(const std::string& why) { ++losses_[why]; }

    const Layer* layer_of(std::string_view key) const {
        auto it = doc_.layers.find(key);
        return it == doc_.layers.end() ? nullptr : &it->second;
    }

    // Layer the entity displays on, after layer "0" substitution inside blocks.
    std::string display_layer(const Entity& e, const Inherit& in) const {
        std::string l = upper_key(e.str(8, "0"));
        if (l == "0" && !in.layer.empty()) l = in.layer;
        return l;
    }

    Color entity_color(const Entity& e, const Inherit& in, const std::string& layer) const {
        Color c;
        const long long aci = e.integer(62, 256);
        if (e.has(420)) {  // true colour overrides the colour number
            const long long rgb = e.integer(420);
            c = {std::uint8_t(rgb >> 16), std::uint8_t(rgb >> 8), std::uint8_t(rgb), 255};
        } else if (aci == 0) {
            c = in.byblock;
        } else if (aci == 256 || aci == 257) {
            const Layer* l = layer_of(layer);
            if (l && l->true_color >= 0) c = {std::uint8_t(l->true_color >> 16), std::uint8_t(l->true_color >> 8), std::uint8_t(l->true_color), 255};
            else c = aci_color(l ? l->aci : 7);
        } else {
            c = aci_color(int(aci < 0 ? -aci : aci));
        }
        const long long tr = e.integer(440, -1);
        if (tr >= 0 && (tr & 0x02000000)) c.a = std::uint8_t(tr & 0xFF);
        return c;
    }

    // Arrowhead of a LEADER at its first vertex, pointing away from the second: the DIMLDRBLK block
    // (leader override, then its dimension style, then the header) scaled by DIMASZ * DIMSCALE, or
    // AutoCAD's closed filled arrow (tip at the vertex, length 1, base half-width 1/6).
    void leader_arrow(const Entity& e, const Inherit& in, Color color, const std::string& layer, Sink& sink, std::uint32_t node) {
        if (e.integer(71, 1) == 0) return;
        std::vector<Vec3> v;
        for (const Group* g = e.begin; g != e.end; ++g) {
            if (g->code == 10) v.push_back({g->num, 0, 0});
            else if (g->code == 20 && !v.empty()) v.back().y = g->num;
            else if (g->code == 30 && !v.empty()) v.back().z = g->num;
        }
        if (v.size() < 2) return;
        DimStyle d = doc_.header_dim;
        if (!d.ldrblk.empty()) d.ldrblk = "NAME:" + upper_key(d.ldrblk);
        if (auto it = doc_.dim_styles.find(upper_key(e.str(3, "STANDARD"))); it != doc_.dim_styles.end()) {
            if (it->second.asz >= 0) d.asz = it->second.asz;
            if (it->second.scale >= 0) d.scale = it->second.scale;
            if (!it->second.ldrblk.empty()) d.ldrblk = it->second.ldrblk;
        }
        // overrides: XDATA ACAD { DSTYLE 1002 { 1070 code, value ... } }
        for (const Group* g = e.begin; g != e.end; ++g) {
            if (!(g->code == 1000 && g->str == "DSTYLE")) continue;
            for (const Group* q = g + 1; q + 1 < e.end && q->code != 1001; ++q) {
                if (q->code != 1070) continue;
                const long long code = q->ival;
                const Group* val = q + 1;
                if (code == 41 && val->code == 1040) d.asz = val->num;
                else if (code == 40 && val->code == 1040) d.scale = val->num;
                else if (code == 341 && val->code == 1005) d.ldrblk = upper_key(val->str);
                ++q;
            }
        }
        const double size = (d.asz >= 0 ? d.asz : 0.18) * (d.scale > 0 ? d.scale : 1.0);
        if (size <= 0) return;
        const Vec3 tip = v[0], dir = normalized(v[0] - v[1]);
        if (!d.ldrblk.empty()) {
            std::string name;
            if (d.ldrblk.rfind("NAME:", 0) == 0) name = d.ldrblk.substr(5);
            else if (auto it = doc_.block_record.find(d.ldrblk); it != doc_.block_record.end()) name = it->second;
            const Block* b = name.empty() ? nullptr : find_block(name);
            if (!b) { lose("LEADER: arrowhead block not found"); return; }
            const double c = dir.x * size, s = dir.y * size;
            Matrix4 m = kIdentity;
            m[0] = c; m[1] = s; m[4] = -s; m[5] = c; m[10] = size;
            m[12] = tip.x - (c * b->base.x - s * b->base.y);
            m[13] = tip.y - (s * b->base.x + c * b->base.y);
            m[14] = tip.z - size * b->base.z;
            add_child(node, b->name, m, block_node(*b, Inherit{color, layer}));
            (void)in;
            return;
        }
        const Vec3 side{-dir.y * size / 6, dir.x * size / 6, 0};
        const Vec3 base = tip - dir * size;
        sink.tri(tip, base + side, base - side, color);
    }

    // An active XCLIP: extension dictionary (360) -> ACAD_FILTER -> SPATIAL -> SPATIAL_FILTER with
    // clipping on (71 = 1).
    const Entity* clip_filter(const Entity& e) const {
        auto object = [&](std::string_view handle) -> const Entity* {
            auto it = doc_.objects.find(upper_key(handle));
            return it == doc_.objects.end() ? nullptr : &it->second;
        };
        auto entry = [&](const Entity* dict, std::string_view name) -> const Entity* {
            if (!dict || dict->type != "DICTIONARY") return nullptr;
            for (const Group* g = dict->begin; g != dict->end; ++g)
                if (g->code == 3 && g->str == name)
                    for (const Group* h = g + 1; h != dict->end && h->code != 3; ++h)
                        if (h->code == 350 || h->code == 360) return object(h->str);
            return nullptr;
        };
        const Group* x = e.find(360);
        if (!x) return nullptr;
        const Entity* filter = entry(entry(object(x->str), "ACAD_FILTER"), "SPATIAL");
        return filter && filter->type == "SPATIAL_FILTER" && filter->integer(71) == 1 ? filter : nullptr;
    }

    // The block's drawn content in block coordinates, nested references flattened.
    void flatten(std::uint32_t node, const Matrix4& m, Sink& out) const {
        const Node& n = scene_.nodes[node];
        auto at = [&](const float* p) { const auto w = cadcore::apply(m, p); return Vec3{w[0], w[1], w[2]}; };
        for (auto mi : n.meshes) {
            const Mesh& me = scene_.meshes[mi];
            if (me.topology == Topology::Points) {
                for (std::size_t i = 0; i < me.vertex_count(); ++i)
                    out.point(at(&me.positions[3 * i]), me.colors.size() >= 4 * (i + 1) ? Color{me.colors[4 * i], me.colors[4 * i + 1], me.colors[4 * i + 2], me.colors[4 * i + 3]} : me.color);
                continue;
            }
            for (std::size_t t = 0; t < me.triangle_count(); ++t) {
                Color c = me.color;
                for (const auto& g : me.groups)
                    if (t >= g.first && t < g.first + g.count) c = g.color;
                out.tri(at(&me.positions[3 * me.indices[3 * t]]), at(&me.positions[3 * me.indices[3 * t + 1]]), at(&me.positions[3 * me.indices[3 * t + 2]]), c);
            }
        }
        for (auto pi : n.polylines) {
            const Polyline& pl = scene_.polylines[pi];
            Polyline3 l;
            for (std::size_t k = 0; k + 2 < pl.points.size(); k += 3) l.push_back(at(&pl.points[k]));
            out.line(std::move(l), pl.color);
        }
        for (auto c : n.children) flatten(c, cadcore::multiply(m, scene_.nodes[c].transform), out);
    }

    // XCLIP: the block's content inside the filter's boundary, as a node in block coordinates.
    // Boundary points (10/20; two = a rectangle) become block coordinates through the filter's
    // first matrix (12 x 40, row-major 3 x 4: the inverse of the reference transform when the clip
    // was made, so the clip stays with the block when the reference moves). Returns the block
    // node itself (unclipped, with a reported loss) for filter forms not decoded.
    std::uint32_t clipped_node(const Block& b, std::uint32_t block, const Entity& f) {
        const std::string key = std::string(f.str(5)) + '\x1f' + std::to_string(block);
        if (auto it = clip_cache_.find(key); it != clip_cache_.end()) return it->second;
        std::vector<Vec3> pts;
        std::vector<double> m;
        for (const Group* g = f.begin; g != f.end; ++g) {
            if (g->code == 10) pts.push_back({g->num, 0, 0});
            else if (g->code == 20 && !pts.empty()) pts.back().y = g->num;
            else if (g->code == 40) m.push_back(g->num);
        }
        const bool plain_plane = f.num(210) == 0 && f.num(220) == 0 && f.num(230, 1) == 1 && f.num(11) == 0 && f.num(21) == 0 && f.num(31) == 0;
        const bool identity2 = m.size() >= 24 && m[12] == 1 && m[13] == 0 && m[14] == 0 && m[15] == 0 && m[16] == 0 && m[17] == 1 && m[18] == 0 &&
                               m[19] == 0 && m[20] == 0 && m[21] == 0 && m[22] == 1 && m[23] == 0;
        if (pts.size() < 2 || m.size() < 24 || !plain_plane || !identity2 || f.integer(72) == 1 || f.integer(73) == 1) {
            lose("INSERT clipping boundary (XCLIP) form not decoded (front/back clip, or a transformed boundary): the whole block is shown");
            return clip_cache_[key] = block;
        }
        if (pts.size() == 2) pts = {pts[0], {pts[1].x, pts[0].y, 0}, pts[1], {pts[0].x, pts[1].y, 0}};
        Polyline3 boundary;
        for (const auto& p : pts) boundary.push_back({m[0] * p.x + m[1] * p.y + m[3], m[4] * p.x + m[5] * p.y + m[7], 0});
        Sink content;
        flatten(block, kIdentity, content);
        Sink kept = clip_to_boundary(content, boundary);
        const std::uint32_t node = std::uint32_t(scene_.nodes.size());
        scene_.nodes.push_back(named_node(b.name + " (clipped)"));
        emit(kept, node, b.name + " (clipped)");
        return clip_cache_[key] = node;
    }

    Color layer_color(const std::string& layer) const {
        const Layer* l = layer_of(layer);
        if (l && l->true_color >= 0) return {std::uint8_t(l->true_color >> 16), std::uint8_t(l->true_color >> 8), std::uint8_t(l->true_color), 255};
        return aci_color(l ? l->aci : 7);
    }

    Matrix4 insert_matrix(const Entity& e, const Block& b, double col_off, double row_off) const {
        const Ocs o(e.point(210, {0, 0, 1}));
        const Vec3 p = e.point(10);
        const double sx = e.num(41, 1), sy = e.num(42, 1), sz = e.num(43, 1);
        const double r = e.num(50) * kDeg, cr = std::cos(r), sr = std::sin(r);
        // local = Rz * (cell + S * (x - base)); world = OCS * (p + local)
        auto lin = [&](const Vec3& v) {  // Rz * S * v
            const double x = v.x * sx, y = v.y * sy;
            return Vec3{cr * x - sr * y, sr * x + cr * y, v.z * sz};
        };
        const Vec3 cell{cr * col_off - sr * row_off, sr * col_off + cr * row_off, 0};
        const Vec3 t = p + cell - lin(b.base);
        Matrix4 m{};
        const Vec3 cols[3] = {lin({1, 0, 0}), lin({0, 1, 0}), lin({0, 0, 1})};
        for (int c = 0; c < 3; ++c) {
            const Vec3 w = o.to_wcs(cols[c]);
            m[std::size_t(c * 4)] = w.x; m[std::size_t(c * 4 + 1)] = w.y; m[std::size_t(c * 4 + 2)] = w.z;
        }
        const Vec3 wt = o.to_wcs(t);
        m[12] = wt.x; m[13] = wt.y; m[14] = wt.z; m[15] = 1;
        return m;
    }

    const Block* find_block(std::string_view name) const {
        auto it = doc_.blocks.find(upper_key(name));
        return it == doc_.blocks.end() ? nullptr : &it->second;
    }

    // The node holding block `b` decoded for context `in` (built on first use).
    std::uint32_t block_node(const Block& b, const Inherit& in) {
        const CacheKey key{&b, std::uint32_t(in.byblock.r) << 24 | std::uint32_t(in.byblock.g) << 16 | std::uint32_t(in.byblock.b) << 8 | in.byblock.a, in.layer};
        if (auto it = cache_.find(key); it != cache_.end()) return it->second;
        if (std::find(stack_.begin(), stack_.end(), upper_key(b.name)) != stack_.end())
            fail(ErrorCode::Malformed, "DXF block '" + b.name + "' inserts itself", b.offset);
        if (stack_.size() > 128) fail(ErrorCode::Malformed, "DXF block nesting deeper than 128", b.offset);
        stack_.push_back(upper_key(b.name));
        const std::uint32_t node = std::uint32_t(scene_.nodes.size());
        scene_.nodes.push_back(named_node(b.name));
        Sink sink;
        for (const Entity& e : b.entities) process(e, in, sink, node, true);
        emit(sink, node, b.name);
        stack_.pop_back();
        cache_[key] = node;
        return node;
    }

    void add_child(std::uint32_t parent, std::string name, const Matrix4& m, std::uint32_t child) {
        Node n = named_node(std::move(name));
        n.transform = m;
        n.children.push_back(child);
        const std::uint32_t idx = std::uint32_t(scene_.nodes.size());
        scene_.nodes.push_back(std::move(n));
        scene_.nodes[parent].children.push_back(idx);
    }

    void process(const Entity& e, const Inherit& in, Sink& sink, std::uint32_t node, bool in_block) {
        if (e.integer(60) == 1) return;  // invisible
        const std::string layer = display_layer(e, in);
        const Kind kind = classify(e.type);
        // a frozen layer hides everything on it; an off layer hides the entity itself, but a block
        // reference on it still shows the block's entities that are on other layers
        if (const Layer* l = layer_of(layer); l && (l->frozen || (l->off && kind != Kind::Insert && kind != Kind::BlockRef))) {
            ++hidden_;
            return;
        }
        const Color color = entity_color(e, in, layer);
        switch (kind) {
            case Kind::Draw: {
                std::vector<std::string> losses, notes;
                DrawContext ctx = dctx_;
                ctx.color = color;
                ctx.notes = &notes;
                ctx.mline_styles = &doc_.mline_styles;
                ctx.layer_color = layer_color(layer);
                draw_entity(e, ctx, sink, losses);
                if (e.type == "LEADER") leader_arrow(e, in, color, layer, sink, node);
                for (const auto& w : losses) lose(w);
                for (const auto& n : notes) ++notes_[n];
                break;
            }
            case Kind::Insert: {
                const Block* b = find_block(e.str(2));
                if (!b) { lose(std::string(e.type) + " of an undefined block"); break; }
                if (b->flags & 4) { lose("external reference '" + std::string(b->xref_path) + "' not loaded"); break; }
                Inherit child{color, layer};
                std::uint32_t bn = block_node(*b, child);
                if (const Entity* f = clip_filter(e)) bn = clipped_node(*b, bn, *f);
                const long long cols = std::max<long long>(1, e.integer(70, 1)), rows = std::max<long long>(1, e.integer(71, 1));
                if (cols * rows > 1'000'000) { lose("MINSERT array larger than 1,000,000 cells"); break; }
                for (long long r = 0; r < rows; ++r)
                    for (long long c = 0; c < cols; ++c)
                        add_child(node, b->name, insert_matrix(e, *b, double(c) * e.num(44), double(r) * e.num(45)), bn);
                for (const Entity& a : e.followers)
                    if (!(a.integer(70) & 1) && a.integer(60) != 1) lose("ATTRIB: text is not rendered (no font engine)");
                break;
            }
            case Kind::BlockRef: {
                const Block* b = e.has(2) ? find_block(e.str(2)) : nullptr;
                if (!b) { lose(std::string(e.type) + " without its geometry block (would need regenerating)"); break; }
                Inherit child{color, layer};
                const std::uint32_t bn = block_node(*b, child);
                // the block holds the dimension's graphics in its own coordinates (OCS of 210)
                const Ocs o(e.point(210, {0, 0, 1}));
                Matrix4 m = kIdentity;
                if (!o.identity()) {
                    m = {o.ax.x, o.ax.y, o.ax.z, 0, o.ay.x, o.ay.y, o.ay.z, 0, o.az.x, o.az.y, o.az.z, 0, 0, 0, 0, 1};
                }
                add_child(node, b->name, m, bn);
                break;
            }
            case Kind::Repeat: {
                // entities between REPEAT and ENDREP, placed on ENDREP's grid: 70 columns x 71 rows,
                // spacing 40 (x) and 41 (y)
                const std::uint32_t inner = std::uint32_t(scene_.nodes.size());
                scene_.nodes.push_back(named_node("REPEAT"));
                Sink rs;
                for (const Entity& c : e.followers) process(c, in, rs, inner, in_block);
                emit(rs, inner, "REPEAT");
                const long long cols = std::max<long long>(1, e.integer(70, 1)), rows = std::max<long long>(1, e.integer(71, 1));
                if (cols * rows > 1'000'000) { lose("REPEAT array larger than 1,000,000 cells"); break; }
                for (long long r = 0; r < rows; ++r)
                    for (long long c = 0; c < cols; ++c) {
                        Matrix4 m = kIdentity;
                        m[12] = double(c) * e.num(40);
                        m[13] = double(r) * e.num(41);
                        add_child(node, "REPEAT cell", m, inner);
                    }
                break;
            }
            case Kind::Text:
                if (e.type == "ATTDEF" && in_block) { ++attdefs_; break; }
                lose(std::string(e.type) + ": text is not rendered (no font engine)");
                break;
            case Kind::Brep:
                if (proxy(e, color, sink)) lose(std::string(e.type) + ": ACIS solid/surface shown as its stored wireframe (surfaces need a B-rep kernel)");
                else lose(std::string(e.type) + ": ACIS solid/surface not decoded (needs a B-rep kernel)");
                break;
            case Kind::NonGraphic:
                ++non_graphic_[std::string(e.type)];
                break;
            case Kind::Lost: {
                // AutoCAD's own graphics of the entity, when stored; raster and underlay content
                // is not part of them
                const bool drawn = proxy(e, color, sink);
                static const std::set<std::string_view> external{"IMAGE", "PDFUNDERLAY", "DWFUNDERLAY", "DGNUNDERLAY", "OLE2FRAME", "OLEFRAME", "POINTCLOUD", "POINTCLOUDEX"};
                if (e.type == "WIPEOUT" && drawn) lose("WIPEOUT: frame shown; masking of what lies behind it is not applied");
                else if (external.count(e.type)) lose(std::string(e.type) + (drawn ? ": frame shown; its image/underlay content is not decoded" : " entities not decoded"));
                else if (!drawn) lose(std::string(e.type) + " entities not decoded");
                break;
            }
        }
    }

    // Draws the entity's stored proxy graphics; false when it has none or they do not frame.
    bool proxy(const Entity& e, Color color, Sink& sink) {
        const std::vector<std::uint8_t> bytes = proxy_bytes(e);
        if (bytes.empty()) return false;
        std::vector<std::string> losses;
        if (!draw_proxy(bytes, color, sink, losses)) return false;
        for (const auto& w : losses) lose(std::string(e.type) + " " + w);
        return true;
    }

    // Moves the collected primitives into meshes / polylines of `node`. Geometry far from the
    // origin is stored relative to its centre (a child node carries the offset) so that float
    // positions keep their precision.
    void emit(Sink& s, std::uint32_t node, const std::string& name) {
        if (s.empty()) return;
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        auto grow = [&](const Vec3& p) {
            lo[0] = std::min(lo[0], p.x); lo[1] = std::min(lo[1], p.y); lo[2] = std::min(lo[2], p.z);
            hi[0] = std::max(hi[0], p.x); hi[1] = std::max(hi[1], p.y); hi[2] = std::max(hi[2], p.z);
        };
        for (const auto& p : s.tri_pts) grow(p);
        for (const auto& p : s.points) grow(p);
        for (const auto& l : s.lines)
            for (const auto& p : l.pts) grow(p);
        Vec3 origin{};
        const double far = std::max({std::fabs(lo[0]), std::fabs(hi[0]), std::fabs(lo[1]), std::fabs(hi[1]), std::fabs(lo[2]), std::fabs(hi[2])});
        if (far > 1e4) origin = {(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2};
        std::uint32_t target = node;
        if (!(origin == Vec3{})) {
            Matrix4 m = kIdentity;
            m[12] = origin.x; m[13] = origin.y; m[14] = origin.z;
            Node n = named_node(name + " (geometry)");
            n.transform = m;
            target = std::uint32_t(scene_.nodes.size());
            scene_.nodes.push_back(std::move(n));
            scene_.nodes[node].children.push_back(target);
        }
        auto put = [&](std::vector<float>& out, const Vec3& p) {
            out.insert(out.end(), {float(p.x - origin.x), float(p.y - origin.y), float(p.z - origin.z)});
        };
        if (!s.tri_pts.empty()) {
            Mesh m;
            m.name = name;
            m.kind = BodyKind::Surface;
            m.positions.reserve(s.tri_pts.size() * 3);
            for (const auto& p : s.tri_pts) put(m.positions, p);
            m.indices.resize(s.tri_pts.size());
            for (std::uint32_t i = 0; i < m.indices.size(); ++i) m.indices[i] = i;
            for (std::uint32_t t = 0; t < s.tri_colors.size(); ++t) {
                if (!m.groups.empty() && m.groups.back().color == s.tri_colors[t]) ++m.groups.back().count;
                else m.groups.push_back({t, 1, s.tri_colors[t]});
            }
            if (m.groups.size() == 1) { m.color = m.groups[0].color; m.groups.clear(); }
            scene_.nodes[target].meshes.push_back(std::uint32_t(scene_.meshes.size()));
            scene_.meshes.push_back(std::move(m));
        }
        if (!s.points.empty()) {
            Mesh m;
            m.name = name + " (points)";
            m.topology = Topology::Points;
            for (std::size_t i = 0; i < s.points.size(); ++i) {
                put(m.positions, s.points[i]);
                const Color c = s.point_colors[i];
                m.colors.insert(m.colors.end(), {c.r, c.g, c.b, c.a});
            }
            scene_.nodes[target].meshes.push_back(std::uint32_t(scene_.meshes.size()));
            scene_.meshes.push_back(std::move(m));
        }
        for (auto& l : s.lines) {
            Polyline pl;
            pl.color = l.color;
            pl.points.reserve(l.pts.size() * 3);
            for (const auto& p : l.pts) put(pl.points, p);
            scene_.nodes[target].polylines.push_back(std::uint32_t(scene_.polylines.size()));
            scene_.polylines.push_back(std::move(pl));
        }
        s = Sink{};
    }

    const Document& doc_;
    Scene& scene_;
    DrawContext dctx_;
    // a block decoded once per (block, BYBLOCK colour, layer replacing "0")
    struct CacheKey {
        const Block* block;
        std::uint32_t rgba;
        std::string layer;
        bool operator<(const CacheKey& o) const { return std::tie(block, rgba, layer) < std::tie(o.block, o.rgba, o.layer); }
    };
    std::map<CacheKey, std::uint32_t> cache_;
    std::map<std::string, std::uint32_t> clip_cache_;  // (filter handle, block node) -> clipped node
    std::vector<std::string> stack_;
    std::map<std::string, std::size_t> losses_;
    std::map<std::string, std::size_t> non_graphic_;
    std::map<std::string, std::size_t> notes_;
    std::size_t paper_ = 0, hidden_ = 0, attdefs_ = 0;
};

Scene decode(ByteView data, const DecodeOptions& options) {
    const GroupStream gs = read_groups(data, options.threads);
    const Document doc = parse_document(gs);
    Scene scene;
    scene.format = gs.binary ? "dxf-binary" : "dxf-ascii";
    scene.unit_to_mm = unit_to_mm(doc.insunits);
    Builder b(doc, scene);
    b.build_model_space();
    b.report();
    for (const auto& s : doc.unknown_sections) scene.notes.push_back("unknown section " + s + " skipped");
    if (!gs.saw_eof) scene.notes.push_back("no EOF marker");
    return scene;
}

}  // namespace

const Decoder& dxf() {
    static const SimpleDecoder d({"dxf", "AutoCAD Drawing Exchange Format", {"dxf"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
