#include "formats/dxf/document.hpp"

#include <cstdlib>

#include "cadcore/error.hpp"

namespace cadcore::formats::autocad {

const Group* Entity::find(int code) const {
    for (const Group* g = begin; g != end; ++g)
        if (g->code == code) return g;
    return nullptr;
}
double Entity::num(int code, double fallback) const { const Group* g = find(code); return g ? g->num : fallback; }
long long Entity::integer(int code, long long fallback) const { const Group* g = find(code); return g ? g->ival : fallback; }
std::string_view Entity::str(int code, std::string_view fallback) const { const Group* g = find(code); return g ? g->str : fallback; }
Vec3 Entity::point(int code, Vec3 f) const { return {num(code, f.x), num(code + 10, f.y), num(code + 20, f.z)}; }

std::string upper_key(std::string_view s) {
    std::string k(s);
    for (char& c : k)
        if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
    return k;
}

namespace {

bool is0(const Group& g, std::string_view v) { return g.code == 0 && g.str == v; }

// Entities in [p, end) up to (not including) a "0 stop" group; p is left on the stop group.
std::vector<Entity> read_entities(const Group*& p, const Group* end, std::string_view stop) {
    std::vector<Entity> out;
    while (p != end && !is0(*p, stop) && !is0(*p, "ENDSEC")) {
        if (p->code != 0) { ++p; continue; }  // stray group between entities
        if (is0(*p, "REPEAT")) {  // R2.x array: the entities up to ENDREP, repeated by ENDREP's grid
            const std::size_t at = p->offset;
            ++p;
            while (p != end && p->code != 0) ++p;
            Entity rep;
            rep.type = "REPEAT";
            rep.followers = read_entities(p, end, "ENDREP");
            if (p == end || !is0(*p, "ENDREP")) fail(ErrorCode::Malformed, "DXF REPEAT without ENDREP", at);
            rep.begin = ++p;
            while (p != end && p->code != 0) ++p;
            rep.end = p;
            out.push_back(std::move(rep));
            continue;
        }
        Entity e;
        e.type = p->str;
        e.begin = ++p;
        while (p != end && p->code != 0) ++p;
        e.end = p;
        const bool has_followers = e.type == "POLYLINE" || (e.type == "INSERT" && e.integer(66) == 1);
        if (has_followers) {
            const std::string_view child = e.type == "POLYLINE" ? "VERTEX" : "ATTRIB";
            while (p != end && is0(*p, child)) {
                Entity v;
                v.type = p->str;
                v.begin = ++p;
                while (p != end && p->code != 0) ++p;
                v.end = p;
                e.followers.push_back(v);
            }
            if (p != end && is0(*p, "SEQEND")) {
                ++p;
                while (p != end && p->code != 0) ++p;
            }
        }
        out.push_back(std::move(e));
    }
    return out;
}

void read_header(const Group* p, const Group* end, Document& doc) {
    for (; p != end; ++p) {
        if (p->code != 9) continue;
        const std::string_view var = p->str;
        auto value = [&](int code) -> const Group* {
            for (const Group* q = p + 1; q != end && q->code != 9; ++q)
                if (q->code == code) return q;
            return nullptr;
        };
        if (var == "$ACADVER") { if (auto g = value(1)) doc.version = std::string(g->str); }
        else if (var == "$INSUNITS") { if (auto g = value(70)) doc.insunits = int(g->ival); }
        else if (var == "$DIMASZ") { if (auto g = value(40)) doc.header_dim.asz = g->num; }
        else if (var == "$DIMSCALE") { if (auto g = value(40)) doc.header_dim.scale = g->num; }
        else if (var == "$DIMLDRBLK") { if (auto g = value(1)) doc.header_dim.ldrblk = std::string(g->str); }
        else if (var == "$EXTMIN" || var == "$EXTMAX") {
            auto x = value(10), y = value(20), z = value(30);
            if (x && y) {
                Vec3 v{x->num, y->num, z ? z->num : 0};
                (var == "$EXTMIN" ? doc.extmin : doc.extmax) = v;
                doc.has_extents = true;
                doc.extents_z = z != nullptr;
            }
        }
    }
    // AutoCAD writes +/-1e20 extents for an empty drawing
    if (doc.has_extents && (doc.extmin.x > doc.extmax.x || std::abs(doc.extmin.x) >= 1e19)) doc.has_extents = false;
}

void read_tables(const Group* p, const Group* end, Document& doc) {
    while (p != end) {
        if (is0(*p, "DIMSTYLE") || is0(*p, "BLOCK_RECORD")) {
            const bool dim = p->str == "DIMSTYLE";
            Entity e;
            e.begin = ++p;
            while (p != end && p->code != 0) ++p;
            e.end = p;
            if (!e.has(2)) continue;
            if (dim) {
                DimStyle d;
                d.asz = e.num(41, -1);
                d.scale = e.num(40, -1);
                d.ldrblk = upper_key(e.str(341));
                doc.dim_styles[upper_key(e.str(2))] = d;
            } else if (e.has(5)) {
                doc.block_record[upper_key(e.str(5))] = std::string(e.str(2));
            }
            continue;
        }
        if (!is0(*p, "LAYER")) { ++p; continue; }
        const Group* b = ++p;
        while (p != end && p->code != 0) ++p;
        Entity e;
        e.type = "LAYER";
        e.begin = b;
        e.end = p;
        if (!e.has(2)) continue;
        Layer l;
        const long long c = e.integer(62, 7);
        l.off = c < 0;
        l.aci = int(c < 0 ? -c : c);
        l.true_color = e.has(420) ? (e.integer(420) & 0xFFFFFF) : -1;
        l.frozen = (e.integer(70) & 1) != 0;
        doc.layers[upper_key(e.str(2))] = l;
    }
}

void read_objects(const Group* p, const Group* end, Document& doc) {
    for (const Group* q = p; q != end;) {  // index every object by handle
        if (q->code != 0) { ++q; continue; }
        Entity e;
        e.type = q->str;
        e.begin = ++q;
        while (q != end && q->code != 0) ++q;
        e.end = q;
        if (e.has(5)) doc.objects[upper_key(e.str(5))] = e;
    }
    while (p != end) {
        if (!is0(*p, "MLINESTYLE")) { ++p; continue; }
        ++p;
        MlineStyle st;
        std::string name;
        bool in_elements = false;
        for (; p != end && p->code != 0; ++p) {
            switch (p->code) {
                case 2: name = std::string(p->str); break;
                case 70: st.flags = int(p->ival); break;
                case 71: in_elements = true; break;
                case 49: st.elements.push_back({p->num, 256}); break;
                case 62:
                    if (in_elements && !st.elements.empty()) st.elements.back().aci = int(p->ival);
                    else if (!in_elements) st.fill_aci = int(p->ival);
                    break;
                default: break;
            }
        }
        if (!name.empty()) doc.mline_styles[upper_key(name)] = std::move(st);
    }
}

void read_blocks(const Group*& p, const Group* end, Document& doc) {
    while (p != end && !is0(*p, "ENDSEC")) {
        if (!is0(*p, "BLOCK")) { ++p; continue; }
        const std::size_t at = p->offset;
        Entity head;
        head.type = "BLOCK";
        head.begin = ++p;
        while (p != end && p->code != 0) ++p;
        head.end = p;
        Block b;
        b.name = std::string(head.str(2));
        b.base = head.point(10);
        b.flags = int(head.integer(70));
        b.xref_path = head.str(1);
        b.offset = at;
        b.entities = read_entities(p, end, "ENDBLK");
        if (p == end || !is0(*p, "ENDBLK")) fail(ErrorCode::Malformed, "DXF block '" + b.name + "' has no ENDBLK", at);
        ++p;
        while (p != end && p->code != 0) ++p;
        doc.blocks[upper_key(b.name)] = std::move(b);
    }
}

}  // namespace

Document parse_document(const GroupStream& gs) {
    Document doc;
    const Group* p = gs.groups.data();
    const Group* end = p + gs.groups.size();
    while (p != end) {
        if (is0(*p, "EOF")) break;
        if (!is0(*p, "SECTION")) { ++p; continue; }
        const std::size_t at = p->offset;
        ++p;
        if (p == end || p->code != 2) fail(ErrorCode::Malformed, "DXF SECTION without a name", at);
        const std::string_view name = p->str;
        ++p;
        const Group* body = p;
        if (name == "ENTITIES") {
            doc.entities = read_entities(p, end, "ENDSEC");
        } else if (name == "BLOCKS") {
            read_blocks(p, end, doc);
        } else {
            while (p != end && !is0(*p, "ENDSEC")) ++p;
            if (name == "HEADER") read_header(body, p, doc);
            else if (name == "TABLES") read_tables(body, p, doc);
            else if (name == "OBJECTS") read_objects(body, p, doc);
            else if (name != "CLASSES" && name != "THUMBNAILIMAGE" && name != "ACDSDATA")
                doc.unknown_sections.push_back(std::string(name));
        }
        if (p == end) fail(ErrorCode::Truncated, "DXF section " + std::string(name) + " has no ENDSEC", at);
        ++p;  // ENDSEC
    }
    return doc;
}

}  // namespace cadcore::formats::autocad
