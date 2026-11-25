// A DXF file as sections: header variables, layer table, block definitions, entities.
#pragma once

#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "formats/dxf/geometry.hpp"
#include "formats/dxf/groups.hpp"

namespace cadcore::formats::autocad {

// One entity: its groups (after the "0 TYPE" group) and, for POLYLINE and INSERT, the
// VERTEX / ATTRIB entities that follow it up to SEQEND.
struct Entity {
    std::string_view type;
    const Group* begin = nullptr;
    const Group* end = nullptr;
    std::vector<Entity> followers;

    const Group* find(int code) const;
    double num(int code, double fallback = 0) const;
    long long integer(int code, long long fallback = 0) const;
    std::string_view str(int code, std::string_view fallback = {}) const;
    bool has(int code) const { return find(code) != nullptr; }
    // Point from codes c, c+10, c+20 (x, y, z); missing coordinates take `fallback`'s.
    Vec3 point(int code, Vec3 fallback = {}) const;
    std::size_t offset() const { return begin ? begin[-1].offset : 0; }
};

struct Layer {
    int aci = 7;              // colour number (negative in the file = layer off; stored positive here)
    long long true_color = -1;  // 0xRRGGBB or -1
    bool off = false, frozen = false;
};

struct Block {
    std::string name;
    Vec3 base;
    int flags = 0;                 // 1 anonymous, 4 external reference, 8 overlay, 16 dependent, 32 resolved
    std::string_view xref_path;
    std::vector<Entity> entities;
    std::size_t offset = 0;
};

// OBJECTS section MLINESTYLE: element offsets and colours, flags (caps, fill).
struct MlineStyle {
    int flags = 0;            // 1 fill, 16 / 256 square start / end cap, 32 / 512 inner arcs, 64 / 1024 round caps
    int fill_aci = 256;
    struct Element { double offset = 0; int aci = 256; };
    std::vector<Element> elements;
};

// Dimension style values that leaders use: arrow size (DIMASZ), overall scale (DIMSCALE) and the
// arrow block (DIMLDRBLK, a block record handle; empty = closed filled).
struct DimStyle {
    double asz = -1, scale = -1;
    std::string ldrblk;
};

struct Document {
    std::string version;                    // $ACADVER, e.g. "AC1015"; empty when absent
    int insunits = -1;                      // $INSUNITS, -1 when absent
    bool has_extents = false;
    bool extents_z = false;                 // the extents carry z (R11 and later)
    Vec3 extmin, extmax;                    // $EXTMIN / $EXTMAX
    std::map<std::string, Layer, std::less<>> layers;
    std::map<std::string, Block, std::less<>> blocks;
    std::vector<Entity> entities;           // ENTITIES section
    std::vector<std::string> unknown_sections;
    std::map<std::string, MlineStyle, std::less<>> mline_styles;  // by upper-case name
    std::map<std::string, DimStyle, std::less<>> dim_styles;      // by upper-case name
    std::map<std::string, std::string, std::less<>> block_record; // upper-case handle -> block name
    DimStyle header_dim;
    std::map<std::string, Entity, std::less<>> objects;            // OBJECTS section by upper-case handle (5)                                           // $DIMASZ, $DIMSCALE, $DIMLDRBLK (a block name)
};

// Table and block names are case-insensitive: maps are keyed by the ASCII upper-case name.
std::string upper_key(std::string_view name);

// Builds the document from a group stream (which must outlive it). Throws DecodeError on
// structural errors (an unterminated section or block).
Document parse_document(const GroupStream& gs);

}  // namespace cadcore::formats::autocad
