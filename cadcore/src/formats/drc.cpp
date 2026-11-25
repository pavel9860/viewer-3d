// Standalone Draco bitstream (.drc): one mesh decoded by codec::draco.
// Attribute types (Draco GeometryAttribute::Type): 0 POSITION, 1 NORMAL, 2 COLOR, 3 TEX_COORD, 4 GENERIC.
#include <cmath>
#include <cstring>
#include <string>

#include "cadcore/error.hpp"
#include "codec/draco.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

int sniff(ByteView d, std::string_view ext) {
    if (d.size >= 5 && std::memcmp(d.data, "DRACO", 5) == 0) return 2;
    return ext == "drc" ? 1 : 0;
}

std::uint8_t to_u8(float v, bool unit) {
    const float x = unit ? v * 255.f : v;
    return std::uint8_t(x <= 0 ? 0 : x >= 255 ? 255 : std::lround(x));
}

Scene decode(ByteView d, const DecodeOptions&) {
    const codec::draco::Mesh dm = codec::draco::decode(d);
    const codec::draco::Attribute* pos = nullptr;
    const codec::draco::Attribute* nrm = nullptr;
    const codec::draco::Attribute* col = nullptr;
    Scene s;
    s.format = "drc";
    for (const auto& a : dm.attributes) {
        const std::uint32_t id = a.unique_id;
        if (a.type == 0 && !pos) pos = &a;
        else if (a.type == 1 && !nrm && a.components == 3) nrm = &a;
        else if (a.type == 2 && !col && (a.components == 3 || a.components == 4)) col = &a;
        else if (a.type == 3 || a.type == 4) s.notes.push_back("Draco attribute " + std::to_string(id) + " (type " + std::to_string(a.type) + ") is not displayed");
        else s.notes.push_back("extra Draco attribute " + std::to_string(id) + " (type " + std::to_string(a.type) + ") ignored");
    }
    if (!pos || pos->components != 3) fail(ErrorCode::Malformed, "Draco: no 3-component POSITION attribute");
    Mesh m;
    m.positions = pos->values;
    m.indices = dm.faces;
    if (nrm) m.normals = nrm->values;
    if (col) {
        bool unit = true;
        for (float v : col->values) unit &= v <= 1.0f;
        m.colors.reserve(std::size_t(dm.num_points) * 4);
        for (std::uint32_t p = 0; p < dm.num_points; ++p) {
            const float* c = &col->values[std::size_t(p) * std::size_t(col->components)];
            m.colors.push_back(to_u8(c[0], unit));
            m.colors.push_back(to_u8(c[1], unit));
            m.colors.push_back(to_u8(c[2], unit));
            m.colors.push_back(col->components == 4 ? to_u8(c[3], unit) : 255);
        }
    }
    if (m.indices.empty()) m.topology = Topology::Points;
    s.meshes.push_back(std::move(m));
    return s;
}

}  // namespace

const Decoder& drc() {
    static const SimpleDecoder d({"drc", "Draco compressed mesh", {"drc"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
