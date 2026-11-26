// CATIA CGMGeom: varint integers, framing and B-rep decode of corpus parts (skipped when the
// samples folder is absent). Ground truth for the meshes is CATIA's stored tessellation
// (research/cgm/compare_tess.py); these cases pin the decoded structure.
#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

#include "formats/catia/cgm.hpp"
#include "formats/catia/cgm_compact.hpp"
#include "formats/catia/cgm_model.hpp"
#include "harness.hpp"

using namespace cadcore;
namespace cgm = cadcore::formats::catia::cgm;
namespace fs = std::filesystem;

namespace {
std::vector<std::uint8_t> load(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), {}};
}
std::uint32_t var(std::vector<std::uint8_t> b) {
    std::size_t k = 0;
    return cgm::varint(ByteView(b.data(), b.size()), k);
}
}  // namespace

TEST_CASE("CGM varint: small values and byte-presence tags") {
    CHECK(var({0x80}) == 0);
    CHECK(var({0x85}) == 5);
    CHECK(var({0x08, 0xda}) == 0xda);
    CHECK(var({0x10, 0x01}) == 0x100);           // byte 1 only
    CHECK(var({0x18, 0x01, 0x01}) == 0x101);
    CHECK(var({0x30, 0x84, 0x02}) == 0x028400);  // bytes 1 and 2
    CHECK(var({0x38, 0x93, 0x9f, 0x01}) == 0x019f93);
}

TEST_CASE("CGM corpus parts: framing and B-rep decode") {
    const fs::path s(CADCORE_SAMPLES_DIR);
    const fs::path part4 = s / "caa/caadoc/CAAElecDeviceItf.edu/CNext/resources/graphic/CAAElecConnect/Part4.CATPart";
    if (!fs::exists(part4)) return;
    const auto file = load(part4);
    const auto st = cgm::find_stream(ByteView(file.data(), file.size()));
    REQUIRE(st.has_value());
    const cgm::Stream fr = cgm::frame(ByteView(st->data(), st->size()));
    CHECK(fr.dialect == cgm::Dialect::Varint);
    CHECK(fr.records.size() == 221);
    const cgm::Ir ir = cgm::decode(fr);
    CHECK(ir.faces.size() == 8);
    CHECK(ir.skipped.empty());
    std::size_t solids = 0;
    for (const auto& sh : ir.shells) solids += sh.solid;
    CHECK(solids == 1);
    const auto bin = cgm::serialize(ir);
    CHECK(bin.size() > 8);
    CHECK(std::string(reinterpret_cast<const char*>(bin.data()), 5) == "CGMIR");
}

TEST_CASE("CGM compact records: canonical form") {
    // 2011Spring edge (compact integers, no edge-curve reference): v0 v1 p0 p1 as self - target
    const std::vector<std::uint8_t> edge = {0x09, 0x0f, 0x07, 0x0b, 0x21};
    cgm::Record r;
    r.id = 100;
    r.type = 511;
    r.body = ByteView(edge.data(), edge.size());
    const cgm::Canonical c = cgm::canonical(r, cgm::CompactLayout{true, false});
    REQUIRE(c.refs.size() == 4);
    CHECK(c.refs[0] == 98);   // 09: +2
    CHECK(c.refs[1] == 104);  // 0f: -4
    CHECK(c.refs[2] == 102);  // 07: -2
    CHECK(c.refs[3] == 103);  // 0b: -3
    CHECK((c.payload == std::vector<std::uint8_t>{0x21}));
    // shell: one face, bits kind 3 (11), sense +1; closing -1 (raw i16 when integers are not compact)
    const std::vector<std::uint8_t> shell = {0x81, 0x03, 0x13, 0xff, 0xff};
    r.type = 520;
    r.body = ByteView(shell.data(), shell.size());
    const cgm::Canonical s = cgm::canonical(r, cgm::CompactLayout{false, true});
    REQUIRE(s.refs.size() == 1);
    CHECK(s.refs[0] == 101);
    CHECK((s.payload == std::vector<std::uint8_t>{0x81, 3, 0, 1, 0, 0xff, 0xff}));
}

TEST_CASE("CGM compact corpus part (2011Spring): B-rep decode") {
    const fs::path p = fs::path(CADCORE_SAMPLES_DIR) / "caa/caadoc/CAACommonLayoutItf.edu/CNext/resources/graphic/CAACloPspEduIn/GASKET-RF-150-4in.CATPart";
    if (!fs::exists(p)) return;
    const auto file = load(p);
    const auto st = cgm::find_stream(ByteView(file.data(), file.size()));
    REQUIRE(st.has_value());
    const cgm::Stream fr = cgm::frame(ByteView(st->data(), st->size()));
    CHECK(fr.dialect == cgm::Dialect::Compact);
    REQUIRE(cgm::compact_layout(fr).has_value());
    const cgm::Ir ir = cgm::decode(fr);
    CHECK(ir.faces.size() == 8);
    CHECK(ir.skipped.empty());
}

TEST_CASE("CGM compact corpus part (2011Fall): B-rep decode") {
    const fs::path p = fs::path(CADCORE_SAMPLES_DIR) / "caa/caadoc/CAAAssemblyUI.edu/CNext/resources/graphic/Nut.CATPart";
    if (!fs::exists(p)) return;
    const auto file = load(p);
    const auto st = cgm::find_stream(ByteView(file.data(), file.size()));
    REQUIRE(st.has_value());
    const cgm::Stream fr = cgm::frame(ByteView(st->data(), st->size()));
    REQUIRE(fr.features.size() >= 3);
    CHECK(fr.features[0] == 0x06);
    const auto layout = cgm::compact_layout(fr);
    REQUIRE(layout.has_value());
    CHECK(layout->compact_ints);
    CHECK(layout->edge_curve_ref);
    CHECK(!layout->derived);
    const cgm::Ir ir = cgm::decode(fr);
    CHECK(ir.faces.size() == 40);
    CHECK(ir.skipped.empty());
}

TEST_CASE("CGM points: parameters are offsets within a Hermite segment") {
    // Points store (offset, segment number); segment numbers are shifted by the knot vector's offset
    // field (-7 on trimmed curves of this part). Every coedge range must lie on its Hermite curve.
    const fs::path p = fs::path(CADCORE_SAMPLES_DIR) / "nist/nist-ctc-pmi-cad-models/NIST-CTC-PMI/CATIA V5R21/nist_ctc_02_asme1_ct5210_rc.CATPart";
    if (!fs::exists(p)) return;
    const auto file = load(p);
    const auto st = cgm::find_stream(ByteView(file.data(), file.size()));
    REQUIRE(st.has_value());
    const cgm::Ir ir = cgm::decode(cgm::frame(ByteView(st->data(), st->size())));
    std::size_t hermite = 0, outside = 0;
    for (const auto& c : ir.coedges) {
        const auto& pc = ir.pcurves[c.pcurve];
        if (pc.kind != cgm::PcHermite) continue;
        const auto n = std::size_t(pc.p[0]);
        const double k0 = pc.p[1], k1 = pc.p[n], tol = 1e-9 * (1 + std::abs(k1));
        ++hermite;
        if (std::min(c.t0, c.t1) < k0 - tol || std::max(c.t0, c.t1) > k1 + tol) ++outside;
    }
    CHECK(hermite > 100);
    CHECK(outside == 0);
    for (const auto& k : ir.skipped) CHECK(k.why.find("point segment") == std::string::npos);
}
