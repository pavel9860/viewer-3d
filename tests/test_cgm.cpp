// CATIA CGMGeom: varint integers, framing and B-rep decode of corpus parts (skipped when the
// samples folder is absent). Ground truth for the meshes is CATIA's stored tessellation
// (research/cgm/compare_tess.py); these cases pin the decoded structure.
#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

#include "formats/catia/cgm.hpp"
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
