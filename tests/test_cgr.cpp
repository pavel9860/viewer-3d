// CGR decoder. No CATIA-exported counterpart of these CGRs exists in the samples, so the
// ground truth here is the grammar itself: the decoder throws unless every geometry and
// polyline node parses with exactly the bytes its framing declares (catia/tessellation.hpp).
// The geometry still needs the user's visual confirmation (renders/).
#include <filesystem>
#include <random>
#include <string>

#include "assembly_gt.hpp"
#include "cadcore/stats.hpp"
#include "formats/catia/store.hpp"
#include "formats/catia/tessellation.hpp"
#include "harness.hpp"

using namespace cadcore;
namespace fs = std::filesystem;

namespace {

const fs::path kTowplane = fs::path(CADCORE_SAMPLES_DIR) / "local/c00d1110";

}  // namespace

TEST_CASE("every project G CGR decodes node for node, without warnings") {
    std::size_t files = 0, triangles = 0;
    for (const auto& e : fs::directory_iterator(kTowplane)) {
        if (e.path().extension() != ".cgr") continue;
        const auto bytes = gt::read_file(e.path());
        const std::string name = e.path().filename().string();
        CAPTURE(name);
        Scene s;
        REQUIRE_NOTHROW(s = decode(ByteView(bytes.data(), bytes.size()), name));
        CHECK(s.format == "cgr");
        CHECK(s.warnings.empty());
        CHECK(!s.meshes.empty());
        triangles += measure(s).triangles;
        ++files;
    }
    CHECK(files == 24);
    std::printf("  %zu CGRs, %zu triangles\n", files, triangles);
}

TEST_CASE("CFV2 tessellation decoder: truncated and corrupted streams raise errors, never crash") {
    const auto file = gt::read_file(kTowplane / "NonAscii_1.cgr");
    const ByteView v(file.data(), file.size());
    const auto store = formats::catia::read_store(v);
    const formats::catia::StreamInfo* graph = nullptr;
    for (std::size_t i = 1; i < store.streams.size(); ++i)
        if (store.streams[i].name == "MainDataStream" && store.streams[i - 1].name == "SceneGraph") graph = &store.streams[i];
    REQUIRE(graph != nullptr);
    const auto stream = formats::catia::read_stream(v, store, *graph);
    const auto g = formats::catia::Generation::Cfv2;
    REQUIRE(formats::catia::decode_tessellation(ByteView(stream.data(), stream.size()), g).size() == 12);
    auto attempt = [&](const std::vector<std::uint8_t>& b) {
        try {
            (void)formats::catia::decode_tessellation(ByteView(b.data(), b.size()), g);
        } catch (const DecodeError&) {
        }
    };
    for (std::size_t n = 0; n < stream.size(); n += 11) attempt(std::vector<std::uint8_t>(stream.begin(), stream.begin() + std::ptrdiff_t(n)));
    std::mt19937 rng(11);
    for (int k = 0; k < 3000; ++k) {
        auto b = stream;
        b[rng() % b.size()] ^= std::uint8_t(1u << (rng() % 8));
        attempt(b);
    }
    // The whole file through the public API, truncated.
    for (std::size_t n = 0; n < file.size(); n += 97) {
        std::vector<std::uint8_t> b(file.begin(), file.begin() + std::ptrdiff_t(n));
        try {
            (void)decode(ByteView(b.data(), b.size()), "part.cgr");
        } catch (const DecodeError&) {
        }
    }
    CHECK(true);
}
