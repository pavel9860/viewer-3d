// 3DXML with binary (V5_CFV3) representations, against ground truth:
//  * the assembly: occurrence table of the STEP export of the same product (assembly_gt.hpp);
//  * the tessellation: where the export kept CATIA's tessellation, each representation must
//    equal the mesh stored in its source CATPart triangle for triangle, bit for bit. Parts the
//    exporter re-tessellated (different sag) cannot be compared that way; their list is fixed
//    below so that any change in it is noticed.
// Plus hostile-input checks on the bit-stream decoder.
#include <filesystem>
#include <random>
#include <set>
#include <string>

#include "assembly_gt.hpp"
#include "cadcore/compare.hpp"
#include "formats/catia/tessellation.hpp"
#include "formats/catia/store.hpp"
#include "formats/zip.hpp"
#include "harness.hpp"

using namespace cadcore;
namespace fs = std::filesystem;

namespace {

const fs::path kAsmA = fs::path(CADCORE_SAMPLES_DIR) / "local/3fc767bf";

bool empty_product_rep(const std::string& w) { return w == "representation '1_4_0.3DRep' holds no tessellation"; }

Scene decode_file(const fs::path& p) {
    const auto bytes = gt::read_file(p);
    return decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), gt::folder_options(p.parent_path()));
}

// The meshes of one node as a scene of their own.
Scene node_scene(const Scene& s, const Node& n) {
    Scene out;
    Node root;
    for (auto m : n.meshes) {
        root.meshes.push_back(std::uint32_t(out.meshes.size()));
        out.meshes.push_back(s.meshes[m]);
    }
    out.nodes.push_back(root);
    return out;
}

}  // namespace

TEST_CASE("3DXML with binary representations matches the STEP export of the same product") {
    gt::check_assembly(kAsmA / "1.0.3dxml", "AsmA.occurrences.tsv", "3dxml", empty_product_rep);
}

TEST_CASE("3DXML binary representations equal the tessellation of their source CATParts") {
    const Scene s = decode_file(kAsmA / "1.0.3dxml");
    const std::set<std::string> retessellated{"1.1.1", "1.1.4", "1.1.5", "part-8568", "1.2.2", "part-9e0a"};
    std::set<std::string> differing;
    std::size_t compared = 0;
    for (const Node& n : s.nodes) {
        constexpr std::string_view suffix = "_ReferenceRep";
        if (n.name.size() <= suffix.size() || n.name.compare(n.name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;
        const std::string part = n.name.substr(0, n.name.size() - suffix.size());
        const fs::path cat = kAsmA / (part + ".CATPart");
        if (!fs::exists(cat)) continue;  // product-level representation
        CAPTURE(part);
        const auto ours = world_triangles(node_scene(s, n));
        const auto ref = world_triangles(decode_file(cat));
        const MatchReport r = match_triangles(ours, ref, 0.0);
        ++compared;
        if (!r.identical()) differing.insert(part);
        std::printf("  %-8s %7zu triangles, CATPart %7zu: %s\n", part.c_str(), r.a, r.b, r.identical() ? "identical" : "re-tessellated by the export");
    }
    CHECK(compared == 19);  // every CATPart of the product
    CHECK(differing == retessellated);
}

TEST_CASE("CFV3 skeleton decoder: truncated and corrupted streams raise errors, never crash") {
    const auto file = gt::read_file(kAsmA / "1.0.3dxml");
    formats::ZipArchive zip(ByteView(file.data(), file.size()));
    const auto rep = zip.read("1_1_3.3DRep");
    REQUIRE(rep.has_value());
    const ByteView v(rep->data(), rep->size());
    const auto store = formats::catia::read_store(v);
    const auto* sk = store.find("Skeleton");
    REQUIRE(sk != nullptr);
    const auto stream = formats::catia::read_stream(v, store, *sk);
    REQUIRE(formats::catia::decode_tessellation(ByteView(stream.data(), stream.size()), formats::catia::Generation::Cfv3).size() == 1);
    auto attempt = [](const std::vector<std::uint8_t>& b) {
        try {
            (void)formats::catia::decode_tessellation(ByteView(b.data(), b.size()), formats::catia::Generation::Cfv3);
        } catch (const DecodeError&) {
        }
    };
    for (std::size_t n = 0; n < stream.size(); ++n) attempt(std::vector<std::uint8_t>(stream.begin(), stream.begin() + std::ptrdiff_t(n)));
    std::mt19937 rng(7);
    for (int k = 0; k < 4000; ++k) {
        auto b = stream;
        b[rng() % b.size()] ^= std::uint8_t(1u << (rng() % 8));
        attempt(b);
    }
    // The whole store too (catalog and segment table).
    for (std::size_t n = 0; n < rep->size(); n += 7) {
        std::vector<std::uint8_t> b(rep->begin(), rep->begin() + std::ptrdiff_t(n));
        try {
            const auto st = formats::catia::read_store(ByteView(b.data(), b.size()));
            if (const auto* s2 = st.find("Skeleton")) (void)formats::catia::read_stream(ByteView(b.data(), b.size()), st, *s2);
        } catch (const DecodeError&) {
        }
    }
    CHECK(true);
}
