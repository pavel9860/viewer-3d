// Ground-truth tests: canonical parts with analytic expectations, written by the
// independent testgen writers, decoded by cadcore, compared exactly.
#include <cmath>
#include <set>

#include "cadcore/decode.hpp"
#include "cadcore/stats.hpp"
#include "harness.hpp"
#include "testgen/testgen.hpp"

using namespace cadcore;
using namespace cadcore::testgen;

namespace {

Scene decode_files(const Files& files) {
    DecodeOptions opt;
    opt.load_resource = [&](std::string_view path) -> std::optional<std::vector<std::uint8_t>> {
        for (const auto& [name, bytes] : files)
            if (name == path) return bytes;
        return std::nullopt;
    };
    const auto& main = files.front();
    return decode(ByteView(main.second.data(), main.second.size()), main.first, opt);
}

bool near(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

std::set<std::uint32_t> model_colors(const Model& m) {
    std::set<std::uint32_t> s;
    for (const auto& b : m.bodies)
        for (const auto& f : b.faces) s.insert(std::uint32_t(f.color.r) << 24 | f.color.g << 16 | f.color.b << 8 | f.color.a);
    return s;
}

std::set<std::uint32_t> scene_colors(const Scene& sc) {
    std::set<std::uint32_t> s;
    auto add = [&](Color c) { s.insert(std::uint32_t(c.r) << 24 | c.g << 16 | c.b << 8 | c.a); };
    for (const auto& m : sc.meshes) {
        if (!m.colors.empty()) {
            for (std::size_t i = 0; i + 3 < m.colors.size(); i += 4) add({m.colors[i], m.colors[i + 1], m.colors[i + 2], m.colors[i + 3]});
        } else if (!m.groups.empty()) {
            for (const auto& g : m.groups) add(g.color);
        } else {
            add(m.color);
        }
    }
    return s;
}

}  // namespace

TEST_CASE("generated parts decode to their analytic geometry in every format") {
    for (const Model& model : canonical_models()) {
        for (const Variant& v : variants()) {
            const bool catpart = v.id.rfind("catpart", 0) == 0;
            if (catpart && !model.instances.empty()) continue;  // a CATPart has no instances
            if (v.id == "catpart-r17-counts" && model.name == "dense_cylinder35k") continue;  // u16 counts cannot hold a 70002-ref strip
            CAPTURE(model.name);
            CAPTURE(v.id);
            const Files files = v.write(model);
            Scene sc;
            REQUIRE_NOTHROW(sc = decode_files(files));
            const SceneStats st = measure(sc);
            CHECK(st.triangles == model.expected.triangles);
            CHECK(near(st.area, model.expected.area, 1e-5));
            if (!std::isnan(model.expected.volume)) CHECK(near(st.signed_volume, model.expected.volume, 1e-5));
            for (int k = 0; k < 3; ++k) {
                CHECK(near(st.min[k], model.expected.min[k], 1e-5));
                CHECK(near(st.max[k], model.expected.max[k], 1e-5));
            }
            // Face colours survive in every format that stores them (STL stores none).
            if (v.id.rfind("stl", 0) != 0) CHECK(scene_colors(sc) == model_colors(model));
            if (catpart) {
                REQUIRE(sc.meshes.size() == model.bodies.size());
                for (std::size_t b = 0; b < model.bodies.size(); ++b)
                    CHECK(sc.meshes[b].kind == (model.bodies[b].solid ? BodyKind::Solid : BodyKind::Surface));
            }
            CHECK(sc.warnings.empty());
        }
    }
}

TEST_CASE("format detection picks the right decoder from content") {
    const Model cube = canonical_models().front();
    for (const Variant& v : variants()) {
        const Files f = v.write(cube);
        const auto& bytes = f.front().second;
        // Content-identifiable formats must be recognised even with a misleading file name.
        if (v.id == "obj" || v.id == "3mf") continue;
        const Decoder* d = detect(ByteView(bytes.data(), bytes.size()), "unnamed.bin");
        CAPTURE(v.id);
        REQUIRE(d != nullptr);
        const std::string family = v.id.substr(0, v.id.find('-'));
        const std::string expected = family == "glb" ? "gltf" : family;
        CHECK(std::string(d->info().id) == expected);
    }
}
