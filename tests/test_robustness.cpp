// Hostile-input tests: truncated and corrupted files must either decode or raise
// DecodeError, never crash, hang or throw anything else. Run under the sanitizer
// build (CADCORE_SANITIZE=ON) to also catch memory errors.
#include <cstring>
#include <random>

#include "cadcore/decode.hpp"
#include "harness.hpp"
#include "testgen/testgen.hpp"

using namespace cadcore;
using namespace cadcore::testgen;

namespace {

void decode_or_error(const std::vector<std::uint8_t>& bytes, const std::string& name) {
    try {
        Scene s = decode(ByteView(bytes.data(), bytes.size()), name);
        (void)s;
    } catch (const DecodeError&) {
        // expected for broken input
    }
}

}  // namespace

TEST_CASE("truncated files never crash") {
    const auto models = canonical_models();
    for (const Variant& v : variants()) {
        const Files f = v.write(models[1]);  // cylinder32: strips, fans, normals
        const auto& full = f.front().second;
        CAPTURE(v.id);
        for (int k = 0; k <= 97; ++k) {
            std::size_t cut = full.size() * std::size_t(k) / 97;
            std::vector<std::uint8_t> part(full.begin(), full.begin() + std::ptrdiff_t(cut));
            CHECK_NOTHROW(decode_or_error(part, f.front().first));
        }
    }
}

TEST_CASE("byte flips never crash") {
    std::mt19937 rng(12345);
    const auto models = canonical_models();
    for (const Variant& v : variants()) {
        const Files f = v.write(models[0]);
        CAPTURE(v.id);
        for (int trial = 0; trial < 300; ++trial) {
            auto bytes = f.front().second;
            const int flips = 1 + int(rng() % 8);
            for (int i = 0; i < flips; ++i) bytes[rng() % bytes.size()] ^= std::uint8_t(1u << (rng() % 8));
            CHECK_NOTHROW(decode_or_error(bytes, f.front().first));
        }
    }
}

TEST_CASE("empty and tiny inputs are rejected cleanly") {
    for (const char* name : {"a.stl", "a.obj", "a.ply", "a.off", "a.gltf", "a.glb", "a.3mf", "a.CATPart", "a.unknown"}) {
        std::vector<std::uint8_t> empty;
        CHECK_NOTHROW(decode_or_error(empty, name));
        std::vector<std::uint8_t> one{0x42};
        CHECK_NOTHROW(decode_or_error(one, name));
    }
}

TEST_CASE("decoders do not depend on the alignment of the caller's buffer") {
    const auto models = canonical_models();
    for (const Variant& v : variants()) {
        const Files f = v.write(models[0]);
        if (f.size() > 1) continue;  // multi-file variants need a resource loader
        const auto& bytes = f.front().second;
        std::vector<std::uint8_t> shifted(bytes.size() + 1);
        std::memcpy(shifted.data() + 1, bytes.data(), bytes.size());
        CAPTURE(v.id);
        Scene s;
        REQUIRE_NOTHROW(s = decode(ByteView(shifted.data() + 1, bytes.size()), f.front().first));
        CHECK(s.triangle_count() == 12);
    }
}
