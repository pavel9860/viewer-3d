// CATIA object-stream grammar: relation headers of OSM streams across versions, one case per
// form seen in the corpus (formats/catia/osm.cpp, relation_header).
#include <cstdint>
#include <vector>

#include "cadcore/error.hpp"
#include "formats/catia/osm.hpp"
#include "harness.hpp"

using namespace cadcore;
using namespace cadcore::formats::catia;

namespace {

// Header tokens of a relation header `bytes` in a stream of `version`.
std::vector<HeaderToken> tokens(std::vector<std::uint8_t> bytes, int version) {
    static std::vector<std::vector<std::uint8_t>> keep;  // the ByteViews below point into these
    keep.push_back(std::move(bytes));
    Object o;
    o.osm_version = version;
    o.relation_header = ByteView(keep.back().data(), keep.back().size());
    return relation_header(o);
}

std::vector<std::uint8_t> link_blob() {  // e5 <len 18> 11 03 <uuid 16>
    std::vector<std::uint8_t> b{0xe5, 18, 0, 0, 0, 0x11, 0x03};
    for (int i = 0; i < 16; ++i) b.push_back(std::uint8_t(i));
    return b;
}

bool numbers(const std::vector<HeaderToken>& t, std::vector<std::uint64_t> want) {
    if (t.size() != want.size()) return false;
    for (std::size_t i = 0; i < t.size(); ++i)
        if (!t[i].is_number() || t[i].number != want[i]) return false;
    return true;
}

}  // namespace

TEST_CASE("OSM relation headers from version 35: kind byte, class names") {
    CHECK(numbers(tokens({0x14, 0x01, 0x83, 0x83}, 35), {2, 2}));  // bit 2: kind 01, then tokens
    CHECK(numbers(tokens({0x04, 0x01, 0x82}, 38), {1}));
    CHECK(numbers(tokens({0x01, 0x34}, 35), {}));                // form 1: a class name
    auto b = link_blob();
    std::vector<std::uint8_t> h{0x1c, 0x00};                     // kind 00: link blob, index
    h.insert(h.end(), b.begin(), b.end());
    h.push_back(0x83);
    const auto t = tokens(h, 36);
    REQUIRE(t.size() == 2);
    CHECK(t[0].is_blob);
    CHECK((t[1].is_number() && t[1].number == 2));
    std::vector<std::uint8_t> h5{0x05, 0x38};                    // form 5: class name, then tokens
    h5.insert(h5.end(), b.begin(), b.end());
    CHECK(tokens(h5, 35).size() == 1);
    bool threw = false;                                          // kind other than 00 / 01
    try { tokens({0x14, 0x07, 0x83}, 35); } catch (const DecodeError&) { threw = true; }
    CHECK(threw);
}

TEST_CASE("OSM relation headers before version 35: leading token, inline class, end marker") {
    CHECK(numbers(tokens({0x24, 0x83, 0xfe}, 30), {2}));          // tokens, closing fe
    CHECK(numbers(tokens({0x12, 0x83, 0x84}, 28), {3}));          // bit 1: one leading token
    auto b = link_blob();
    std::vector<std::uint8_t> h{0x16, 0x82};                     // form 6: leading token, blob, n, parent
    h.insert(h.end(), b.begin(), b.end());
    h.insert(h.end(), {0x85, 0x82});
    const auto t = tokens(h, 30);
    REQUIRE(t.size() == 3);
    CHECK(t[0].is_blob);
    CHECK((t[1].number == 4 && t[2].number == 1));
    CHECK(numbers(tokens({0x04, 0xfe, 0xfe, 0x3f, 0xfe, 0xfe}, 28), {}));  // inline class definition
    const auto f = tokens({0x04, 0xfd, 0x91, 0xbd}, 28);         // fd: startup index token
    REQUIRE(f.size() == 2);
    CHECK((f[0].is_index && f[0].number == 16));
}
