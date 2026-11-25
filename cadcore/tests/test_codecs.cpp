// Codec tests: DEFLATE, CRC-32 and the zip reader, against known content.
#include <iterator>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <cstring>
#include <charconv>
#include "core/text.hpp"
#include <random>
#include <string>

#include "cadcore/error.hpp"
#include "codec/inflate.hpp"
#include "formats/zip.hpp"
#include "harness.hpp"
#include "testgen/zipwriter.hpp"

using namespace cadcore;

namespace {

std::string sample_text(std::size_t n, unsigned seed, bool repetitive) {
    std::mt19937 rng(seed);
    std::string s;
    s.reserve(n);
    const char* words[] = {"vertex ", "facet ", "normal ", "0.125 ", "-3.5e2 ", "endloop\n", "<triangle v1=\"", "\"/>"};
    while (s.size() < n) {
        if (repetitive) s += words[rng() % 8];
        else s += char(rng() & 0xFF);
    }
    s.resize(n);
    return s;
}

}  // namespace

TEST_CASE("CRC-32 matches the standard check value") {
    const std::string check = "123456789";
    CHECK(codec::crc32(reinterpret_cast<const std::uint8_t*>(check.data()), check.size()) == 0xCBF43926u);
    CHECK(codec::crc32(nullptr, 0) == 0u);
}

TEST_CASE("deflate round trip through our encoder and decoder") {
    for (std::size_t n : {0u, 1u, 2u, 3u, 100u, 4096u, 70000u, 300000u}) {
        for (bool rep : {true, false}) {
            CAPTURE(n);
            CAPTURE(rep);
            const std::string src = sample_text(n, unsigned(n) + rep, rep);
            const auto packed = testgen::deflate_fixed(reinterpret_cast<const std::uint8_t*>(src.data()), src.size());
            const auto out = codec::inflate(ByteView(packed.data(), packed.size()), src.size(), src.size());
            REQUIRE(out.size() == src.size());
            CHECK(std::equal(out.begin(), out.end(), reinterpret_cast<const std::uint8_t*>(src.data())));
        }
    }
}

TEST_CASE("stored DEFLATE blocks decode") {
    // final stored block: header bits (1, 00), aligned, LEN=5, NLEN=~5, "hello"
    const std::uint8_t raw[] = {0x01, 0x05, 0x00, 0xFA, 0xFF, 'h', 'e', 'l', 'l', 'o'};
    const auto out = codec::inflate(ByteView(raw, sizeof raw));
    CHECK(std::string(out.begin(), out.end()) == "hello");
}

TEST_CASE("zip reader returns entries byte-exact and verifies CRC-32") {
    const std::string a = sample_text(50000, 7, true), b = sample_text(3000, 9, false);
    auto zip = testgen::write_zip({{"dir/A.txt", a, true}, {"b.bin", b, false}});
    formats::ZipArchive arc(ByteView(zip.data(), zip.size()));
    auto ra = arc.read("DIR/a.TXT");  // part names compare case-insensitively
    REQUIRE(ra.has_value());
    CHECK(std::string(ra->begin(), ra->end()) == a);
    auto rb = arc.read("/b.bin");
    REQUIRE(rb.has_value());
    CHECK(std::string(rb->begin(), rb->end()) == b);
    CHECK(!arc.read("missing").has_value());

    // Flip one byte of the stored entry's data: the CRC check must catch it.
    const std::size_t pos = std::string(zip.begin(), zip.end()).find(b.substr(0, 16));
    REQUIRE(pos != std::string::npos);
    zip[pos + 5] ^= 0x01;
    formats::ZipArchive bad(ByteView(zip.data(), zip.size()));
    bool threw = false;
    try {
        bad.read("b.bin");
    } catch (const DecodeError& e) {
        threw = std::string(e.what()).find("CRC") != std::string::npos;
    }
    CHECK(threw);
}

TEST_CASE("truncated DEFLATE streams raise errors, never crash") {
    const std::string src = sample_text(20000, 3, true);
    const auto packed = testgen::deflate_fixed(reinterpret_cast<const std::uint8_t*>(src.data()), src.size());
    for (std::size_t cut = 0; cut < packed.size(); cut += 97) {
        bool ok = true;
        try {
            auto out = codec::inflate(ByteView(packed.data(), cut), src.size(), src.size());
            ok = out.size() <= src.size();
        } catch (const DecodeError&) {
        }
        CHECK(ok);
    }
}

TEST_CASE("text::parse_real gives exactly std::from_chars results") {
    std::uint64_t st = 0x9E3779B97F4A7C15ull;
    auto rnd = [&] { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; };
    long checked = 0, mismatches = 0;
    auto check = [&](const std::string& t) {
        float f1 = 0, f2 = 0;
        double d1 = 0, d2 = 0;
        const char* b = t.data();
        const char* e = b + t.size();
        const char* p1 = text::parse_real(b, e, f1);
        const auto r1 = std::from_chars(b, e, f2, std::chars_format::general);
        const char* p2 = text::parse_real(b, e, d1);
        const auto r2 = std::from_chars(b, e, d2, std::chars_format::general);
        ++checked;
        const bool okf = (p1 == nullptr) == (r1.ec != std::errc()) && (!p1 || (p1 == r1.ptr && std::memcmp(&f1, &f2, 4) == 0));
        const bool okd = (p2 == nullptr) == (r2.ec != std::errc()) && (!p2 || (p2 == r2.ptr && std::memcmp(&d1, &d2, 8) == 0));
        if (!okf || !okd) { ++mismatches; CAPTURE(t); CHECK(okf); CHECK(okd); }
    };
    for (const char* t : {"0", "-0", "1", "0.1", ".5", "5.", "1e5", "1E-5", "1e", "-.25e+3", "inf", "-nan", "x", "", "-", ".", "123456789012345678901234",
                          "0.000000000000000000000001", "3.4028235e38", "3.4028236e38", "1e-45", "1.17549435e-38", "16777217", "0.1e400"})
        check(t);
    for (int i = 0; i < 300000; ++i) {
        std::string t = rnd() & 1 ? "-" : "";
        const int nd = int(rnd() % 20) + 1, dot = int(rnd() % (nd + 2));
        for (int k = 0; k < nd; ++k) { if (k == dot) t += '.'; t += char('0' + rnd() % 10); }
        if (rnd() % 3 == 0) t += "e" + std::to_string(int(rnd() % 90) - 45);
        check(t);
    }
    // exact float midpoints and their neighbours, written with 17 digits
    for (int i = 0; i < 100000; ++i) {
        std::uint32_t fb = std::uint32_t(rnd()) & 0x7f7fffff;
        float f;
        std::memcpy(&f, &fb, 4);
        if (!std::isfinite(f) || f == 0) continue;
        const double mid = (double(f) + double(std::nextafter(f, 1e30f))) / 2;
        for (double v : {mid, std::nextafter(mid, 0.0), std::nextafter(mid, 1e300)}) {
            char buf[40];
            std::snprintf(buf, sizeof buf, "%.17g", v);
            check(buf);
        }
    }
    CHECK(mismatches == 0);
    CHECK(checked > 600000);
}

TEST_CASE("crc32_combine joins the CRCs of two parts at every split") {
    std::vector<std::uint8_t> d(3000);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = std::uint8_t(i * 131 + (i >> 3));
    const std::uint32_t whole = codec::crc32(d.data(), d.size());
    bool all = true;
    for (std::size_t cut = 0; cut <= d.size(); cut += 7)
        all &= codec::crc32_combine(codec::crc32(d.data(), cut), codec::crc32(d.data() + cut, d.size() - cut), d.size() - cut) == whole;
    CHECK(all);
    std::vector<std::uint8_t> big(9 << 20);
    for (std::size_t i = 0; i < big.size(); ++i) big[i] = std::uint8_t(i ^ (i >> 11));
    CHECK(codec::crc32_parallel(big.data(), big.size(), 6) == codec::crc32(big.data(), big.size()));
}

TEST_CASE("inflate_parallel equals inflate on large real streams, also when corrupted") {
    // every DEFLATE entry of at least 4 MB in the corpus archives (3MF, 3DXML)
    namespace fs = std::filesystem;
    std::vector<std::pair<std::vector<std::uint8_t>, std::size_t>> streams;  // raw stream, declared size
    for (const char* root : {CADCORE_SAMPLES_DIR}) {
        std::error_code ec;
        for (auto it = fs::recursive_directory_iterator(root, ec); it != fs::recursive_directory_iterator(); it.increment(ec)) {
            if (ec || !it->is_regular_file()) continue;
            const std::string ext = it->path().extension().string();
            if (ext != ".3mf" && ext != ".3dxml") continue;
            std::ifstream f(it->path(), std::ios::binary);
            std::vector<std::uint8_t> z((std::istreambuf_iterator<char>(f)), {});
            // local file headers: PK 03 04 ... method 8, compressed size, name, extra
            for (std::size_t p = 0; p + 30 <= z.size(); ++p) {
                if (!(z[p] == 'P' && z[p + 1] == 'K' && z[p + 2] == 3 && z[p + 3] == 4)) continue;
                auto u16 = [&](std::size_t o) { return std::size_t(z[o] | z[o + 1] << 8); };
                auto u32 = [&](std::size_t o) { return std::size_t(u16(o) | u16(o + 2) << 16); };
                const std::size_t method = u16(p + 8), csize = u32(p + 18), usize = u32(p + 22), data = p + 30 + u16(p + 26) + u16(p + 28);
                if (method == 8 && csize >= (4u << 20) && data + csize <= z.size())
                    streams.push_back({std::vector<std::uint8_t>(z.begin() + std::ptrdiff_t(data), z.begin() + std::ptrdiff_t(data + csize)), usize});
            }
        }
    }
    REQUIRE(!streams.empty());
    auto outcome = [](auto&& fn) -> std::pair<Buffer, std::string> {
        try { return {fn(), ""}; } catch (const DecodeError& e) { return {{}, e.what()}; }
    };
    std::uint64_t seed = 12345;
    auto rnd = [&] { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return seed; };
    int compared = 0, same = 0;
    for (auto& [raw, usize] : streams) {
        for (int variant = 0; variant < 4; ++variant) {
            std::vector<std::uint8_t> r = raw;
            if (variant > 0) r[rnd() % r.size()] ^= std::uint8_t(1u << (rnd() % 8));  // one flipped bit
            const ByteView v(r.data(), r.size());
            const auto seq = outcome([&] { return codec::inflate(v, usize, usize); });
            const auto par = outcome([&] { return codec::inflate_parallel(v, usize, usize, 6); });
            ++compared;
            same += seq == par;
        }
    }
    CAPTURE(compared);
    CHECK(same == compared);
}
