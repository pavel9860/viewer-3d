#include "formats/catia/data.hpp"

#include <cstring>

#include "cadcore/error.hpp"
#include "formats/catia/store.hpp"

namespace cadcore::formats::catia {
namespace {

[[noreturn]] void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CATIA Data stream: " + what, at); }


// NUL-terminated string at `at` within `b`, or the rest of `b`.
std::string cstr(ByteView b, std::size_t& at) {
    std::size_t e = at;
    while (e < b.size && b[e] != 0) ++e;
    std::string s(reinterpret_cast<const char*>(b.data + at), e - at);
    at = e + 1;
    return s;
}

}  // namespace

std::optional<DataStream> read_data_stream(ByteView file) {
    if (file.size < 8 || std::memcmp(file.data, "V5_CFV", 6) != 0) return std::nullopt;
    const Store store = read_store(file);
    const StreamInfo* info = store.find("Data");
    if (!info) return std::nullopt;
    DataStream out;
    out.bytes = std::make_shared<const std::vector<std::uint8_t>>(read_stream(file, store, *info));
    const std::uint8_t* d = out.bytes->data();
    const std::size_t n = out.bytes->size();
    if (n == 0) return out;
    if (n < 16) bad("shorter than a section header", 0);
    out.big = d[8] == 0 && d[9] == 1;
    const bool big = out.big;
    for (std::size_t a = 0; a < n;) {
        if (a + 16 > n) bad("section header past the stream end", a);
        const std::uint64_t end = a + 8 * std::uint64_t(load_u32(d + a + 4, big));
        if (end < a + 16 || end > n) bad("section length out of range", a);
        if (load_u16(d + a + 8, big) != 1) bad("section without its container header", a + 8);
        DataSection s;
        s.offset = a;
        s.kind = load_u16(d + a + 10, big);
        s.index = load_u32(d + a + 12, big);
        for (std::size_t q = a + 16; q < end;) {
            if (q + 8 > end) bad("item header past its section", q);
            const std::uint64_t len = 8 * std::uint64_t(load_u32(d + q + 4, big));
            if (len < 8 || q + len > end) bad("item length out of range", q);
            DataItem it{load_u16(d + q, big), q + 8, ByteView(d + q + 8, std::size_t(len - 8))};
            if (it.tag == 2 && s.cls.empty()) {
                std::size_t k = 0;
                s.cls = cstr(it.body, k);
                if (k < it.body.size) s.super = cstr(it.body, k);
            }
            s.items.push_back(it);
            q += std::size_t(len);
        }
        out.sections.push_back(std::move(s));
        a = std::size_t(end);
    }
    return out;
}

}  // namespace cadcore::formats::catia
