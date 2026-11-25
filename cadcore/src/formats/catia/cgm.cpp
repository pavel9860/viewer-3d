#include "formats/catia/cgm.hpp"

#include <cstdio>
#include <cstring>

#include "cadcore/error.hpp"
#include "formats/catia/data.hpp"
#include "formats/catia/store.hpp"

namespace cadcore::formats::catia::cgm {
namespace {

[[noreturn]] void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CGM stream: " + what, at); }

std::string hex_name(const std::uint8_t* g) {
    std::string n = "Cont";
    char buf[16];
    for (int i = 0; i < 4; ++i) {
        std::snprintf(buf, sizeof buf, "_%x", load_u32(g + 4 * i, true));  // the GUID bytes in file order
        n += buf;
    }
    return n;
}

// varint dialect: records from the directory at `dir` (e5 0d 06 01 00 u32 level u32 size).
void frame_varint(ByteView s, std::size_t dir, bool big, Stream& out) {
    out.dialect = Dialect::Varint;
    out.level = load_u32(s.data + dir + 5, big);
    const std::size_t size = load_u32(s.data + dir + 9, big);
    std::size_t k = dir + 13;
    const std::size_t end = k + size;
    if (end > s.size) bad("directory size past the stream end", dir);
    while (k < end) {
        if (k + 13 > end || s[k] != 0xe5 || (s[k + 1] != 0x0d && s[k + 1] != 0x0f)) bad("record does not start with e5 0d|0f", k);
        const std::size_t hl = s[k + 1] == 0x0d ? 13 : 15;
        Record r;
        r.offset = k;
        r.flags = s[k + 2];
        r.type = load_u16(s.data + k + 3, big);
        const std::size_t len = load_u32(s.data + k + 5, big);
        r.id = load_u32(s.data + k + 9, big);
        if (k + hl + len > end) bad("record past the directory size", k);
        r.body = ByteView(s.data + k + hl, len);
        out.records.push_back(r);
        k += hl + len;
    }
    out.end = end;
}

// Class number of a compact type code; 0 when not established. Mapped from the reference graph
// (which codes each code's records reference, at which positions) and the payload lengths of the
// 65 2009Fall streams, whose reference lists still use the varint encoding (FINDINGS.md).
std::uint16_t compact_class(std::uint8_t code) {
    switch (code) {
        case 2: return 10;    case 4: return 12;    case 5: return 13;    case 6: return 14;    case 14: return 100;
        case 24: return 150;  case 25: return 151;  case 26: return 152;  case 32: return 160;
        case 35: return 191;  case 36: return 192;  case 37: return 193;
        case 29: return 155;
        case 39: return 200;  case 40: return 201;  case 41: return 202;  case 42: return 203;  case 43: return 204;
        case 44: return 210;  case 45: return 211;  case 46: return 212;  case 48: return 214;  case 50: return 216;
        case 56: return 242;  case 57: return 243;  case 59: return 244;  case 64: return 251;
        case 52: return 231;  case 55: return 241;
        case 91: return 500;  case 92: return 502;  case 93: return 510;  case 94: return 511;
        case 95: return 512;  case 96: return 513;  case 97: return 520;  case 98: return 521;
        case 100: return 524; case 101: return 525;
        default: return 0;
    }
}

// compact dialect: directory a8 06 7a u32 level u32 size at `dir`.
void frame_compact(ByteView s, std::size_t dir, Stream& out) {
    out.dialect = Dialect::Compact;
    out.level = load_u32(s.data + dir + 3, false);
    out.varint_bodies = false;
    const std::size_t size = load_u32(s.data + dir + 7, false);
    std::size_t k = dir + 11;
    const std::size_t end = k + size;
    if (end > s.size) bad("directory size past the stream end", dir);
    std::uint32_t prev = 0;
    while (k < end) {
        const std::uint8_t m = s[k];
        const unsigned hi = m >> 4, lo = m & 15;
        const unsigned lw = hi == 0xb ? 1 : hi == 0xa ? 4 : 0;
        if (!lw || lo <= lw || k + 3 + lo > end) bad("unknown compact record marker", k);
        const unsigned iw = lo - lw;
        Record r;
        r.offset = k;
        r.flags = s[k + 1];
        r.code = s[k + 2];
        r.type = compact_class(r.code);
        std::size_t len = 0;
        for (unsigned i = 0; i < lw; ++i) len |= std::size_t(s[k + 3 + i]) << (8 * i);
        if (out.level <= 82) {
            // the id field holds the signed delta from the previous record's id (compact integer)
            const std::uint8_t x = s[k + 3 + lw];
            std::int64_t d;
            if ((x & 3) == 1) d = x >> 2;
            else if ((x & 3) == 3) d = -std::int64_t(x >> 2) - 1;
            else {
                const unsigned n = x >> 2 & 3;  // bytes that follow; x's high nibble above them
                if (n + 1 > iw) bad("compact id field", k);
                std::uint64_t v = 0;
                for (unsigned i = 0; i < n; ++i) v |= std::uint64_t(s[k + 4 + lw + i]) << (8 * i);
                v |= std::uint64_t(x >> 4) << (8 * n);
                d = (x & 3) == 0 ? std::int64_t(v) : -std::int64_t(v) - 1;
            }
            prev = std::uint32_t(std::int64_t(prev) + d);
            r.id = prev;
        } else {
            for (unsigned i = 0; i < iw && i < 4; ++i) r.id |= std::uint32_t(s[k + 3 + lw + i]) << (8 * i);
        }
        const std::size_t body = k + 3 + lo;
        if (body + len > end) bad("record past the directory size", k);
        r.body = ByteView(s.data + body, len);
        out.records.push_back(r);
        k = body + len;
    }
    out.end = end;
}

std::size_t find_from(ByteView s, std::size_t from, const char* needle, std::size_t n) {
    const std::size_t at = find_bytes(s.data + from, s.size - from, reinterpret_cast<const std::uint8_t*>(needle), n);
    return at == s.size - from ? ByteView::npos : from + at;
}

}  // namespace

const char* dialect_name(Dialect d) {
    switch (d) {
        case Dialect::Varint: return "varint";
        case Dialect::Compact: return "compact";
        case Dialect::Fixed: return "fixed";
    }
    return "?";
}

std::optional<std::vector<std::uint8_t>> find_stream(ByteView file) {
    const auto data = read_data_stream(file);
    if (!data) return std::nullopt;
    for (const DataSection& sec : data->sections) {
        if (sec.cls != "CGMGeom") continue;
        for (const DataItem& it : sec.items) {
            if (it.tag != 3 || it.body.size < 16) continue;
            const std::string name = hex_name(it.body.data);
            const Store store = read_store(file);
            if (const StreamInfo* st = store.find(name)) return read_stream(file, store, *st);
        }
    }
    return std::nullopt;
}

std::vector<std::vector<std::uint32_t>> see_body_lists(ByteView file) {
    std::vector<std::vector<std::uint32_t>> out;
    const auto data = read_data_stream(file);
    if (!data) return out;
    for (const DataSection& sec : data->sections) {
        if (sec.cls != "CATSeeBodyCont") continue;
        for (const DataItem& it : sec.items) {
            if (it.tag != 3 || it.body.size < 16) continue;
            const Store store = read_store(file);
            const StreamInfo* st = store.find(hex_name(it.body.data));
            if (!st) continue;
            const auto b = read_stream(file, store, *st);
            std::vector<std::vector<std::uint32_t>> lists;
            for (std::size_t k = 0; k < b.size();) {
                if (k + 4 > b.size()) bad("CATSeeBodyCont: list count past the end", k);
                const std::size_t n = load_u32(b.data() + k, data->big);
                k += 4;
                if (n > (b.size() - k) / 4) bad("CATSeeBodyCont: list past the end", k);
                std::vector<std::uint32_t> l(n);
                for (auto& x : l) { x = load_u32(b.data() + k, data->big); k += 4; }
                lists.push_back(std::move(l));
            }
            return lists;
        }
    }
    return out;
}

std::uint32_t varint(ByteView b, std::size_t& k) {
    if (k >= b.size) fail(ErrorCode::Truncated, "CGM varint past the end", k);
    const std::uint8_t t = b[k++];
    if (t >= 0x80) return t & 0x7f;
    std::uint32_t v = 0;
    for (int j = 0; j < 4; ++j)
        if (t >> (3 + j) & 1) {
            if (k >= b.size) fail(ErrorCode::Truncated, "CGM varint past the end", k);
            v |= std::uint32_t(b[k++]) << (8 * j);
        }
    return v;
}

Stream frame(ByteView s) {
    if (s.size < 9 || s[0] > 1) bad("no stream header", 0);
    Stream out;
    out.big = s[0] == 0;
    const std::size_t n = load_u32(s.data + 5, out.big);
    if (9 + n > s.size) bad("version string past the end", 5);
    out.version.assign(reinterpret_cast<const char*>(s.data + 9), n);
    // The directory follows the version string(s) and the root ids; the earliest marker wins.
    const std::size_t from = 9 + n;
    const std::size_t c = find_from(s, from, "\xa8\x06\x7a", 3);
    const std::size_t v = find_from(s, from, "\xe5\x0d\x06\x01\x00", 5);
    const std::size_t vb = find_from(s, from, "\xe5\x0d\x06\x00\x01", 5);  // big-endian type
    const std::size_t vv = std::min(v, vb);
    if (c != ByteView::npos && c < vv) {
        // feature bytes: after the second version string (u32 length + text), up to the directory
        if (from + 4 <= c) {
            const std::size_t n2 = load_u32(s.data + from, out.big);
            if (from + 4 + n2 <= c) out.features.assign(s.data + from + 4 + n2, s.data + c);
        }
        frame_compact(s, c, out);
        return out;
    }
    if (vv != ByteView::npos) {
        std::size_t k = from;
        if (k + 2 <= vv && s[k] == 3 && s[k + 1] == 0) k += 2;
        varint(s, k);
        if (k < vv && s[k] != 0x17) {
            const std::uint32_t n = varint(s, k);
            for (std::uint32_t i = 0; i < n && k < vv; ++i) out.roots.push_back(varint(s, k));
        }
        if (k >= vv || s[k] != 0x17) out.roots.clear();  // not the documented header form
        frame_varint(s, vv, out.big, out);
        return out;
    }
    if (s.find("\xb3", from, std::min(s.size, from + 64)) != ByteView::npos) {
        out.dialect = Dialect::Fixed;
        fail(ErrorCode::Unsupported, "CGM stream: fixed-width dialect (CXR4/5, CGA-V2R4) not framed", from);
    }
    bad("no record directory", from);
}

}  // namespace cadcore::formats::catia::cgm
