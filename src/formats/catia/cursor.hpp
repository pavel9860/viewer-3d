// Cursor over the token encoding shared by CATIA V5 object streams (OSM records and the older
// SPEC containers): compact numbers, number tokens, length-prefixed blobs and GUIDs.
#pragma once

#include <cstdint>
#include <cstring>
#include <string>

#include "cadcore/bytes.hpp"
#include "cadcore/error.hpp"
#include "formats/catia/osm.hpp"

namespace cadcore::formats::catia {

[[noreturn]] inline void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CATIA object stream: " + what, at); }

// Runs `f`; a DecodeError from it is rethrown with "[<context()>]" appended (context() is only
// called on error, so it may format freely).
template <class F, class Context>
decltype(auto) in_context(F&& f, Context&& context) {
    try {
        return f();
    } catch (const DecodeError& e) {
        fail(e.code(), e.message() + " [" + context() + "]", e.offset());
    }
}

// Cursor over one record body; `base` maps positions back to file offsets for messages.
struct Cur {
    ByteView v;
    std::size_t k = 0, base = 0;
    // Byte order of the writer's platform for u32 / f64 fields (compact 80 form, blob lengths,
    // long names, numbers): little-endian except in files written on big-endian UNIX systems.
    bool big = false;

    bool end() const { return k >= v.size; }
    std::uint8_t peek(std::size_t ahead = 0) const {
        if (k + ahead >= v.size) bad("record ends early", base + k);
        return v[k + ahead];
    }
    std::uint8_t u8() { const std::uint8_t x = peek(); ++k; return x; }
    std::uint32_t u32() {
        if (k + 4 > v.size) bad("record ends early", base + k);
        k += 4;
        return load_u32(v.data + k - 4, big);
    }
    // A u32 that is little-endian whatever the writer's platform (the object id after ea).
    std::uint32_t le32() {
        if (k + 4 > v.size) bad("record ends early", base + k);
        k += 4;
        return load_u32(v.data + k - 4, false);
    }
    double f64() {
        if (k + 8 > v.size) bad("record ends early", base + k);
        std::uint64_t u = 0;
        for (int i = 0; i < 8; ++i) u |= std::uint64_t(v.data[k + std::size_t(i)]) << (8 * (big ? 7 - i : i));
        k += 8;
        double x;
        std::memcpy(&x, &u, 8);
        return x;
    }
    void expect(std::uint8_t b, const char* what) {
        if (u8() != b) bad(std::string("expected ") + what, base + k - 1);
    }

    // Compact unsigned number: 80 u32 | 81..d0 (x-0x81) | d1..ef nn ((x-0xd1)*256+nn).
    std::uint64_t compact() {
        const std::uint8_t x = u8();
        if (x == 0x80) return u32();
        if (x >= 0x81 && x <= 0xd0) return x - 0x81u;
        if (x >= 0xd1 && x <= 0xef) return std::uint64_t(x - 0xd1) * 256 + u8();
        bad("invalid compact number", base + k - 1);
    }
    double number() {
        switch (u8()) {
            case 0xe6: return f64();
            case 0xe7: return 1.0;
            case 0xe8: return 0.0;
            case 0xe9: return -1.0;
            default: bad("invalid number token", base + k - 1);
        }
    }
    ByteView guid() {  // f6 <16 bytes>
        expect(0xf6, "guid marker");
        if (k + 16 > v.size) bad("guid past record end", base + k);
        const ByteView g = v.sub(k, 16);
        k += 16;
        return g;
    }
    ByteView blob() {
        expect(0xe5, "blob marker");
        const std::uint32_t n = u32();
        if (n > v.size - k) bad("blob past record end", base + k);
        const ByteView b = v.sub(k, n);
        k += n;
        return b;
    }
};

// List: compact count, then (when not empty) the separator 81 and `count` elements.
template <class F>
void counted(Cur& c, F&& one) {
    const std::uint64_t n = c.compact();
    if (n > c.v.size) bad("list count exceeds the record", c.base + c.k);
    if (n == 0) return;
    c.expect(0x81, "list separator");
    for (std::uint64_t i = 0; i < n; ++i) one();
}

// Chunk values and parameter blocks; `name(c)` reads a name token and returns its table index.
inline LinkValue link_value(Cur& c) {
    LinkValue l;
    l.a = c.compact();
    const std::uint8_t x = c.peek();
    if (x == 0xe5) {
        l.blob = c.blob();
        l.has_blob = true;
    } else if (x >= 0xe0 && x <= 0xef) {  // 3-byte form, seen only in catalog links
        l.number = (std::uint64_t(c.u8()) << 16);
        l.number |= std::uint64_t(c.u8()) << 8;
        l.number |= c.u8();
    } else {
        l.number = c.compact();
    }
    return l;
}

// `element_values`: link lists carry a compact after each link (OSM records; SPEC has none).
template <class Name>
Value chunk_value(Cur& c, std::uint8_t flags, std::uint8_t code, Name&& name, bool element_values) {
    Value v;
    if (c.peek() == 0xfe) {  // no inline value (kept in the relations, or none)
        ++c.k;
        return v;
    }
    const bool array = flags & 0x20;
    switch (code) {
        case 0x12:
        case 0x1d:
            if (array) {  // list of names: count, separator, name references
                v.kind = Value::Kind::Names;
                counted(c, [&] { v.names.push_back(name(c)); });
                return v;
            }
            v.kind = Value::Kind::Name;
            v.i = name(c);
            return v;
        case 0x1c:
        case 0x08:
            if (array) {
                v.kind = Value::Kind::Ints;
                counted(c, [&] { v.ints.push_back(c.compact()); });
                return v;
            }
            v.kind = Value::Kind::Int;
            v.i = c.compact();
            return v;
        case 0x07:
            if (array) {
                v.kind = Value::Kind::Doubles;
                counted(c, [&] { v.d.push_back(c.number()); });
            } else {
                v.kind = Value::Kind::Double;
                v.d.push_back(c.number());
            }
            return v;
        case 0x0a:  // octets: count, 81, count bytes
            if (array) {
                v.kind = Value::Kind::Octets;
                const std::uint64_t n = c.compact();
                if (n > c.v.size - c.k) bad("octet array past the record", c.base + c.k);
                if (n) {
                    c.expect(0x81, "list separator");
                    if (n > c.v.size - c.k) bad("octet array past the record", c.base + c.k);
                    v.octets = c.v.sub(c.k, std::size_t(n));
                    c.k += std::size_t(n);
                }
                return v;
            }
            bad("scalar octet value", c.base + c.k);
        case 0x18:
            if (array) {
                v.kind = Value::Kind::Links;
                counted(c, [&] {
                    v.links.push_back(link_value(c));
                    if (element_values) c.compact();  // per-element value
                });
            } else {
                v.kind = Value::Kind::Link;
                v.links.push_back(link_value(c));
            }
            return v;
        default:
            bad("value of unknown type 0x" + std::to_string(code), c.base + c.k);
    }
}

// Parameter (literal) block, present when the e9 descriptor's extra byte has bit 0x01:
//   88 <n> <n> <n> <n> <m> [<name> <name> when m == 0xfffffffe] <count> <number>*count <n>
// (seen on user properties and knowledgeware parameters; the fields are not interpreted).
template <class Name>
void literal_block(Cur& c, Name&& name) {
    c.expect(0x88, "parameter block marker");
    for (int i = 0; i < 4; ++i) c.compact();
    if (c.compact() == 0xfffffffeu) {
        name(c);
        name(c);
    }
    const std::uint64_t n = c.compact();
    if (n > c.v.size) bad("parameter block count exceeds the record", c.base + c.k);
    for (std::uint64_t i = 0; i < n; ++i) c.number();
    c.compact();
}

}  // namespace cadcore::formats::catia
