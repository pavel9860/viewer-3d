#include "formats/catia/spec.hpp"

#include <limits>
#include <memory>

#include "formats/catia/cursor.hpp"
#include "formats/catia/store.hpp"

namespace cadcore::formats::catia {
namespace {

constexpr std::uint32_t kEmpty = std::numeric_limits<std::uint32_t>::max();  // the empty name, numbered last

std::uint32_t le32(const std::uint8_t* p) {
    std::uint32_t x;
    std::memcpy(&x, p, 4);
    return x;
}

// Names written inline on first use and numbered in order; references as in OSM pools.
struct Names {
    std::vector<std::string>& table;
    std::uint32_t operator()(Cur& c) const {
        const std::uint8_t x = c.u8();
        if (x == 0x32) {
            const std::uint32_t n = c.u32();
            if (n == 0 || n > table.size()) bad("name reference outside the table", c.base + c.k - 4);
            return n - 1;
        }
        if (x >= 0x34 && x <= 0xfd) {
            if (x - 0x34u >= table.size()) bad("name reference outside the table", c.base + c.k - 1);
            return x - 0x34u;
        }
        if (x > 0x31) bad("invalid name token", c.base + c.k - 1);
        const std::size_t n = x == 0 ? c.u32() : x - 1u;
        if (n > c.v.size - c.k) bad("name past the container end", c.base + c.k);
        if (n == 0) return kEmpty;
        table.emplace_back(reinterpret_cast<const char*>(c.v.data + c.k), n);
        c.k += n;
        return std::uint32_t(table.size() - 1);
    }
};

// Fields that differ between SPEC versions (determined on versions 5..21 of the CAA samples).
struct Layout {
    int stamps;         // compacts after the handle (time stamps): 1 up to version 12, then 2
    bool kind2;         // a second kind byte (version 19+)
    bool startup_names; // startup container names after the id (version 10+)
    bool chunk_end;     // fe closing each chunk (version 10+)
    bool header_fe;     // fe in the container header (version 19+)
    bool param_count;   // a third compact in parameter values (version 13+)
    explicit Layout(int v)
        : stamps(v <= 12 ? 1 : 2), kind2(v >= 19), startup_names(v >= 10), chunk_end(v >= 10), header_fe(v >= 19), param_count(v >= 13) {}
};

// `first`: the container's first object (its catalog manager), whose chunk list has no closing fe
// from version 10 on (found in all 79 containers of versions 10..21 of the corpus; OSM catalog
// records also lack it).
void object(Cur& c, const Layout& L, Names& name, Object& o, bool first) {
    o.offset = c.base + c.k;
    o.flag = c.u8();
    if (o.flag != 0 && o.flag != 2) bad("invalid object flag", o.offset);
    o.name = name(c);
    o.cls = kEmpty;  // set by the relation block when the class is defined there
    for (int i = 0; i <= L.stamps; ++i) c.compact();  // handle, time stamps
    const std::uint8_t id_form = c.u8();
    if (id_form == 0xfe) {
        c.expect(0xf6, "guid marker");
        if (c.k + 16 > c.v.size) bad("object guid past the container end", c.base + c.k);
        c.k += 16;
    } else if (id_form == 0x80) {
        if (c.k + 5 > c.v.size) bad("object id past the container end", c.base + c.k);
        c.k += 5;
    } else {
        bad("invalid object identifier form", c.base + c.k - 1);
    }
    o.kind[0] = c.u8();
    if (L.kind2) o.kind[1] = c.u8();
    if (o.kind[0] & 1) c.compact();
    c.expect(0xea, "object id marker");
    o.id = c.u32();
    if (L.startup_names) {
        const std::uint64_t n = c.compact();
        if (n > c.v.size) bad("container name count exceeds the container", c.base + c.k);
        for (std::uint64_t i = 0; i < n; ++i) name(c);  // container of the startup
    }
    if (c.peek() == 0xe5) {
        o.startup.blob = c.blob();
        o.startup.has_blob = true;
        o.startup.number = c.compact();
    } else if (c.peek() == 0xfd) {
        ++c.k;
        o.startup.number = c.compact();
        o.startup_by_index = true;
    } else {
        o.startup.number = c.compact();
    }
    o.record_offset = c.base + c.k;
    while (c.peek() != 0xfe) {
        if (first && L.chunk_end && !o.chunks.empty() && (c.peek() == 0x00 || c.peek() == 0x02 || c.peek() == 0xde)) break;  // next object / relations
        Chunk ch;
        ch.attr = name(c);
        ch.head = c.compact();
        const std::uint8_t t = c.u8();
        if (t != 0xe8 && t != 0xe9) bad("invalid type descriptor", c.base + c.k - 1);
        ch.flags = c.u8();
        ch.code = c.u8();
        if (t == 0xe9) ch.extra = c.u8();
        ch.qualifier = name(c);
        ch.value = chunk_value(c, ch.flags, ch.code, name, false);
        if (ch.extra & 0x01) literal_block(c, name);
        if (L.chunk_end) c.expect(0xfe, "chunk end");
        o.chunks.push_back(std::move(ch));
    }
    if (c.peek() == 0xfe) ++c.k;
    if (o.flag == 2) {  // parameter value
        c.compact();
        const std::uint64_t type = c.compact();
        if (L.param_count) c.compact();
        if (c.peek() == 0xfe) ++c.k;
        else name(c);  // unit
        if (type == 0x07 || type == 0x15) c.number();
        else if (type == 0x12) name(c);
        else if (type == 0x11 || type == 0x1c) c.compact();
        else bad("parameter of unknown type " + std::to_string(type), c.base + c.k);
    }
}

constexpr std::uint64_t kNoObject = std::numeric_limits<std::uint64_t>::max();

// compact (object id) | e5 <blob> <compact> (link) | 7f <a> <n> <id>*n <index>*(n+1) (a path into
// another container, not followed here)
std::uint64_t relation_value(Cur& c, std::vector<LinkValue>& links) {
    if (c.peek() == 0x7f) {
        ++c.k;
        c.compact();
        const std::uint64_t n = c.compact();
        if (n > 64) bad("pointer path too long", c.base + c.k);
        for (std::uint64_t i = 0; i < 2 * n + 1; ++i) c.compact();
        return kNoObject;
    }
    if (c.peek() != 0xe5) return c.compact();
    LinkValue l;
    l.blob = c.blob();
    l.has_blob = true;
    l.number = c.compact();
    links.push_back(l);
    return kNoObject;
}

void relation_block(Cur& c, Names& name, Object& o) {
    const std::uint8_t head = c.u8();
    if ((head & 7) != 0 && c.peek() == 0xfe) {  // class defined here: fe (fe | link) <name> fe (as in OSM before version 35)
        ++c.k;
        if (c.peek() == 0xfe) {
            ++c.k;
        } else {
            c.blob();
            c.compact();
        }
        o.cls = name(c);
        c.expect(0xfe, "class definition end");
    }
    while (c.peek() != 0xfe) {
        HeaderToken t;
        if (c.peek() == 0xe5) {
            t.is_blob = true;
            t.blob = c.blob();
        } else if (c.peek() == 0xfd) {  // startup by index
            ++c.k;
            t.is_index = true;
            t.number = c.compact();
        } else {
            t.number = c.compact();
        }
        o.header.push_back(t);
    }
    ++c.k;
    while (c.peek() != 0xfe) {
        Relation r;
        r.attr = name(c);
        const Chunk* chunk = nullptr;
        for (const Chunk& ch : o.chunks)
            if (ch.attr == r.attr) chunk = &ch;
        if (!chunk) bad("relation for an attribute the object does not declare", c.base + c.k);
        auto one = [&] {
            const std::uint64_t v = relation_value(c, r.links);
            if (v != kNoObject) r.objects.push_back(std::uint32_t(v));
        };
        if (chunk->array()) {  // count, an index into the list (< count; not interpreted), elements
            const std::uint64_t n = c.compact();
            if (n > c.v.size) bad("relation count exceeds the container", c.base + c.k);
            if (n && c.compact() >= n) bad("relation list index outside the list", c.base + c.k);
            for (std::uint64_t i = 0; i < n; ++i) one();
        } else {
            one();
        }
        o.rels.push_back(std::move(r));
    }
    ++c.k;
    o.preparsed = true;
}

void property_tables(Cur& c, Names& name) {
    while (c.peek() != 0xfe) {
        const std::uint64_t n = c.compact();
        if (n > c.v.size) bad("property count exceeds the container", c.base + c.k);
        for (std::uint64_t i = 0; i < n; ++i) {
            name(c);
            switch (c.compact()) {
                case 2:
                case 3: c.compact(); break;
                case 0:
                case 14: break;
                case 13: {
                    const std::uint8_t len = c.u8();
                    if (len < 0xe8 || len > 0xef) bad("invalid property octet length", c.base + c.k - 1);
                    c.u8();  // form (80, 81, 83, 84 seen)
                    if (len - 0xe7u > c.v.size - c.k) bad("property past the container end", c.base + c.k);
                    c.k += len - 0xe7u;
                    break;
                }
                default: bad("property of unknown type", c.base + c.k - 1);
            }
        }
    }
    ++c.k;
}

// One SPEC container (the bytes after its item header).
Stream container(ByteView v, std::size_t base, std::shared_ptr<const std::vector<std::uint8_t>> storage) {
    Stream s;
    s.offset = base;
    s.storage = std::move(storage);
    Cur c{v, 4, base};
    Names name{s.names};
    // "SPEC" <version> <compact> [fe from version 19]
    // "OSMB" <version> <compact> <name>* fe          "OSMN" <name>* <version> <compact> fe
    const bool osmb = std::memcmp(v.data, "OSMB", 4) == 0, osmn = std::memcmp(v.data, "OSMN", 4) == 0;
    if (osmn)
        while (c.peek() < 0x81) name(c);
    const std::uint64_t version = c.compact();
    if (version == 0 || version > 1000) bad("invalid SPEC version", base + 4);
    s.spec_version = int(version);
    c.compact();
    const Layout layout(s.spec_version);
    if (osmb)
        while (c.peek() != 0xfe) name(c);
    if (osmb || osmn || layout.header_fe) c.expect(0xfe, "container header end");
    const std::uint64_t highest = c.compact();
    while (c.peek() != 0xde) {
        Object o;
        object(c, layout, name, o, s.objects.empty());
        if (!s.objects.empty() && o.id <= s.objects.back().id) bad("object ids do not increase", o.offset);
        if (o.id > highest) bad("object id above the container's highest id", o.offset);
        s.objects.push_back(std::move(o));
    }
    if (s.objects.empty()) bad("SPEC container without objects", base);
    ++c.k;
    for (Object& o : s.objects) {
        o.relation_offset = c.base + c.k;
        relation_block(c, name, o);
    }
    property_tables(c, name);
    if (v.size - c.k >= 8) bad("bytes after the SPEC container", c.base + c.k);  // < 8: padding
    // the empty name gets the last number
    const std::uint32_t empty = std::uint32_t(s.names.size());
    s.names.emplace_back();
    auto fix = [&](std::uint32_t& n) { if (n == kEmpty) n = empty; };
    for (Object& o : s.objects) {
        fix(o.name);
        fix(o.cls);
        for (Chunk& ch : o.chunks) {
            fix(ch.attr);
            fix(ch.qualifier);
            if (ch.value.kind == Value::Kind::Name) {
                std::uint32_t n = std::uint32_t(ch.value.i);
                fix(n);
                ch.value.i = n;
            }
            for (auto& n : ch.value.names) fix(n);
        }
        for (Relation& r : o.rels) fix(r.attr);
    }
    s.records_parsed = true;
    return s;
}

}  // namespace

std::vector<Stream> parse_spec_streams(ByteView file) {
    std::vector<Stream> out;
    if (file.size < 8 || std::memcmp(file.data, "V5_CFV", 6) != 0) return out;
    const Store store = read_store(file);
    const StreamInfo* info = store.find("Data");
    if (!info) return out;
    auto bytes = std::make_shared<const std::vector<std::uint8_t>>(read_stream(file, store, *info));
    const std::uint8_t* d = bytes->data();
    const std::size_t n = bytes->size();
    std::size_t base = 0;  // offsets in messages are positions in the Data stream
    for (std::size_t a = 0; a < n;) {
        if (a + 16 > n) bad("Data section header past the stream end", base + a);
        const std::uint64_t end = a + 8 * std::uint64_t(le32(d + a + 4));
        if (end <= a + 8 || end > n) bad("Data section length out of range", base + a);
        if ((le32(d + a + 8) & 0xffff) != 1) bad("Data section without its container header", base + a + 8);
        for (std::size_t q = a + 16; q < end;) {
            if (q + 8 > end) bad("Data item header past its section", base + q);
            const std::uint16_t tag = std::uint16_t(d[q] | d[q + 1] << 8);
            const std::uint64_t len = 8 * std::uint64_t(le32(d + q + 4));
            if (len < 8 || q + len > end) bad("Data item length out of range", base + q);
            const bool spec = std::memcmp(d + q + 8, "SPEC", 4) == 0 || std::memcmp(d + q + 8, "OSMB", 4) == 0 || std::memcmp(d + q + 8, "OSMN", 4) == 0;
            if ((tag == 4 || tag == 5) && len >= 12 && spec)
                out.push_back(container(ByteView(d + q + 8, std::size_t(len - 8)), base + q + 8, bytes));
            q += std::size_t(len);
        }
        a = std::size_t(end);
    }
    return out;
}

}  // namespace cadcore::formats::catia
