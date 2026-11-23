#include "formats/catia/osm.hpp"

#include "formats/catia/cursor.hpp"

#include <cstring>

#include "cadcore/error.hpp"

namespace cadcore::formats::catia {
namespace {

// Interned name reference: 34..fd (x-0x34) | 32 u32 (n-1).
std::uint32_t pool_name(Cur& c) {
    const std::uint8_t x = c.u8();
    if (x == 0x32) {
        const std::uint32_t n = c.u32();
        if (n == 0) bad("zero escaped name id", c.base + c.k - 4);
        return n - 1;
    }
    if (x >= 0x34 && x <= 0xfd) return x - 0x34u;
    bad("invalid name reference", c.base + c.k - 1);
}

// 7c 07 body: chunk* [fe]; returns chunks. The closing fe is absent on catalog records.
std::vector<Chunk> attribute_record(ByteView body, std::size_t base) {
    Cur c{body, 0, base};
    std::vector<Chunk> out;
    while (!c.end() && c.peek() != 0xfe) {
        Chunk ch;
        ch.attr = pool_name(c);
        ch.head = c.compact();
        const std::uint8_t t = c.u8();
        if (t != 0xe8 && t != 0xe9) bad("invalid type descriptor", base + c.k - 1);
        ch.flags = c.u8();
        ch.code = c.u8();
        if (t == 0xe9) ch.extra = c.u8();
        ch.qualifier = pool_name(c);  // usually pool entry 3 (the empty name)
        ch.value = chunk_value(c, ch.flags, ch.code, pool_name, true);
        if (ch.extra & 0x01) literal_block(c, pool_name);
        c.expect(0xfe, "chunk end");
        out.push_back(std::move(ch));
    }
    if (!c.end()) {
        ++c.k;  // record terminator
        if (!c.end()) bad("bytes after the attribute record", base + c.k);
    }
    return out;
}

std::uint32_t le32(ByteView f, std::size_t at) {
    if (at + 4 > f.size) fail(ErrorCode::Truncated, "CATIA object stream ends early", at);
    std::uint32_t x;
    std::memcpy(&x, f.data + at, 4);
    return x;
}

// Bounded tagged record: 7c <tag> <u32 len> (len counts from the 7c).
ByteView tagged(ByteView f, std::size_t at, std::uint8_t tag, std::size_t limit) {
    if (at + 6 > limit || f[at] != 0x7c || f[at + 1] != tag) bad("expected record 7c " + std::to_string(tag), at);
    const std::uint32_t len = le32(f, at + 2);
    if (len < 6 || len > limit - at) bad("record length out of range", at);
    return f.sub(at, len);
}

bool object_header_valid(ByteView f, std::size_t at) {
    if (at + 13 > f.size || f[at] != 0x7c || f[at + 1] != 0x05) return false;
    const std::uint32_t len = le32(f, at + 2);
    if (len < 20 || len > f.size - at) return false;
    if (f[at + 6] != 0 && f[at + 6] != 2) return false;
    if (f[at + 7] != 0x7c || f[at + 8] != 0x06) return false;
    const std::uint32_t hl = le32(f, at + 9);
    if (hl < 12 || 7 + std::size_t(hl) + 6 > len) return false;
    return f[at + 7 + hl] == 0x7c && f[at + 7 + hl + 1] == 0x07;
}

Stream parse_stream(ByteView f, std::size_t at) {
    Stream s;
    s.offset = at;
    std::size_t q = at + 6;
    while (q + 2 <= f.size && f[q] == 0x7c && f[q + 1] == 0x05) {
        const ByteView o = tagged(f, q, 0x05, f.size);
        Object ob;
        ob.offset = q;
        ob.flag = o[6];
        if (ob.flag == 8) {  // compact object form (outside product streams): class, name, opaque rest
            Cur h{o.sub(7, o.size - 7), 0, q + 7};
            ob.cls = pool_name(h);
            ob.name = pool_name(h);
            ob.refs = h.v.sub(h.k, h.v.size - h.k);
            s.objects.push_back(std::move(ob));
            q += o.size;
            continue;
        }
        if (ob.flag != 0 && ob.flag != 2 && ob.flag != 7) bad("invalid object flag", q + 6);
        const ByteView hdr = tagged(f, q + 7, 0x06, q + o.size);
        if (hdr.size < 12) bad("object header too short", q + 7);
        Cur h{hdr.sub(6, hdr.size - 6), 0, q + 13};
        if (h.peek() == 0x32 || h.peek() >= 0x34) {
            // Layout B: <class> <handle> <stamp> <stamp> (fe f6 <guid 16> | 80 <5 bytes>) <kind 2>
            // [<compact> when the kind is odd] ea <id> ...  (the object header of SPEC containers)
            ob.cls = pool_name(h);
            ob.name = ob.cls;
            for (int i = 0; i < 3; ++i) h.compact();
            const std::uint8_t id_form = h.u8();
            if (id_form == 0xfe) {
                h.expect(0xf6, "guid marker");
                if (h.k + 16 > h.v.size) bad("object guid past header end", h.base + h.k);
                h.k += 16;
            } else if (id_form == 0x80) {
                if (h.k + 5 > h.v.size) bad("object id past header end", h.base + h.k);
                h.k += 5;
            } else {
                bad("invalid object identifier form", h.base + h.k - 1);
            }
            ob.kind[0] = h.u8();
            ob.kind[1] = h.u8();
            if (ob.kind[0] & 1) h.compact();  // odd kinds carry a number (as in SPEC headers)
        } else {  // Layout A: <kind 2> <class> <name> ea <id> ...
            ob.kind[0] = h.u8();
            ob.kind[1] = h.u8();
            ob.cls = pool_name(h);
            ob.name = pool_name(h);
        }
        h.expect(0xea, "object id marker");
        ob.id = h.u32();
        ob.refs = h.v.sub(h.k, h.v.size - h.k);
        const std::size_t r = q + 7 + hdr.size;
        const ByteView rec = tagged(f, r, 0x07, q + o.size);
        ob.record = rec.sub(6, rec.size - 6);
        ob.record_offset = r + 6;
        ob.trailer = f.sub(r + rec.size, q + o.size - (r + rec.size));
        s.objects.push_back(std::move(ob));
        q += o.size;
    }
    if (s.objects.empty()) bad("stream without objects", at);

    // Relations: one 7c 09 block per object.
    if (q + 7 > f.size || f[q] != 0xde) bad("missing relations marker", q);
    ++q;
    if (f[q] != 0x7c || f[q + 1] != 0x08) bad("missing relations record", q);
    q += 6;
    for (Object& ob : s.objects) {
        ob.relation_block = tagged(f, q, 0x09, f.size);
        ob.relation_offset = q;
        q += ob.relation_block.size;
    }
    if (q + 2 > f.size || f[q] != 0x7c || f[q + 1] == 0x09) bad("more relation blocks than objects", q);

    // 7c 0b record, one byte, then the name table.
    const ByteView b0b = tagged(f, q, 0x0b, f.size);
    q += b0b.size + 1;
    const ByteView pool = tagged(f, q, 0x02, f.size);
    Cur p{pool.sub(6, pool.size - 6), 0, q + 6};
    const std::uint64_t count = p.compact();
    if (count > p.v.size) bad("name count exceeds the table", q + 6);
    s.names.reserve(std::size_t(count));
    for (std::uint64_t i = 0; i < count; ++i) {
        std::uint32_t len = p.u8();
        std::size_t n;
        if (len == 0) {
            n = p.u32();
        } else {
            n = len - 1;
        }
        if (n > p.v.size - p.k) bad("name past table end", p.base + p.k);
        s.names.emplace_back(reinterpret_cast<const char*>(p.v.data + p.k), n);
        p.k += n;
    }
    if (!p.end()) bad("bytes after the name table", p.base + p.k);

    for (const Object& ob : s.objects)
        if (ob.cls >= s.names.size() || ob.name >= s.names.size()) bad("object class/name outside the name table", ob.offset);
    return s;
}

[[noreturn]] void unsupported_object(const Object& o) {
    fail(ErrorCode::Unsupported, "CATIA object stream: compact (flag 8) object in a stream that must be read", o.offset);
}

}  // namespace

void parse_records(Stream& s) {
    if (s.records_parsed) return;
    for (Object& ob : s.objects) {
        if (ob.flag == 8) unsupported_object(ob);
        try {
            ob.chunks = attribute_record(ob.record, ob.record_offset);
        } catch (const DecodeError& e) {
            fail(e.code(), e.message() + " [object " + std::to_string(ob.id) + " of class " + s.name(ob.cls) + "]", e.offset());
        }
        for (const Chunk& c : ob.chunks) {
            if (c.attr >= s.names.size()) bad("attribute outside the name table", ob.offset);
            if (c.value.kind == Value::Kind::Name && c.value.i >= s.names.size()) bad("name value outside the name table", ob.offset);
        }
        // Relation block: header (grammar in relation_header), then the 7c 0a body, which runs to the
        // end of the block: the body is the one 7c 0a record whose length reaches the block end.
        const ByteView blk = ob.relation_block;
        std::size_t body_at = 0;
        for (std::size_t q = 7; q + 6 <= blk.size; ++q)
            if (blk[q] == 0x7c && blk[q + 1] == 0x0a && le32(blk, q + 2) == blk.size - q) {
                if (body_at) bad("relation block with two candidate bodies", ob.offset);
                body_at = q;
            }
        if (!body_at) bad("relation block without a body", ob.offset);
        ob.osm_version = s.osm_version;
        const ByteView body = tagged(blk, body_at, 0x0a, blk.size);
        if (body_at + body.size != blk.size) bad("bytes after a relation body", ob.offset);
        ob.relation_header = blk.sub(6, body_at - 6);
        ob.relation_body = body.sub(6, body.size - 6);
    }
    s.records_parsed = true;
}

const std::string& Stream::name(std::uint32_t index) const {
    if (index >= names.size()) bad("name index out of range", offset);
    return names[index];
}

const Chunk* Stream::find(const Object& o, std::string_view attribute) const {
    for (const Chunk& c : o.chunks)
        if (names[c.attr] == attribute) return &c;
    return nullptr;
}

const Object* Stream::by_id(std::uint32_t id) const {
    // Ids increase along the stream (gaps allowed): binary search.
    std::size_t lo = 0, hi = objects.size();
    while (lo < hi) {
        const std::size_t mid = (lo + hi) / 2;
        if (objects[mid].id < id) lo = mid + 1;
        else hi = mid;
    }
    return lo < objects.size() && objects[lo].id == id ? &objects[lo] : nullptr;
}

std::vector<Stream> parse_streams(ByteView f) {
    std::vector<Stream> out;
    for (std::size_t at = f.find("\x7c\x04"); at != ByteView::npos; at = f.find("\x7c\x04", at + 1)) {
        if (!object_header_valid(f, at + 6)) continue;
        out.push_back(parse_stream(f, at));
        // the stream's OSMX header ("OSMX" <version compact>) precedes its records
        for (std::size_t m = at; m-- > 0;)
            if (m + 5 <= f.size && std::memcmp(f.data + m, "OSMX", 4) == 0) {
                Cur v{f, m + 4, 0};
                out.back().osm_version = int(v.compact());
                break;
            }
        if (!out.back().osm_version) bad("object stream without an OSMX header", at);
        const Object& last = out.back().objects.back();
        at = last.offset;  // continue after this stream's objects
    }
    for (const Stream& s : out) {
        std::uint32_t last = 0;
        for (const Object& o : s.objects) {
            if (o.flag == 8) continue;
            if (o.id <= last) bad("object ids do not increase", o.offset);
            last = o.id;
        }
    }
    return out;
}

std::vector<Relation> relations(const Stream& s, const Object& o) {
    if (o.preparsed) return o.rels;
    Cur c{o.relation_body, 0, 0};
    std::vector<Relation> out;
    while (c.peek() != 0xfe) {
        Relation r;
        r.attr = pool_name(c);
        const Chunk* chunk = nullptr;
        for (const Chunk& ch : o.chunks)
            if (ch.attr == r.attr) chunk = &ch;
        if (!chunk) bad("relation for an attribute the object does not declare", o.offset);
        if (chunk->array()) {
            counted(c, [&] {
                r.objects.push_back(std::uint32_t(c.compact()));
                const std::uint8_t t = c.u8();
                if (t == 0x80) c.u32();
                else if (t != 0x81) bad("invalid relation tag", o.offset);
            });
        } else {
            r.objects.push_back(std::uint32_t(c.compact()));
        }
        out.push_back(std::move(r));
    }
    ++c.k;
    if (!c.end()) bad("bytes after the relations of an object", o.offset);
    (void)s;
    return out;
}

// <head>, form f = head & 7 (tabulated over all relation blocks of the corpus):
//   f = 1, and f = 5 from version 35: a class name;  other f with bit 2 from version 35: a reference
//   kind byte (00 | 01);  before version 35: head bit 1 (f = 2, 6, 7) adds one leading token that is
//   not a header token, and the tokens may open with an inline class definition
//   fe (fe | e5 <blob> <compact>) <name> fe.
// Then token* [fe], token: e5 <blob> | f6 <guid> | fd <compact> (startup index) | <compact>.
std::vector<HeaderToken> relation_header(const Object& o) {
    if (o.preparsed) return o.header;
    Cur c{o.relation_header, 0, o.relation_offset + 6};
    const std::uint8_t head = c.u8();
    const int form = head & 7;
    const bool v35 = o.osm_version >= 35;
    if (form == 1 || (v35 && form == 5)) {
        pool_name(c);
    } else if (v35 && (head & 4)) {
        const std::uint8_t kind = c.u8();
        if (kind > 1) bad("relation header: invalid reference kind", c.base + c.k - 1);
    }
    std::vector<HeaderToken> out;
    auto token = [&](HeaderToken& t) {
        if (c.peek() == 0xe5) {
            t.is_blob = true;
            t.blob = c.blob();
        } else if (c.peek() == 0xf6) {
            t.is_guid = true;
            t.blob = c.guid();
        } else if (c.peek() == 0xfd) {
            ++c.k;
            t.is_index = true;
            t.number = c.compact();
        } else {
            t.number = c.compact();
        }
    };
    if (!v35 && (head & 2)) {
        HeaderToken lead;
        token(lead);
    }
    if (!v35 && form != 0 && !c.end() && c.peek() == 0xfe && c.k + 1 < c.v.size) {  // class defined here
        ++c.k;
        if (c.peek() == 0xfe) {
            ++c.k;
        } else {
            c.blob();
            c.compact();
        }
        pool_name(c);
        c.expect(0xfe, "class definition end");
    }
    while (!c.end()) {
        if (c.peek() == 0xfe && c.k + 1 == c.v.size) break;  // header end
        HeaderToken t;
        token(t);
        out.push_back(t);
    }
    return out;
}

}  // namespace cadcore::formats::catia
