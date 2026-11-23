// CATIA V5 object stream ("OSMX" storage): the serialized feature/product objects of one
// container, their deferred relations and the interned name table.
//
// Layout (all records "7c <tag> <u32le length>", length counted from the 7c):
//   7c 04 ...                       stream header
//   7c 05 len {00|02}               object
//     7c 06 len <u8><u8> <class> <name> ea <u32 id> <refs...>
//     7c 07 len <chunk>* [fe]       attribute record
//     <trailer>
//   de 7c 08 <u32>                  relations: one block per object, in object order
//     7c 09 len <header> 7c 0a len <relation>* fe
//   7c 0b len ... <u8>
//   7c 02 len <compact count> <name>*   interned names; attribute/class ids index it
//
// Every structure is validated; any residue or unknown form throws DecodeError(Malformed).
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::formats::catia {

struct LinkValue {
    std::uint64_t a = 0;
    ByteView blob;              // e5-prefixed blob (document link), empty if none
    std::uint64_t number = 0;   // compact value when there is no blob
    bool has_blob = false;
};

struct Value {
    enum class Kind { Empty, Name, Names, Int, Ints, Double, Doubles, Link, Links, Octets } kind = Kind::Empty;
    std::uint64_t i = 0;              // Name (pool index) / Int
    std::vector<std::uint32_t> names; // Names
    std::vector<std::uint64_t> ints;  // Ints
    std::vector<double> d;            // Double (1 value) / Doubles
    std::vector<LinkValue> links;     // Link (1) / Links
    ByteView octets;                  // Octets
};

struct Chunk {
    std::uint32_t attr = 0;           // pool index of the attribute name
    std::uint64_t head = 0;
    std::uint8_t flags = 0, code = 0; // type descriptor: e8 <flags> <code> <qualifier> | e9 <flags> <code> <extra> <qualifier>
    std::uint8_t extra = 0;
    std::uint32_t qualifier = 0;      // name reference closing the descriptor
    bool array() const { return flags & 0x20; }
    Value value;
};

// Deferred relation values of one object (after parse_records): attribute -> object ids.
struct Relation {
    std::uint32_t attr = 0;
    std::vector<std::uint32_t> objects;  // object ids of this stream
    std::vector<LinkValue> links;        // values that link to another document (SPEC streams)
};
// Relation block header after its head byte and class: link blobs (e5), GUIDs (f6 + 16 bytes)
// and compact numbers in order.
struct HeaderToken {
    bool is_blob = false, is_guid = false;
    bool is_index = false;  // SPEC "fd <n>": a startup named by index (number = n)
    ByteView blob;  // blob or GUID bytes
    bool is_number() const { return !is_blob && !is_guid && !is_index; }
    std::uint64_t number = 0;
};
struct Object {
    std::size_t offset = 0;
    std::uint8_t flag = 0;            // 0, or 2 when data follows the record
    std::uint8_t kind[2] = {0, 0};
    std::uint32_t cls = 0, name = 0;  // pool indices
    std::uint32_t id = 0;
    ByteView refs;                    // header handles (not interpreted)
    ByteView record;                  // 7c 07 body
    std::size_t record_offset = 0, relation_offset = 0;
    ByteView trailer;
    ByteView relation_block;          // 7c 09 block (whole)
    // Filled by parse_records():
    std::vector<Chunk> chunks;
    ByteView relation_header, relation_body;  // header: the tokens after the head byte and class
    int osm_version = 0;  // of the stream (relation header grammar depends on it)
    // SPEC streams (see spec.hpp) carry these already decoded:
    bool preparsed = false;
    std::vector<HeaderToken> header;
    std::vector<Relation> rels;
    LinkValue startup;                // the feature's startup: blob + index, or index alone
    bool startup_by_index = false;    // "fd <index>": a startup named by an index without a blob
};

struct Stream {
    std::size_t offset = 0;
    int spec_version = 0;  // 0: OSM stream; otherwise the version of a SPEC container (spec.hpp)
    int osm_version = 0;   // OSM streams: the version in their OSMX header
    std::shared_ptr<const std::vector<std::uint8_t>> storage;  // owns the bytes of a reassembled stream
    std::vector<Object> objects;
    std::vector<std::string> names;
    bool records_parsed = false;

    const std::string& name(std::uint32_t index) const;
    const Chunk* find(const Object& o, std::string_view attribute) const;
    const Object* by_id(std::uint32_t id) const;
};

// Every object stream of the file with its framing validated: object, header, record and
// relation-block records, the name table, and every class/name reference.
std::vector<Stream> parse_streams(ByteView file);

// Parses the attribute records and relation blocks of every object of `s` (strict).
void parse_records(Stream& s);

std::vector<Relation> relations(const Stream& s, const Object& o);

std::vector<HeaderToken> relation_header(const Object& o);

}  // namespace cadcore::formats::catia
