// DXF group-code stream: the (code, value) pairs every DXF file is made of.
//
// ASCII: two lines per group, the code (an integer, may be padded) and the value.
// Binary: the 22-byte sentinel "AutoCAD Binary DXF\r\n\x1a\0", then per group a code
//   (u16le from R13 on; R12 and older use one byte, 255 escaping a following u16le)
//   and a value whose type follows from the code: null-terminated string, f64le, i16le,
//   i32le, i64le, a one-byte bool (290-299), or a binary chunk (310-319, 1004: one length
//   byte then the bytes).
#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::formats::autocad {

enum class ValueType : std::uint8_t { String, Double, Int, Bool, Binary };

// The value type of a group code (DXF reference, "Group Code Value Types").
ValueType value_type(int code);

struct Group {
    int code = 0;
    bool hex = false;      // a binary chunk written as hex text (ASCII files)
    std::string_view str;  // the value text (ASCII, and binary strings); hex text for ASCII binary chunks
    double num = 0;        // Double values; Int/Bool values converted
    long long ival = 0;    // Int/Bool values
    std::size_t offset = 0;
};

// Splits a whole DXF file into groups. Throws DecodeError on a broken stream (a non-numeric
// code line, a truncated binary value, a value that does not parse as its code's type).
// Groups after the first "0 EOF" are ignored; a missing EOF is accepted (`saw_eof` false).
struct GroupStream {
    std::vector<Group> groups;
    bool binary = false;
    bool saw_eof = false;
};
// Large ASCII files are read on up to `threads` threads.
GroupStream read_groups(ByteView data, unsigned threads = 1);

bool is_binary_dxf(ByteView data);
// The file starts with "0 / SECTION" (after optional 999 comments).
bool looks_like_ascii_dxf(ByteView data);

}  // namespace cadcore::formats::autocad
