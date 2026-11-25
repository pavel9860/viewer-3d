// CATIA V5 document identity: the 16-byte UUID that assembly links use to find a document.
//
// It is the GUID of the "CATProdCont" row of the container-name table (u32 fields in the writer's
// byte order: big-endian in files from big-endian UNIX systems). A row is
//   S+0  u32 stream id   S+4 u32 weight   S+8 u32 table id   S+12 u32 index
//   S+16 u16 1, u16 0x6c  S+20 u32 2     S+24 u32 ?         S+28 u32 0
//   S+32 u16 2, u16 0x2081  S+36 u32 f10 = ceil((|n1|+|n2|+2)/8)+1
//   S+40 n1 00 n2 00 (zero padded to a multiple of 8)   then u16 3, u16 0xf7, u32 3, GUID[16]
#pragma once

#include <array>
#include <cstdint>
#include <optional>

#include "cadcore/bytes.hpp"

namespace cadcore::formats::catia {

using Uuid = std::array<std::uint8_t, 16>;

// The document's UUID; nullopt when the file has no (or more than one) CATProdCont row.
std::optional<Uuid> document_uuid(ByteView file);

// The target document of an assembly link blob (11 03 <uuid> ...).
Uuid link_target(ByteView blob, std::size_t offset_for_errors);

}  // namespace cadcore::formats::catia
