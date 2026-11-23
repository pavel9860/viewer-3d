// CATIA V5 helpers for hosts that resolve assembly links themselves.
#pragma once

#include <optional>
#include <string>

#include "cadcore/bytes.hpp"

namespace cadcore {

// The document id that CATProduct links use to reference this CATPart/CATProduct
// (32 lowercase hex digits), or nullopt when the file carries none. Reads at most the
// first 2 MB.
std::optional<std::string> catia_document_id(ByteView file);

}  // namespace cadcore
