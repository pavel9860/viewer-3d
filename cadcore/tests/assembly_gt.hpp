// Assembly ground truth shared by the CATProduct and 3DXML tests: the occurrence tables of
// STEP exports of the same assemblies (tests/fixtures/*.occurrences.tsv, written by
// tools/step_assembly_gt.py). Every STEP occurrence must exist in the decoded tree with the
// same placement, and the STEP product definitions must map one-to-one onto the decoded
// parent/child nodes.
#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "cadcore/decode.hpp"

namespace gt {

std::vector<std::uint8_t> read_file(const std::filesystem::path& p);

// Options resolving resources (linked documents) inside `dir`.
cadcore::DecodeOptions folder_options(const std::filesystem::path& dir);

// Decodes `file` and checks it against `fixture` (file name in tests/fixtures). Every scene
// warning must satisfy `allowed`. Instance nodes are the children of a reference node that
// themselves have exactly one child (the referenced node); childless children are the
// reference's own representations.
void check_assembly(const std::filesystem::path& file, const std::string& fixture, const std::string& format,
                    const std::function<bool(const std::string&)>& allowed);

}  // namespace gt
