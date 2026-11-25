#include "formats/registry.hpp"

namespace cadcore::formats {

// Order matters only for equal sniff scores: signature-based formats first.
std::vector<const Decoder*> all_decoders() { return {&catpart(), &catproduct(), &threedxml(), &cgr(), &drc(), &dxf(), &gltf(), &threemf(), &ply(), &off(), &stl(), &obj()}; }

}  // namespace cadcore::formats
