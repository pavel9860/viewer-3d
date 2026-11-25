// Clipping drawn geometry to a closed 2D boundary (xy; even-odd inside), for XCLIP.
#pragma once

#include <vector>

#include "formats/dxf/draw.hpp"

namespace cadcore::formats::autocad {

// Sink content of `in` that lies inside `boundary` (a closed loop in xy): lines split at the
// boundary, triangles intersected with the region (z from each triangle's plane), points kept
// when inside.
Sink clip_to_boundary(const Sink& in, const Polyline3& boundary);

}  // namespace cadcore::formats::autocad
