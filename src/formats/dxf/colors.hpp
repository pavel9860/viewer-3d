// AutoCAD Color Index (ACI) palette.
#pragma once

#include "cadcore/scene.hpp"

namespace cadcore::formats::autocad {

// RGB of ACI colour 1..255 (0 and 256 are BYBLOCK / BYLAYER and resolved by the caller).
// 1-9 are fixed; 10-249 are 24 hues x 5 brightness levels x full / half saturation;
// 250-255 are greys.
Color aci_color(int index);

}  // namespace cadcore::formats::autocad
