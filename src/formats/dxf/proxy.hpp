// AutoCAD proxy graphics: the display list AutoCAD stores with entities other applications may
// not know (groups 92/160 = byte count, 310 = data). It is AutoCAD's own rendering of the entity.
//
//   i32 total_size  i32 command_count, then per command  i32 size (incl. these 8 bytes)  i32 type  data
// Little-endian; points are 3 x f64. Every command is checked against its size. Geometry:
//    2 circle (centre, radius, normal)              4 arc (centre, radius, normal, start vector, sweep, type)
//    3 circle through 3 points                       5 arc through 3 points (+ type)
//    6 polyline (i32 n, n points)                    7 polygon (i32 n, n points; filled when fill is on)
//    9 shell (i32 nv, points, i32 list size, faces: n then n indices; edge/face data follow)
//   12 xline / 13 ray (2 points)                    32 polyline with normal (i32 n, points, normal)
//   44 elliptic arc (centre, normal, major and minor radius, start, end, major-axis angle, type)
// State: 14 colour index, 22 true colour, 20 fill on/off, 29 push matrix (16 f64, row-major), 31 pop.
// 16 layer, 18 linetype, 19 marker, 23 lineweight, 24 linetype scale, 25 thickness, 51 (4 bytes), 34 (8 bytes):
// attributes without geometry. 11 / 38 (text) are reported as text.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "formats/dxf/draw.hpp"

namespace cadcore::formats::autocad {

// The proxy graphics bytes of an entity (its 310 groups after a 92 or 160 byte count), or empty.
std::vector<std::uint8_t> proxy_bytes(const Entity& e);

// Draws the display list into `out` (entity colour `color` until a colour command). Returns false
// when the data does not frame (then nothing is drawn); content it cannot show goes to `losses`.
bool draw_proxy(const std::vector<std::uint8_t>& data, Color color, Sink& out, std::vector<std::string>& losses);

}  // namespace cadcore::formats::autocad
