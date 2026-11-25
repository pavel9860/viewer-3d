// Geometry of the DXF entities that draw directly (curves, faces, meshes, hatches), in the
// coordinates of the block or space that holds them.
#pragma once

#include <string>
#include <vector>

#include "cadcore/scene.hpp"
#include "formats/dxf/document.hpp"

namespace cadcore::formats::autocad {

// Collects the drawn primitives of one block (or model space) for one display context.
struct Sink {
    struct Line { Polyline3 pts; Color color; };
    std::vector<Line> lines;
    std::vector<Vec3> tri_pts;   // 3 per triangle
    std::vector<Color> tri_colors;  // 1 per triangle
    std::vector<Vec3> points;
    std::vector<Color> point_colors;

    void line(Polyline3 p, Color c) { if (p.size() >= 2) lines.push_back({std::move(p), c}); }
    void tri(const Vec3& a, const Vec3& b, const Vec3& c, Color col);
    void quad(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, Color col) { tri(a, b, c, col); tri(a, c, d, col); }
    void point(const Vec3& p, Color c) { points.push_back(p); point_colors.push_back(c); }
    bool empty() const { return lines.empty() && tri_pts.empty() && points.empty(); }
};

struct DrawContext {
    Color color;
    bool has_extents = false;      // drawing extents, for infinite lines (RAY, XLINE)
    Vec3 extmin, extmax;
    std::vector<std::string>* notes = nullptr;  // informational remarks (no content lost)
    const std::map<std::string, MlineStyle, std::less<>>* mline_styles = nullptr;
    Color layer_color{255, 255, 255, 255};       // the entity's layer colour (MLINE elements BYLAYER)
};

// Draws `e` when it is one of the directly drawn entity types; returns false for any other type.
// Content of a drawn entity that cannot be shown is appended to `losses` (one short reason each).
bool draw_entity(const Entity& e, const DrawContext& ctx, Sink& out, std::vector<std::string>& losses);

}  // namespace cadcore::formats::autocad
