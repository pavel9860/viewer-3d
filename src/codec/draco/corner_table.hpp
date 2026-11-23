// Corner tables for Draco connectivity: the position table built by EdgeBreaker and the
// per-attribute tables that cut it along attribute seams.
//
// Corner c belongs to face c / 3; Next/Previous rotate inside the face. Opposite(c) is the
// corner facing c across the edge opposite to c, or -1 on a boundary.
#pragma once

#include <cstdint>
#include <vector>

#include "codec/draco/bitstream.hpp"

namespace cadcore::codec::draco {

using IVec = std::vector<std::int32_t>;
constexpr std::int32_t kInvalid = -1;

inline std::int32_t next_corner(std::int32_t c) { return c < 0 ? c : (c % 3 == 2 ? c - 2 : c + 1); }
inline std::int32_t prev_corner(std::int32_t c) { return c < 0 ? c : (c % 3 == 0 ? c + 2 : c - 1); }

class CornerTable {
public:
    IVec opposite;   // per corner
    IVec vertex;     // corner -> vertex
    IVec left_most;  // vertex -> its left-most corner (-1: isolated)

    void reset(std::size_t num_faces) {
        opposite.assign(num_faces * 3, kInvalid);
        vertex.assign(num_faces * 3, kInvalid);
        left_most.clear();
    }
    std::int32_t num_faces() const { return std::int32_t(vertex.size() / 3); }
    std::int32_t num_corners() const { return std::int32_t(vertex.size()); }
    std::int32_t num_vertices() const { return std::int32_t(left_most.size()); }
    std::int32_t add_vertex() { left_most.push_back(kInvalid); return std::int32_t(left_most.size()) - 1; }

    static std::int32_t Next(std::int32_t c) { return next_corner(c); }
    static std::int32_t Previous(std::int32_t c) { return prev_corner(c); }
    std::int32_t Opposite(std::int32_t c) const { return c < 0 ? kInvalid : opposite[std::size_t(c)]; }
    std::int32_t Vertex(std::int32_t c) const { return c < 0 ? kInvalid : vertex[std::size_t(c)]; }
    std::int32_t LeftMostCorner(std::int32_t v) const { return v < 0 ? kInvalid : left_most[std::size_t(v)]; }
    std::int32_t SwingLeft(std::int32_t c) const { return Next(Opposite(Next(c))); }
    std::int32_t SwingRight(std::int32_t c) const { return Previous(Opposite(Previous(c))); }
    std::int32_t GetLeftCorner(std::int32_t c) const { return c < 0 ? kInvalid : Opposite(Previous(c)); }
    std::int32_t GetRightCorner(std::int32_t c) const { return c < 0 ? kInvalid : Opposite(Next(c)); }
    bool IsOnBoundary(std::int32_t v) const { return SwingLeft(LeftMostCorner(v)) == kInvalid; }
};

// Attribute connectivity: the position table with seam edges removed and vertices split
// where the attribute is discontinuous.
class AttributeCornerTable {
public:
    const CornerTable* pos = nullptr;
    std::vector<bool> edge_on_seam;    // per corner: the edge opposite the corner is a seam
    std::vector<bool> vertex_on_seam;  // per position vertex
    IVec vertex;                       // corner -> attribute vertex
    IVec left_most;                    // attribute vertex -> left-most corner

    void init(const CornerTable& p) {
        pos = &p;
        edge_on_seam.assign(std::size_t(p.num_corners()), false);
        vertex_on_seam.assign(std::size_t(p.num_vertices()), false);
        vertex.assign(std::size_t(p.num_corners()), kInvalid);
        left_most.clear();
    }
    void add_seam_edge(std::int32_t c);
    bool recompute_vertices();  // false on a corrupt (cyclic) seam configuration

    std::int32_t num_faces() const { return pos->num_faces(); }
    std::int32_t num_vertices() const { return std::int32_t(left_most.size()); }
    static std::int32_t Next(std::int32_t c) { return next_corner(c); }
    static std::int32_t Previous(std::int32_t c) { return prev_corner(c); }
    bool IsCornerOppositeToSeamEdge(std::int32_t c) const { return edge_on_seam[std::size_t(c)]; }
    std::int32_t Opposite(std::int32_t c) const { return c < 0 || edge_on_seam[std::size_t(c)] ? kInvalid : pos->Opposite(c); }
    std::int32_t Vertex(std::int32_t c) const { return c < 0 ? kInvalid : vertex[std::size_t(c)]; }
    std::int32_t LeftMostCorner(std::int32_t v) const { return v < 0 ? kInvalid : left_most[std::size_t(v)]; }
    std::int32_t SwingLeft(std::int32_t c) const { return Next(Opposite(Next(c))); }
    std::int32_t SwingRight(std::int32_t c) const { return Previous(Opposite(Previous(c))); }
    std::int32_t GetLeftCorner(std::int32_t c) const { return Opposite(Previous(c)); }
    std::int32_t GetRightCorner(std::int32_t c) const { return Opposite(Next(c)); }
    bool IsOnBoundary(std::int32_t v) const {
        const std::int32_t c = LeftMostCorner(v);
        return c == kInvalid || SwingLeft(c) == kInvalid;
    }
    bool IsCornerOnSeam(std::int32_t c) const { return vertex_on_seam[std::size_t(pos->Vertex(c))]; }
};

// Visits the corners around the vertex of `start`: swings left until a boundary or a full
// turn, then continues right from `start` ("VertexCornersIterator").
template <class Table, class F>
void for_each_vertex_corner(const Table& t, std::int32_t start, F&& f) {
    std::int32_t c = start;
    bool left = true;
    std::size_t guard = 0;
    const std::size_t limit = std::size_t(t.num_faces()) * 3 + 1;
    while (c != kInvalid) {
        if (++guard > limit) bad("corner cycle around a vertex");
        f(c);
        if (left) {
            c = t.SwingLeft(c);
            if (c == kInvalid) {
                c = t.SwingRight(start);
                left = false;
            } else if (c == start) {
                c = kInvalid;
            }
        } else {
            c = t.SwingRight(c);
        }
    }
}

}  // namespace cadcore::codec::draco
