// The format-independent result of decoding: meshes, edges, an instance tree.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cadcore {

struct Color {
    std::uint8_t r = 200, g = 200, b = 200, a = 255;
    bool operator==(const Color&) const = default;
};

enum class BodyKind : std::uint8_t { Unknown, Solid, Surface };

// Triangles: `indices` holds 3 entries per triangle. Points: a point cloud, `indices` is empty
// and every vertex is one point.
enum class Topology : std::uint8_t { Triangles, Points };

// Triangle range [first, first + count) that came from one source face / group.
struct FaceGroup {
    std::uint32_t first = 0;
    std::uint32_t count = 0;
    Color color;
};

struct Mesh {
    std::string name;
    std::vector<float> positions;       // x, y, z per vertex
    std::vector<float> normals;         // optional, parallel to positions
    std::vector<std::uint8_t> colors;   // optional per-vertex RGBA
    std::vector<std::uint32_t> indices; // 3 per triangle
    std::vector<FaceGroup> groups;      // optional per-face colours / picking ranges
    Color color;                        // default colour when no groups
    BodyKind kind = BodyKind::Unknown;
    Topology topology = Topology::Triangles;

    std::size_t vertex_count() const { return positions.size() / 3; }
    std::size_t triangle_count() const { return indices.size() / 3; }
};

struct Polyline {
    std::vector<float> points;  // x, y, z per point
    Color color{0, 0, 0, 255};
};

using Matrix4 = std::array<double, 16>;  // column-major, like glTF
constexpr Matrix4 kIdentity{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

struct Node {
    std::string name;
    Matrix4 transform = kIdentity;          // relative to the parent
    std::vector<std::uint32_t> meshes;      // indices into Scene::meshes
    std::vector<std::uint32_t> polylines;   // indices into Scene::polylines
    std::vector<std::uint32_t> children;    // indices into Scene::nodes
    // Named non-geometric data of the node as decoded (e.g. a CATIA product's type and
    // application attributes such as a piping specification), in file order.
    std::vector<std::pair<std::string, std::string>> properties;
};

struct Scene {
    std::string format;                 // decoder name, e.g. "stl-binary"
    std::vector<Mesh> meshes;
    std::vector<Polyline> polylines;
    std::vector<Node> nodes;            // nodes[0] is the root
    double unit_to_mm = 1.0;            // multiply coordinates by this to get millimetres
    std::vector<std::string> warnings;  // content that could not be shown (a loss): always user-visible
    std::vector<std::string> notes;     // informational, no geometry lost (e.g. unsupported material extension)

    std::size_t triangle_count() const;
    std::size_t vertex_count() const;
    // Root node referencing every mesh/polyline, for flat formats.
    void make_flat_root(const std::string& name = "root");
};

}  // namespace cadcore
