// Canonical test parts with analytic expectations, and writers for every supported format.
//
// The writers are independent of the decoders (they share only the Scene types),
// so encode -> decode -> compare against the analytic values tests each decoder
// against ground truth.
#pragma once

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "cadcore/scene.hpp"

namespace cadcore::testgen {

enum class PrimKind { Triangles, Strip, Fan };

struct Prim {
    PrimKind kind;
    std::vector<std::uint32_t> refs;
};

struct Face {
    std::vector<Prim> prims;
    Color color;
    bool has_normal = false;
    std::array<float, 3> normal{0, 0, 0};
};

struct Body {
    std::string name;
    std::vector<std::array<float, 3>> points;
    std::vector<Face> faces;
    bool solid = true;
};

struct Instance {
    std::uint32_t body;
    Matrix4 transform = kIdentity;
};

struct Expected {
    std::size_t triangles;  // world triangles (instances expanded)
    double area;
    double volume;          // NaN when the model is not closed
    std::array<double, 3> min, max;
};

struct Model {
    std::string name;
    std::vector<Body> bodies;
    std::vector<Instance> instances;  // empty = each body once, untransformed
    Expected expected;
};

// Triangles of a face under the shared primitive rules (strips alternate, fans pivot on the first ref).
std::vector<std::array<std::uint32_t, 3>> triangulate(const Face& f);

std::vector<Model> canonical_models();

// Scalable models for complexity tests.
// grid_model: (k+1)^2 points, 2k^2 triangles; one strip face per row, or one face per quad.
Model grid_model(int k, bool face_per_quad);
// instances_model: one cube placed `count` times (node/instance handling).
Model instances_model(int count);

// A written test file: name -> bytes (several entries for multi-file formats; the first is the main file).
using Files = std::vector<std::pair<std::string, std::vector<std::uint8_t>>>;

struct CatPartOptions {
    int encoding = 0;        // 0 = prefix table, 1 = raw, 2 = texcoords + table
    bool force_u16 = false;  // write u16be references even when u8 would fit
    bool lod = false;        // add one LOD mesh before the face list
    bool u16_counts = false; // loop counts as u16be (older CATIA releases, e.g. V5R17); implies u16 references
};

Files write_stl_binary(const Model&);
Files write_stl_ascii(const Model&);
Files write_obj(const Model&);
Files write_ply(const Model&, int encoding /*0 ascii, 1 le, 2 be*/);
Files write_off(const Model&, bool binary);
Files write_gltf(const Model&, bool glb);
Files write_3mf(const Model&);
Files write_catpart(const Model&, const CatPartOptions&);
Files write_dxf(const Model&, int format /*0 ascii, 1 binary R13+, 2 binary R12*/);

// Every (format variant, writer) pair; `catpart` variants only for models without instances.
struct Variant {
    std::string id;  // e.g. "ply-be"
    Files (*write)(const Model&);
};
const std::vector<Variant>& variants();

}  // namespace cadcore::testgen
