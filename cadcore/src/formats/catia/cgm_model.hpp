// CGM exact B-rep -> flat exchange IR (varint dialect).
//
// The IR keeps CGM's own parametrizations; the kernel (viewer_3d/kernel, brep/src/cgm.rs) evaluates
// them. Only parameter bookkeeping happens here (coedge ranges through the edge curves' parameter
// maps), no geometry evaluation.
//
// Record grammar (varint dialect; every body starts with a reference list: varint count, varint ids):
//   topology
//     500 body    refs 525 (solid parts), 520 (sheets), 502 (wires), 524 (points)
//     525 -> 513  refs volumes
//     513 volume  refs shells (its shells bound a solid; other shells are sheets)
//     520 shell   refs faces            varint n, i16 kind, i16 sense[n], i16
//     512 face    refs surface, loops   i16 sense [u8]
//     521 loop    refs (pcurve, edge)*n [, surface (V5R8 and later)]
//                                       varint n, i16 a, i16 b, n x (i16, i16 sense, i16), i16
//                                       [, i16 c (with the surface); c = 3: a 73-byte block follows]
//                 the senses chain the loop up to one direction for the whole loop
//     511 edge    refs edge-curve, vertex, vertex, point, point   u8 flags
//     14  point   refs curves           count, count x (f64 parameter, i32)
//   edge curves 191/192/193: refs p-curves (all sharing the parameter of the first)
//     count n, then n x { i16 map, f64 t0 t1 tol on the first p-curve, map data }, [f64 t0 t1 tol]
//     map 0: none; 1: identity (f64 t0 t1 tol); 2: affine between the ranges (f64 b0 b1);
//     3: t_first = s * t + off (f64 b0 b1 s off); 5: piecewise (f64 b0 b1, count N, f64 A[N], B[N]);
//     the closing range is present when some map is 2, 3 or 5; then 4 bytes.
//   surfaces
//     200 plane    o U V scale[2] box[4] (15 f64), or o, u8 zero-mask of U V, the non-zero
//                  components, scale[2], box[4]          p = o + u U + v V
//     201 cylinder o X Y r u0 u1 v0 v1 su sv [ref]       angle = u / (r su), height = v along X x Y
//     202 cone     o X Y Z a r0 u0 u1 v0 v1 su sv ...   (23 f64) angle u / su, slant length v
//     210 extrusion refs profile (p-curve, or edge curve of one); d[3] v0 v1 ... (10 f64)
//     203 sphere   o A B C (axes times r) r ...          angles u / r, v / r
//     204 torus    o X Y Z R r ...                       angles u / R, v / r
//     231 NURBS    knot vectors u, v; i16 rational; poles (v fastest); [weights]
//     212, 241     refs[0] is the carrier surface (241: fillet whose carrier kind is its first i16)
//   p-curves (refs: the surface)
//     150 line     o d t0 t1 (6 f64), or o, unused, d, t0, t1, scale (8 f64)
//     151 circle   c, i32 i32, r, t0, t1, s [, ref]      angle = t / (r s)
//     160 Hermite  knot vector, i32 n, P[2][n], D1[2][n], i16 has-second, [D2[2][n]], t0 t1
//   knot vector: i32 degree, i32, i32, i32 n, i32, f64 knots[n], i32 multiplicities[n]
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "formats/catia/cgm.hpp"

namespace cadcore::formats::catia::cgm {

// Surface kinds and their parameters in the IR.
enum SurfaceKind : std::uint32_t {
    SurfPlane = 1,     // o[3] U[3] V[3]
    SurfCylinder = 2,  // o[3] X[3] Y[3] r su      p = o + r(cos a X + sin a Y) + v (X x Y), a = u / su
    SurfSphere = 3,    // o[3] A[3] B[3] C[3] k    p = o + cos b (cos a A + sin a B) + sin b C, a = u/k, b = v/k
    SurfTorus = 4,     // o[3] X[3] Y[3] Z[3] R r  p = o + (R + r cos b)(cos a X + sin a Y) + r sin b Z, a = u/R, b = v/r
    SurfNurbs = 5,
    SurfCone = 6,
    SurfExtrusion = 7, // pcurve d[3]               p = profile(u) + v d (profile: a p-curve on another surface)      // o[3] X[3] Y[3] Z[3] a r0 su  p = o + (r0 + v cos a)(cos t X + sin t Y) + v sin a Z, t = u / su     // du dv nu nv rational | knots u (nu+du+1) | knots v | poles[nu][nv][3] | weights
};
// P-curve kinds.
enum PCurveKind : std::uint32_t {
    PcLine = 1,     // o[2] d[2]                  (u, v) = o + t d
    PcCircle = 2,   // c[2] r s                   (u, v) = c + r (cos a, sin a), a = t / (r s)
    PcHermite = 3,  // n | knots[n] | P[2][n] | D1[2][n] | has2 | D2[2][n]   C2 quintic (cubic if !has2)
};

struct Ir {
    struct Surface { std::uint32_t kind = 0, id = 0; std::vector<double> p; };
    struct PCurve { std::uint32_t kind = 0, surface = 0, id = 0; std::vector<double> p; };
    struct Edge { std::uint32_t v0 = 0, v1 = 0, id = 0; };
    struct Coedge { std::uint32_t pcurve = 0, edge = 0; std::int32_t sense = 1; double t0 = 0, t1 = 0; };
    struct Loop { std::uint32_t first = 0, count = 0; };
    struct Face { std::uint32_t surface = 0; std::int32_t sense = 1; std::uint32_t first_loop = 0, loops = 0, id = 0; };
    struct Shell { std::uint32_t id = 0; bool solid = false; std::vector<std::uint32_t> faces; std::vector<std::int32_t> senses; };  // solid: bounds a volume (513)
    struct Skip { std::uint32_t id = 0, type = 0; std::string why; };  // a face (or shell) left out

    std::vector<Surface> surfaces;
    std::vector<PCurve> pcurves;
    std::uint32_t vertices = 0;
    std::vector<Edge> edges;
    std::vector<Coedge> coedges;
    std::vector<Loop> loops;
    std::vector<Face> faces;
    std::vector<Shell> shells;
    // a CGM body (500): its shells (indices into `shells`); flags bit 0: named by the stream header's
    // roots; bit 1 + j: listed in list j of the part's CATSeeBodyCont
    struct Body { std::uint32_t id = 0, flags = 0; std::vector<std::uint32_t> shells; };
    std::vector<Body> bodies;
    std::vector<Skip> skipped;
};

// The B-rep of a framed varint-dialect stream. Faces whose geometry or grammar is not decoded are
// left out and listed in `skipped` with the reason. Throws Unsupported for other dialects.
Ir decode(const Stream& stream, const std::vector<std::vector<std::uint32_t>>& visible = {});

// Binary IR (little-endian):
//   "CGMIR\0\0\2", u32 surfaces pcurves vertices edges coedges loops faces shells bodies skipped,
//   surfaces { u32 kind id n, f64[n] }, pcurves { u32 kind surface id n, f64[n] },
//   edges { u32 v0 v1 id }, coedges { u32 pcurve edge, i32 sense, u32 0, f64 t0 t1 },
//   loops { u32 first count }, faces { u32 surface, i32 sense, u32 first_loop loops id },
//   shells { u32 id solid n, n x (u32 face, i32 sense) }, bodies { u32 id flags n, u32 shell[n] },
//   skipped { u32 id type n, char[n] }.
std::vector<std::uint8_t> serialize(const Ir& ir);

}  // namespace cadcore::formats::catia::cgm
