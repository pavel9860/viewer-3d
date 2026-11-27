// CGM exact B-rep -> flat IR. The IR keeps CGM's parametrizations; the kernel (kernel/brep/src/cgm.rs)
// evaluates them. Only parameter bookkeeping happens here (coedge ranges through edge-curve maps).
//
// Varint-dialect record grammar (every body starts with a reference list: varint count, varint ids):
//   500 body     refs 525 (solids), 520 (sheets), 502 (wires), 524 (points)
//   525 -> 513   volumes; a 513's shells bound a solid
//   520 shell    refs faces; varint n, i16 kind, i16 sense[n], i16
//   512 face     refs surface, loops; i16 sense [u8]
//   521 loop     refs (pcurve, edge)*n [, surface]; varint n, i16 a, i16 b, n x (i16, i16 sense, i16),
//                i16 [, i16 c; c != 0: a 73-byte block]
//   511 edge     refs [edge curve,] vertex vertex point point; u8
//   14 point     refs curves; count, (f64 offset in segment, i32 segment)[count]
//   191-193      edge curve: refs p-curves sharing the first's parameter; count n, n x {i16 map, f64 t0 t1
//                tol, data}, [f64 t0 t1 tol], 4 bytes. Maps: 0 none; 1 identity (t0 t1 tol); 2 affine
//                (b0 b1); 3 (q - off) / s (b0 b1 s off); 5 table (b0 b1, count N, A[N], B[N])
//   200 plane    o U V scale[2] box[4], or o, u8 zero mask, non-zero components, scale[2], box[4]
//   201 cylinder o X Y r u0 u1 v0 v1 su sv; 202 cone o X Y Z a r0 u0 u1 v0 v1 su ...
//   203 sphere   o A B C r ...; 204 torus o X Y Z R r ...; 231 NURBS knots u v, i16 rational, poles, [weights]
//   210 extrusion refs profile; d[3] v0 v1 ...; 211 revolution refs profile; o X Y Z, box, ints, k ...
//   214 offset   refs carrier, base; f64 d ...; 212, 241, 242, 244: carriers of their refs[0]
//   216 fillet, 243 two-curve sections: i32 n, knot vector, P D1 D2 [n][10 | 6], box ...; 251 refs 243
//   150 line     o d t0 t1 (6 f64) or o, unused, d, t0, t1, scale; 151 circle c, i32 i32, r t0 t1 s, offset
//   160 Hermite  knot vector, i32 n, P[2][n], D1[2][n], i16 has2, [D2[2][n]], t0 t1
//   170 B-spline i16 0, knot vector, poles[2], tail; 130 (3D) knot vector, i16 rational, poles, tail
//   155 expr     f64 t0 t1 a b, items u v: i16 code, [1, 4: i16 inputs, i16 0], values, children
//   knot vector  i32 degree, i32, i32 x, i32 n, i32, f64 knots[n], i32 mults[n]
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "formats/catia/cgm.hpp"

namespace cadcore::formats::catia::cgm {

// IR surface kinds: parameters and evaluation.
enum SurfaceKind : std::uint32_t {
    SurfPlane = 1,      // o[3] U[3] V[3]              p = o + u U + v V
    SurfCylinder = 2,   // o[3] X[3] Y[3] r su         p = o + r (cos a X + sin a Y) + v (X x Y), a = u / su
    SurfSphere = 3,     // o[3] A[3] B[3] C[3] k       p = o + cos b (cos a A + sin a B) + sin b C, a = u / k, b = v / k
    SurfTorus = 4,      // o[3] X[3] Y[3] Z[3] R r     p = o + (R + r cos b)(cos a X + sin a Y) + r sin b Z, a = u / R, b = v / r
    SurfNurbs = 5,      // du dv nu nv rational | knots u | knots v | poles[nu][nv][3] | [weights]
    SurfCone = 6,       // o[3] X[3] Y[3] Z[3] a r0 su p = o + (r0 + v cos a)(cos t X + sin t Y) + v sin a Z, t = u / su
    SurfExtrusion = 7,  // pcurve d[3]                 p = C(u) + v d, C: a p-curve on another surface
    SurfOffset = 8,     // base d                      p = base(u, v) + d unit(Su x Sv)
    SurfRevolution = 9,       // o[3] X[3] Y[3] Z[3] k | profile   p = o + Rot_Z(v / k)(C(u) - o), X toward Y
    SurfFillet = 10,          // n | knots[n] | S[n][10] | S'[n][10] | S''[n][10] | u0 u1 v0 v1
                              //   section S = (p1, p2, c, th), quintic Hermite in u;
                              //   p = c + (sin(th - f)(p1 - c) + sin f (p2 - c)) / sin th, f = th (v - v0) / (v1 - v0)
    SurfRuled = 11,           // n | knots[n] | S[n][6] | S'[n][6] | S''[n][6] | u0 u1 v0 v1
                              //   S = (A, B), quintic Hermite in u; p = A + v (B - A)
    SurfExtrusionCurve = 12,  // d[3] | profile                     p = C(u) + v d
    // profile C: 1 pcurve | 2 o[3] d[3] (line) | 3 degree nk knots[nk] np rational poles[np][3] [weights]
};
enum PCurveKind : std::uint32_t {
    PcLine = 1,     // o[2] d[2]           (u, v) = o + t d
    PcCircle = 2,   // c[2] r s            (u, v) = c + r (cos a, sin a), a = t / (r s)
    PcHermite = 3,  // n | knots[n] | P[2][n] | D1[2][n] | has2 | D2[2][n]   quintic (cubic without D2)
    PcBSpline = 4,  // degree nk | knots[nk] | np rational | poles[np][2] | [weights]
    PcExpr = 5,     // t0 t1 a b | node(u) | node(v); node = code nvals vals nkids nodes, leaves read t:
                    //   1: v0 + v1 t + v2 t^2; 2: v3 + v0 cos(v1 - v2 t); 3: k0 / k1; 4: v0 k0 + v1 k1 + v2;
                    //   7: v4 + v0 atan(v1 cos(v2 t + v3))
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

struct DecodeOptions {
    bool pending_kinds = false;  // emit kinds the kernel does not evaluate yet (none at present)
};
// The B-rep of a framed stream (varint dialect, or compact through cgm_compact.hpp). Faces whose
// geometry or grammar is not decoded are left out and listed in `skipped` with the reason.
Ir decode(const Stream& stream, const std::vector<std::vector<std::uint32_t>>& visible = {}, const DecodeOptions& options = {});

// Binary IR (little-endian):
//   "CGMIR\0\0\2", u32 surfaces pcurves vertices edges coedges loops faces shells bodies skipped,
//   surfaces { u32 kind id n, f64[n] }, pcurves { u32 kind surface id n, f64[n] },
//   edges { u32 v0 v1 id }, coedges { u32 pcurve edge, i32 sense, u32 0, f64 t0 t1 },
//   loops { u32 first count }, faces { u32 surface, i32 sense, u32 first_loop loops id },
//   shells { u32 id solid n, n x (u32 face, i32 sense) }, bodies { u32 id flags n, u32 shell[n] },
//   skipped { u32 id type n, char[n] }.
std::vector<std::uint8_t> serialize(const Ir& ir);

}  // namespace cadcore::formats::catia::cgm
