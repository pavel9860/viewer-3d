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
//                                       [, i16 c (with the surface); c != 0: a 73-byte block follows]
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
//     214 offset   refs carrier, base; f64 d, i16 carrier kind, 4 f64   (the base's parameters)
//     210 extrusion refs profile (p-curve, or edge curve of one); d[3] v0 v1 ... (10 f64)
//     203 sphere   o A B C (axes times r) r ...          angles u / r, v / r
//     204 torus    o X Y Z R r ...                       angles u / R, v / r
//     231 NURBS    knot vectors u, v; i16 rational; poles (v fastest); [weights]
//     212, 241, 242, 244 refs[0] is the carrier, with its parameters (241, 244: fillets whose first
//                  i16 names the carrier kind; 242: the carrier with a parameter box)
//   p-curves (refs: the surface)
//     150 line     o d t0 t1 (6 f64), or o, unused, d, t0, t1, scale (8 f64)
//     151 circle   c, i32 i32, r, t0, t1, s [, ref]      angle = t / (r s)
//     160 Hermite  knot vector, i32 n, P[2][n], D1[2][n], i16 has-second, [D2[2][n]], t0 t1
//     211 revolution refs profile (edge curve, p-curve, 100 line, 130 B-spline curve); o X Y Z (12 f64),
//                  box (4 f64), i32, i32, k, f64, f64, f64, i32, f64
//     216 fillet   i32 n, knot vector (degree 5, n knots), f64 P[n][10], D1[n][10], D2[n][10], u0 u1 v0 v1,
//                  more (not decoded)
//   3D curves: 100 line o[3] d[3] ...; 130 B-spline knot vector, i16 rational, poles[3], [weights], ...
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
    SurfOffset = 8,    // base d                     p = base(u, v) + d unit(Su x Sv)
    SurfExtrusion = 7, // pcurve d[3]               p = profile(u) + v d (profile: a p-curve on another surface)
    // Emitted only with DecodeOptions::pending_kinds (the kernel does not evaluate them yet;
    // research/cgm/KERNEL_SPEC.md):
    SurfRevolution = 9, // o[3] X[3] Y[3] Z[3] k | profile      p(t, s) = o + Rot_Z(s / k) (C(t) - o),
                        //   Rot_Z turning X toward Y; t: the profile's parameter (first p-curve coordinate),
                        //   s: second coordinate; profile C: 1 pcurve-index (a p-curve on its surface),
                        //   2 o[3] d[3] (3D line C = o + t d), 3 B-spline (degree, nk, knots[nk] expanded,
                        //   np, rational, poles[np][3], weights[np] when rational)
    SurfRuled = 11,     // n | knots[n] | P[n][6] | D1[n][6] | D2[n][6] | u0 u1 v0 v1
                        //   (A, B)(u) = C2 quintic Hermite of the 6 values (A[3] B[3]), as SurfFillet;
                        //   p(u, v) = A(u) + v (B(u) - A(u)) (CGM 243, and 251 through its 243)
    SurfExtrusionCurve = 12, // d[3] | profile       p = C(u) + v d (CGM 210 whose profile is a 3D curve;
                             //   profile as SurfRevolution's, forms 2 and 3)
    SurfFillet = 10,    // n | knots[n] | P[n][10] | D1[n][10] | D2[n][10] | u0 u1 v0 v1
                        //   section(u) = C2 quintic Hermite of the 10 values (p1[3] p2[3] c[3] angle), as
                        //   PcHermite; p(u, v) = c + (sin(th - f) (p1 - c) + sin f (p2 - c)) / sin th,
                        //   th = section angle, f = th (v - v0) / (v1 - v0)      // o[3] X[3] Y[3] Z[3] a r0 su  p = o + (r0 + v cos a)(cos t X + sin t Y) + v sin a Z, t = u / su     // du dv nu nv rational | knots u (nu+du+1) | knots v | poles[nu][nv][3] | weights
};
// P-curve kinds.
enum PCurveKind : std::uint32_t {
    PcLine = 1,     // o[2] d[2]                  (u, v) = o + t d
    PcCircle = 2,   // c[2] r s                   (u, v) = c + r (cos a, sin a), a = t / (r s)
    PcHermite = 3,  // n | knots[n] | P[2][n] | D1[2][n] | has2 | D2[2][n]   C2 quintic (cubic if !has2)
    PcExpr = 5,     // t0 t1 a b | tree(u) | tree(v); a node: code nvals vals[nvals] nkids kids...;
                    //   leaves take the curve parameter t (cgm_model.cpp expr()); CGM class 155, emitted
                    //   only with DecodeOptions::pending_kinds
    PcBSpline = 4,  // degree nk | knots[nk] (expanded) | np rational | poles[np][2] | weights[np if rational]
                    //   (u, v) = sum N_i,degree(t) P_i; t is in the knot range (CGM class 170; emitted
                    //   only with DecodeOptions::pending_kinds until the kernel evaluates it)
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
struct DecodeOptions {
    bool pending_kinds = false;  // emit surface kinds the kernel does not evaluate yet (SurfRevolution, SurfFillet)
};
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
