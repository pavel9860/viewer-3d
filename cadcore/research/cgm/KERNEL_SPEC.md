# Pending IR kinds (cadcore → kernel)

These kinds are decoded by cadcore, but the kernel's IR reader (`kernel/brep/src/cgm.rs`, `surface()`) does
not evaluate them yet. An unknown kind makes the kernel reject the whole IR file. To keep that from
happening, cadcore emits these kinds only on request:

- `cadcli cgm PART OUT.cgmir --pending`
- `cgm::DecodeOptions::pending_kinds`

Without the option, their faces stay in `skipped` with the reason "awaits kernel support".
Once the kernel evaluates them, cadcore's default can flip.

The IR layout is unchanged (CGMIR v2). Each surface is written as `kind, id, n, f64 p[n]`.

## 9: revolution (CGM class 211)

```
p = o[3] X[3] Y[3] Z[3] k | profile
p(t, s) = o + Rot_Z(s / k) (C(t) - o)
Rot_Z turns X toward Y about the axis through o along Z: with w = C(t) - o, x = w.X, y = w.Y, z = w.Z, a = s / k:
p = o + (x cos a - y sin a) X + (x sin a + y cos a) Y + z Z
```

- `t` is the profile's own parameter, which is the **first** coordinate of the face's p-curves on this surface.
  `s` is the second coordinate.
- The profile is `p[13]` followed by its data:
  - `1, pcurve`: an IR p-curve index. C is that p-curve evaluated on its surface (like the extrusion's profile).
  - `2, o[3], d[3]`: a 3D line, C(t) = o + t d.
  - `3, degree, nk, knots[nk], np, rational, poles[np][3], weights[np if rational]`: a 3D B-spline. The
    knots are expanded (already repeated by multiplicity).
- Evidence (research scripts, ground truth = CATIA's stored tessellation):
  - The axis is Z. On CAAGsiStartVolume, 288 stored vertices lie on the surface with axis Z and none with
    axis X or Y.
  - The rotation from the profile equals s / k, in the X→Y sense. On PINZA_DERECHA, the stored vertices
    cover 0.0067..1.7568 rad against s/k = 0..1.7568. Cilindro Cabeca, Tibia and assembly A part-8568 agree.
  - On hair-dryer, the shared edges between pairs of 211 surfaces agree to 2.5e-14.
  - Every partial turn in the corpus starts at a multiple of 2π (or at 0). An origin offset is therefore
    indistinguishable from the data, and it does not affect full turns.

## 10: fillet (CGM class 216, rolling ball)

```
p = n | knots[n] | P[n][10] | D1[n][10] | D2[n][10] | u0 u1 v0 v1
section(u): 10 values (p1[3] p2[3] c[3] th) interpolated by the C2 quintic Hermite of P, D1, D2 over the
            knots, the same basis as PcHermite (h0..h5 in cgm.rs PCurve::uv)
p(u, v) = c + (sin(th - f) (p1 - c) + sin f (p2 - c)) / sin th,   f = th (v - v0) / (v1 - v0)
```

- p1 and p2 are the contact points and c is the centre; |p1 - c| = |p2 - c| = r. th is the arc angle.
- v = v0 is on p1 and v = v1 is on p2. v1 is the radius in every sample seen.
- The face's p-curves are on a 241 carrier, whose parameters are the fillet's (cadcore resolves carriers).
- Evidence:
  - On plug.CATPart, every fillet has between 93 and 109 stored vertices on its surface (tolerance 1e-4 of
    the part's size).
  - Fillet points along each constant-v edge lie on the neighbouring plane or cylinder within 3.8e-6, and
    the swapped orientation (1 - v) fails.
  - Between knots, |p1 - c| and |p2 - c| agree, and the p1-c-p2 angle matches th, to 1.6e-5.
- Not verified: whether f is linear in v. The surface's shape does not depend on that choice, but interior
  trimming curves in v would. No case in the corpus could distinguish it.

## P-curve kind 4: B-spline (CGM class 170, compact code 33)

```
p = degree nk | knots[nk] | np rational | poles[np][2] | weights[np if rational]
(u, v)(t) = sum_i N_i,degree(t) P_i        (rational: sum w_i N_i P_i / sum w_i N_i)
```

- The knots are expanded (already repeated by multiplicity), so nk = np + degree + 1. Every corpus
  instance is non-rational. The flag is there so the layout matches the 3D B-spline profile of kind 9.
- t is the knot parameter. Coedge ranges lie in [knots[degree], knots[nk - degree - 1]], which need
  not start at 0. 59 coedges go past it by less than 1e-3 of the range (vertex tolerance), so
  evaluate slightly past the ends by extending the end spans, not by clamping.
- Revolution (kind 9) profiles of form 1 may now point at a kind-4 p-curve.
- Evidence: on 34 parts, edges whose other coedge decodes meet within 3.3e-15 of the part size (one
  edge excepted, an open line-parameter case in cadcore; FINDINGS.md).
- Emitted only with `--pending`, like kinds 9 and 10.

## 11: ruled two-curve surface (CGM class 243, and 251 through its 243)

```
p = n | knots[n] | P[n][6] | D1[n][6] | D2[n][6] | u0 u1 v0 v1
(A, B)(u): the 6 values (A[3] B[3]) by the C2 quintic Hermite of P, D1, D2 over the knots (as kind 10)
p(u, v) = A(u) + v (B(u) - A(u))
```

- The face p-curves use these (u, v) directly. v often runs past [0, 1] (boxes [-1, 2], [-2, 2],
  [-4.5, 1]), so the surface extends linearly.
- Evidence: stored-mesh vertices lie on these surfaces (6 parts, ground truth = CATIA's stored
  tessellation); p-curve line ends land on stored vertices within 2e-8; coedges meet their
  neighbours within 1.1e-7 (4 parts, all three stored forms).
- Emitted only with `--pending`.

## 12: extrusion of a 3D curve (CGM class 210 whose profile is a 3D curve)

```
p = d[3] | profile          p(u, v) = C(u) + v d
```

- profile: as kind 9's, form 2 (3D line) or form 3 (3D B-spline, knots expanded); u is the curve's
  parameter (the knot parameter for form 3).
- Evidence: on assembly A part-9e0a, 398 of 401 coedges meet their neighbours within 1e-4. The other 3 have
  vertex parameters that come through an edge-curve table map (map 5) and are off; that is open in
  cadcore, not in the surface. Emitted only with `--pending`.

## Still open (not emitted)

