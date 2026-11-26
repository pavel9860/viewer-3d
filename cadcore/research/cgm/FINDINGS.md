# CGMGeom decode — findings so far (paused 2026-11-23)
Prototypes copied (uncommitted) to cadcore/research/cgm/spike/.
Framing (verified corpus-wide, 1,534 streams):
- varint dialect V5R8–R19 (480): header `01|u32 len, u32 strlen+ver, root-id varints, 17 xx`; dir `e5 0d 06 01 00 u32 level u32 size`;
  record `e5 0d|0f fl u16 type u32 len u32 id [i16 if 0f] body[len]`; chains exactly to dir size; then trailer (tolerances, levels). Big-endian variant: lead byte 00 (CXR6 Unix, 2 files).
- compact dialect 2009Fall+ (1,022): dir `a8 06 7a u32 level u32 size`; record `m fl type8 len id`, m high nibble b→len 1B, a→len 4B; low nibble = lenw+idw. Type byte is per-file (not global); body int/double compression differs 2009Fall (old varints) / 2011 / 2012+. Not decoded yet.
- fixed dialect CXR4/5, CGA-V2R4 (30): `u32 id, u16 type, b3 00, u32 0, u32 len, body`, fixed-width fields.
Varint ints: byte>=0x80 → value&0x7f; else bits3–6 = presence mask of bytes 0–3 (`30 84 02`=0x028400). R8–R10/CXR: some counts are raw i32.
Every record body starts with a ref list (count + refs); all refs resolve.
Types (varint dialect): 500 body→525→513→520 shell→512 face(surface, loops 521)→521 loop = (pcurve, edge511)* + surface + i16 senses;
511 edge = [edgecurve 191/192/193, v0, v1, pt14, pt14]+flag; 510 vertex→13→14 (point on curves: n×(f64 param,i32)), 10 3D point, 12 point on surface.
Surfaces: 200 plane (origin,U,V,2 scales,uv box; 6-bit zero-mask form when len%8==1), 201 cylinder (o,X,Y,r,u0,u1,v0,v1,s1,s2[,x]; u=arc length, θ=u/r), 204 sphere?, 202/203 cone/torus?, 212/214/216/231/241/242/244 TBD.
P-curves: 150 line (o2,d2,t0,t1; old form 8 doubles), 151 circle (c2, 2×i32, r, t0, t1, scale, offset; θ=(t-off)/(r·s)), 160 nurbs TBD.
Edge curve 192 payload: varint n, i16 flags, [t0,t1,tol], per extra pcurve: flag1 identity [b0,b1,tol], flag2 linear [b0,b1,s,off] (t0=s·t1+off), flag4 table [b0,b1, N, N params, N params]; then [c0,c1,tol]+4 bytes — parser still fails on many (ecfail).
B-rep-only parts: 319 (cadcli): 227 varint, 56 compact, 30 fixed, 6 no CGM stream.
Next: finish 192/193 grammar, remaining surfaces/pcurves, then C++ module + cadcli + kernel import.

## Status 2026-11-23 (second pass)
Implemented: cadcore `src/formats/catia/cgm.{hpp,cpp}` (location, framing varint+compact),
`cgm_model.{hpp,cpp}` (varint-dialect decode to the flat IR; grammar in its header), `cadcli cgm`;
kernel `brep/src/cgm.rs`, `cgm_mesh.rs`, bin `cgm`. Check: `research/cgm/compare_tess.py`.
Facts established (corpus-verified):
- Edge-curve maps (191/192/193): all 71,407 parse exactly; maps nest.
- Cylinder/sphere/torus/circle parameters are arc lengths; the trailing "offset" double is only a
  reference parameter (0 mismatches over all old-dialect vertices).
- 160 = C2 quintic Hermite p-curve; 231 = NURBS surface; 212/241 resolve to their carrier.
- Loop: senses at i16 index 3+3j; the stored senses chain up to one global sign per loop.
- Orientation: outward normal = -(shell sense) * (Su x Sv); the face record's i16 does not enter.
Results vs CATIA's stored tessellation (old dialect, parts with all faces decoded): 150/192 closed
and within 1% volume, 16 more within 5%.
Open: which body CATIA shows (needs the feature tree -> CGM body link; header roots only partly
identify it); surfaces 202, 210, 211, 214, 216 (fillet), 242, 244; p-curves 152, 155, 170;
compact dialect bodies (56 B-rep-only parts); fixed dialect framing (30); sheet-only parts.

## Status 2026-11-25
- Shown bodies: the part's CATSeeBodyCont (u32-counted id lists, 3 or 5 lists in all 1,500 parts
  that have one) lists the CGM bodies CATIA displays: lists 0 and 1 (5-list form: origin planes +
  results in list 0; 3-list form: results in list 1, planes in list 2). IR v2 has a body table
  (bodies share unchanged shells across versions).
- New surfaces: cone 202 (angle u/su, slant v), extrusion 210 (profile p-curve + direction),
  offset 214 (refs[1] base + d along unit(Su x Sv); all 614 checked points), 242/244 carriers.
- Compact dialect: type codes are one fixed enumeration (cgm.cpp compact_class). 2009Fall bodies
  keep varint reference lists but pack scalars (ints as (v << 2) | 1 bytes, masked doubles);
  2011+ pack further. Not decoded.
- vs stored tessellation (old dialect): 171 of 192 fully decoded parts closed and within 1%.
- B-rep-only: 134 closed solids, 21 open (216 fillet, 155 composite p-curve, 211, 251, 152, 170),
  71 without a solid (surface parts; sheets not shown), 90 not meshed (compact dialect 55, fixed
  dialect 28, no CGM stream located 5, damaged store 2).

## Compact dialect bodies (2011+), in progress
Spike: `spike/cpt.py` (records, ids), `cgraph.py` (reference graph), `pairs.py` / `pairs2.py`
(old <-> compact record pairs of the same part saved by two releases).
- Signed integer (ids, references, small fields): first byte b, low 2 bits:
  01 -> +(b >> 2); 11 -> -((b >> 2) + 1); 00 -> (b >> 2) bytes follow, little-endian, positive;
  10 -> the same, negative (-v - 1).
- Record ids: the id field holds the delta from the previous record's id (1 = `05`); monotone in
  940 of 957 streams (17 to check); 91 id fields carry extra bytes after the delta.
- References: the signed difference self - target (`03` = next record, `05` = previous).
- Variable-length reference lists keep a varint count (0x80 | n); fixed-arity records omit it.
- Two families. Levels 79/80 (2011Spring): ids equal the old releases' ids, so 14 parts pair
  record by record with their V5R15-R17 saves (`pairs.pkl`). Levels 74-78, 81 (2012-2017, all of the
  user's own parts): ids renumbered and records reordered, pairing needs content matching.
- 2011Spring findings: edge (94) = v0, v1, p0, p1 references + flag (its edge curve is dropped;
  the edge curve (36) gains a reference to its edge after its p-curve list); vertex (93) = one
  reference + flag; point (6) = count + references, then per curve f64 parameter and a small int
  (the old i32); line p-curve (24) = surface reference, form int (0 general: 6 f64; 1 vertical:
  x t0 t1; 2 horizontal: y t0 t1); edge-curve maps: flag as small int, zero tolerances and the
  4 closing bytes dropped.
- 2012+: line p-curves use other form codes (1, 2, 8, 13, 72, 136, ...) with 0-5 f64; point
  parameters that equal a curve end become 1-byte codes (0x84 / 0x88 seen).
- 2012+ pairing by content (`spike/pairs3.py`: a compact record is the old record of the same
  class whose 8-byte doubles it contains): of the 7 parts saved both ways, 6 have both files here,
  giving 4,174 line p-curves, 1,215 planes, 615 cylinders, 386 Hermite p-curves (points and
  topology carry no doubles; they must pair through references next).
- 2012+ line p-curve (24) after its surface reference: a flag byte. Bit 0: d = (0, 1), stores
  o.x t0 t1 (o.y = 0); bit 1: d = (1, 0), stores o.y t0 t1 (o.x = 0); bit 6: d = (0, -1);
  bit 7: d = (-1, 0); bit 3: symmetric range, only t0 stored (t1 = -t0); with no direction bit
  the direction is stored. Bits 2 and 4 make values derived, not stored (form 6 keeps only the
  range length pi/4; form 17 stores 1.646, 2.646 for o.y 1.646, range [0, 0.785]); not solved.

## Status 2026-11-25 (compact bodies, 2011Spring family decoded)
Implemented: `cgm_compact.{hpp,cpp}` (compact record -> canonical varint-dialect record; grammar in
its header), `compact_layout`, compact ids delta-decoded in `frame`, `Stream.features`; the decoder
runs unchanged on canonical records. Reference transcoder: `spike/c11.py`.
- The stream header's bytes after the second version string are per-file features, not the
  directory level: `bf 77`, `b5 77`, `b5 75` (2011Spring), `06 b5 75 ..` (2011Fall), `08 b5 75 ..`
  (2012-2022, GUIDs follow). Byte 0 bit 3: small integers as compact signed integers (else raw
  i16 / i32 of the varint dialect); byte 1 bit 1: edges drop their edge-curve reference.
- Only `bf 77` parts have old twins (14, id-paired): every paired record of the decoded classes
  transcodes byte-exact (topology, points, edge curves, lines, circles, Hermite, planes,
  cylinders, cones, spheres, tori). `b5 75` / `b5 77`: grammar-exact over the corpus.
- Class codes (from reference arity vs the old classes): 2 = 10, 4 = 12, 14 = 100, 26 = 152,
  29 = 155, 35 = 191 (not 193), 37 = 193, 41 = 202, 44 = 210, 45 = 211, 46 = 212, 48 = 214,
  50 = 216, 52 = 231, 55 = 241, 56 = 242, 57 = 243, 59 = 244, 64 = 251.
- Knot vectors: deg, n, form byte (4: field x = 0; 0: x follows), knots, 1 = standard
  multiplicities (ends deg + 1, interior deg - 2: holds for every NURBS of degree 3 and 5 and the
  quintic Hermite curves). NURBS: then rational, poles, weights, a 152/154-byte block.
- Objects of one class can use two layouts in one file (circles with / without the two i32;
  edge-curve map 1 with / without its identity range): the record length decides, exactly one
  layout fits every instance. Edge curves may list their curves twice (A B A B, header 08, 264
  records): the second copy shares the maps. Loop blocks follow any non-zero closing i16 (1 in
  compact streams, 3 in the varint dialect).
- Results (292 2011Spring compact parts): 82,509 of 86,340 faces decode (95.6%); the rest are
  classes the varint-dialect decoder lacks too (216, 211, 243 surfaces; 152, 155, 170 p-curves;
  codes 33, 28, 22, 16). vs stored tessellation: 209 complete parts, 120 closed within 1%, 54 more
  within 5%. Old dialect unchanged (172 closed within 1% of 218 complete).
- Compact integer, corrected: first byte b, low bits 00/10: n = bits 2-3 bytes follow and b's high
  nibble tops the value (v = (b >> 4) << 8n | bytes). The old reading (n = b >> 2) only held for
  b < 0x10; the 2-byte plane headers (`b4 d4`, `e4 59`, `d4 ed`, ...) are this form.
- 2012+ (`08 b5 75`) decoded (cgm_compact.hpp has the grammar): topology as 2011 with compact
  integers and the edge-curve reference kept; values that equal parameters of other records are
  omitted, resolved in passes (surface boxes, p-curve ranges, edge-curve ranges, points).
  Verified by closure on the first streams: every p-curve pair of an edge curve meets in 3D (2,058
  checks, 0 apart) and every vertex parameter lands on all its curves (1,386 checks, 0 apart).
  - plane: compact-integer header v (bits above) then the stored doubles.
  - line: flag byte (directions in bits 0, 1, 6, 7; +directions store [c, t0, t1], bit 4: c = box
    start of the other parameter, bit 2: t0 = box start, bit 3: t1 = box end or -t0; -directions and
    general lines store the range first, then o (and d)).
  - edge curve: byte 0x80 | maps << 2 | bit 0 / bit 1: the closing range's start / end equal the first range's;
    for one curve with map 0 or 1 the first range is the p-curve's own range.
  - point: code byte; bit 7: refs 0 and 1 share the parameter; bits 2-3: refs[0] at the start (2)
    or end (1) of its edge curve's range, or a stored f64 (0). The other refs' fields are skipped
    (their parameters follow from the edge-curve maps); the trailing int's presence depends on the
    curve classes (not needed by the decoder).
  - knot vectors: form-byte bit 3 omits the multiplicity byte (standard multiplicities).
  - surface boxes: cylinder u0 u1 v0 v1; cone the same after its scalars; sphere and torus store
    angles (scaled by r / R to the arc-length parameters); NURBS knot ends; carriers their first
    reference's, offsets their base's; extrusions u = the profile's range, v from the payload.
- Edge-curve byte x = 0x80 | map count << 2 | closing bits (no separate count in 2012+).
- Open (2012+): classes 152,
  155, 251 and code 16 (as in the varint dialect), the meaning of line bit 5 on +directions.
- 2012+ fixes found by meshing (an edge's coedges / loops in UV, not covered by the closure checks):
  map-1 edge curves with two curves take b0 / b1 from a0 / a1 by the x bits (the single-curve rule
  applies only to one-curve edge curves); point start / end codes use the edge curve's range, or the
  first curve's range when that is the infinite default; points keep every ref whose code is start
  / end; general and -direction lines: bit 2 = range derived (stored: origin [and direction]), bits
  4 / 5 = origin on the box start / end of the other parameter.
- Results 2012+ (657 streams, vs stored tessellation): 514 parts fully decoded; 187 closed within
  1%, 128 more within 5%, 199 off. Open lead: a full circle split into two edges on one edge curve,
  where an "end" code resolves to the full range end for the first half. 2011Spring and the varint
  dialect unchanged (120 / 172 within 1%).
- B-rep-only compact parts (55): the 14 of the 2011Spring family decode completely but hold no solid
  (9 only the origin planes, the rest sheets); the 30 of 2012+ decode (meshes as above); 11 are 2009Fall (bodies not decoded).
- Points (14), all dialects: per curve (f64 offset, i32 segment number). The parameter is the start of
  that segment plus the offset; on a Hermite p-curve segment s starts at knot[s + x - 1], x being the
  knot vector's third int (0, or negative on curves trimmed from a longer one: -7 seen), and an edge
  curve's parameter is its first curve's. Segment 1 on every other curve. Ignoring the segment cut
  Hermite edges short (40 of 148 stored vertices off the surfaces on one part, hex-nut chamfers on
  others). Corpus vs stored tessellation: within 1% 484 -> 501, off 240 -> 226, open edges -1348;
  31 partial parts gain open edges while their stored vertices now lie on the surfaces (mesher).
- Volume vs stored tessellation is biased by CATIA's coarse mesh (a 10-gon around a cylinder: +6.9%;
  thin walls: up to 10%). Stored vertices lying on the decoded surfaces is the sharper check: of 49
  meshing 2011Spring parts marked off, 40 had every stored vertex on the surfaces.

## 2011Fall compact layout (feature bytes 06 b5 75)
- The 2011Spring body grammar, with compact integers and edges keeping their edge-curve reference;
  no derived values (those start with 08, 2012+). Of the 8 layout combinations, only this one decodes
  any face: all faces of the 4 distinct parts (8 files) decode, with no skipped records.
- vs stored tessellation: all 4 are closed (no open edges), and every stored vertex is within 1.4e-4 of
  our mesh (well under 1e-3 of the diagonal). The volume differences (0.02% to 3.4%) come from
  CATIA's coarse stored mesh.

## 2009Fall compact layout (directory levels 83-85)
- Bodies use the varint dialect's grammar: references are varint counts and absolute varint ids.
  Header feature byte 0 is `0b` (bit 3: scalar integers as compact signed integers) or `01` (raw
  i16 / i32, as in the varint dialect). This is the same bit as in 2011Spring. 63 of 65 parts are
  `0b`; 2 are `01`. The ±1 sequences of shells and loops are raw i16 in both forms.
- Class differences from the varint dialect:
  - Lines take an int form (0 general, 1 / 2 axis-parallel, as in 2011Spring).
  - Hermite curves: after the value count, compact streams repeat the knot vector's fifth field. That
    field is b, the number of knots before parameter 0. After t0 t1 comes an int −2.
  - Edge curves end with an i32. Edge curves that list their curves twice (A B A B) add two ints
    (not decoded).
  - Loops always carry the closing int c; c ≠ 0 adds the 73-byte block (4 f64, 4 ints, u8, 24 bytes).
- Every record of the corpus's 65 2009Fall parts (219,050) is consumed exactly by this grammar.
  Records of 4 parts saved by V5R10 / V5R17 show the same structure. The twins were edited between
  saves, so their doubles differ in low bits and pairing is not byte-exact.
- Results: 9,955 of 11,107 faces decode. The rest are classes not decoded in any dialect (211, 216,
  251, code 33). vs stored tessellation: 41 complete parts (25 within 1%, 9 within 5%), 13 partial,
  11 without a stored mesh. On 39 parts every stored vertex lies on our surfaces. Polia inferior's far
  vertices lie on its 6 skipped fillet (216) faces.

## Classes 211 (revolution) and 216 (fillet) — decoded, pending kernel support (KERNEL_SPEC.md)
- 211: refs the profile (edge curve, p-curve, 100 3D line or 130 3D B-spline). Then o X Y Z (12 f64), box
  (4 f64), two ints, k, 3 f64, int, f64. The last f64 is the middle of the s range, where full turns are
  split into two faces. p(t, s) = o + Rot_Z(s/k)(C(t) - o).
- 216: int n, knot vector (degree 5), n sections of 10 values (p1, p2, centre, angle) with their first and
  second derivatives, then a tail starting with the box u0 u1 v0 v1. The face is the circular arc p1→p2
  around the centre; v = v0 is at p1 and v = v1 is at p2.
- Knot vectors in all compact layouts may list their multiplicities. In raw-integer streams the list
  follows a 0 byte; in compact streams it follows directly; a 1 byte means standard multiplicities.
- 130 (3D B-spline curve, varint): knot vector, i16 rational, poles[3], [weights], tail. 100 (3D line):
  o d, 1.0, t0 t1.
- 170 = compact code 33 (B-spline p-curve): i16 0, knot vector, 2D poles, tail (two i16, t0 t1, 1.0,
  0.0, u8). 2009Fall bodies fit; 2012+ bodies have more fields (open).
- With `--pending`: the 108 parts that have 211/216/251 decode 57,211 faces instead of 55,686, and 36 of
  them become complete. The default mode is unchanged across the corpus (477,290 faces).
