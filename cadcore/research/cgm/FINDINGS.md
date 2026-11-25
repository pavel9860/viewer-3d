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
