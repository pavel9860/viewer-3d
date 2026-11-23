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
