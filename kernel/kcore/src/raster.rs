//! Exact CPU triangle rasterizer for kernel meshes (browser viewer without WebGL, and the tested reference path).
//!
//! Coverage is exact: vertices are snapped to 1/16 pixel and edge functions are evaluated on those integers (exact in
//! f64 for every coordinate a front vertex can reach), with a top-left fill rule, so the triangles of a watertight
//! mesh cover every pixel centre exactly once. Depth is 1/z interpolated in screen space (perspective-correct),
//! compared and stored as f32. Triangles crossing the near plane are clipped against it, from a canonical vertex
//! order so both triangles of a shared edge get the identical clip point. Closed meshes skip back faces, decided on
//! the snapped screen area (the same data as coverage).
//!
//! Parallel structure (sort-middle), two phases with a barrier between them:
//!   1. `bin`: each worker takes a contiguous range of triangles ("chunk"), transforms their vertices, clips, culls,
//!      shades, and appends one record per screen triangle to its chunk, indexed into the bands it overlaps.
//!   2. `raster`: each worker owns a set of bands (rows `band*BH .. +BH`, interleaved across workers for balance),
//!      clears them and draws the records binned there, chunk by chunk, in triangle order.
//! Every pixel therefore sees all its triangles in the original order: the frame is identical for any number of
//! workers and bands (the first triangle with the greatest depth wins).

pub const SUB: f64 = 16.0;
pub const NEAR: f64 = 1e-3;
pub const BACKGROUND: u32 = 0xff21_1d1b;
const EDGE: u32 = 0xff10_0f0d;

/// per-frame camera and options, shared by all workers
#[derive(Clone, Copy)]
#[repr(C)]
pub struct Cam {
    /// rows of the world -> camera rotation
    pub r: [f64; 9],
    /// translation after rotation; the camera looks along -z
    pub t: [f64; 3],
    /// focal length in pixels
    pub f: f64,
    pub w: u32,
    pub h: u32,
    /// bit 0: cull back faces, bit 1: edge lines, bit 2: face colours
    pub flags: u32,
    pub frame: u32,
    /// band height in rows and number of bands
    pub bh: u32,
    pub nb: u32,
}

/// the mesh, normalized: positions, triangles, unit normals and base colours per triangle
pub struct Mesh<'a> {
    pub v: &'a [f64],
    pub t: &'a [u32],
    pub n: &'a [f64],
    pub c: &'a [f32],
}

/// per-worker vertex cache entry, valid when `stamp` equals the frame number: snapped screen position (1/16 px)
/// and 1/z, or iz = 0 when the vertex is not in front of the near plane (its camera position is then recomputed
/// for clipping, which is rare)
#[derive(Clone, Copy, Default)]
#[repr(C)]
pub struct Vtx {
    pub sx: i32,
    pub sy: i32,
    pub iz: f64,
    pub stamp: u32,
    pub pad: u32,
}
/// snapped coordinates are clamped to +-2^26 (4M pixels): edge-function products stay exact in f64
const SNAP_MAX: f64 = (1 << 26) as f64;

/// a run of `CL` consecutive triangles: bounding sphere and normal cone (unit axis, cos and sin of the half-angle;
/// cos = -1 when the normals span more than a hemisphere)
#[derive(Clone, Copy, Default)]
#[repr(C)]
pub struct Cluster {
    pub c: [f64; 3],
    pub rad: f64,
    pub axis: [f64; 3],
    pub cos: f64,
    pub sin: f64,
    pub pad: f64,
}
pub const CL: usize = 32;
/// a cluster is culled as back-facing only with this margin (in cos of the view angle) to spare: triangles nearly
/// edge-on are left to the exact per-triangle test on snapped coordinates
const CONE_MARGIN: f64 = 0.05;
/// bounding boxes narrower than this many pixels are drawn by testing every pixel instead of computing spans
const NARROW: i64 = 16;

/// a screen triangle, counter-clockwise (positive area, y down), 1/16-pixel coordinates, 1/z depths
#[derive(Clone, Copy, Default)]
#[repr(C)]
pub struct Rec {
    pub x: [f64; 3],
    pub y: [f64; 3],
    pub z: [f64; 3],
    pub area: f64,
    pub tid: u32,
    pub rgb: u32,
    /// pixels whose centre can be inside, clamped to the frame: columns minx..=maxx, rows miny..=maxy
    pub minx: i32,
    pub maxx: i32,
    pub miny: i32,
    pub maxy: i32,
}

/// one chunk's binning output: records, and per band a list (in record order) of the records overlapping it.
/// The lists are linked through `pool` entries [record, next]; `need` reports the entries the frame required, so
/// the caller can grow `pool` and bin again when it was too small (`need > pool.len() / 2`).
pub struct Bins<'a> {
    pub rec: &'a mut [Rec],
    pub head: &'a mut [u32],
    pub tail: &'a mut [u32],
    pub pool: &'a mut [u32],
    pub need: &'a mut u32,
}
const NIL: u32 = u32::MAX;

pub struct Frame<'a> {
    pub color: &'a mut [u32],
    pub depth: &'a mut [f32],
    pub id: &'a mut [i32],
}

#[inline(always)]
fn floor(x: f64) -> f64 {
    let t = x as i64 as f64;
    if t > x { t - 1.0 } else { t }
}
#[inline(always)]
fn round(x: f64) -> f64 { floor(x + 0.5) }
#[inline(always)]
fn min3(a: f64, b: f64, c: f64) -> f64 { a.min(b).min(c) }
#[inline(always)]
fn max3(a: f64, b: f64, c: f64) -> f64 { a.max(b).max(c) }

/// records needed per triangle in the worst case (a near-clipped quad gives two)
pub const RECS_PER_TRI: usize = 2;

/// camera-space position of vertex i (z = distance in front of the eye)
#[inline(always)]
fn cam_pos(m: &Mesh, cam: &Cam, i: usize) -> (f64, f64, f64) {
    let (r, t) = (&cam.r, &cam.t);
    let (x, y, z) = (m.v[3 * i], m.v[3 * i + 1], m.v[3 * i + 2]);
    (r[0] * x + r[1] * y + r[2] * z + t[0], r[3] * x + r[4] * y + r[5] * z + t[1], -(r[6] * x + r[7] * y + r[8] * z + t[2]))
}

#[inline(always)]
fn snap(x: f64) -> i32 { round(x).max(-SNAP_MAX).min(SNAP_MAX) as i32 }

#[inline(always)]
fn vertex(m: &Mesh, cam: &Cam, c: &mut [Vtx], i: usize) -> Vtx {
    let c = &mut c[i];
    if c.stamp == cam.frame { return *c; }
    c.stamp = cam.frame;
    let (cx, cy, cz) = cam_pos(m, cam, i);
    if cz > NEAR {
        c.iz = 1.0 / cz;
        let s = cam.f * c.iz;
        c.sx = snap((cam.w as f64 / 2.0 + s * cx) * SUB);
        c.sy = snap((cam.h as f64 / 2.0 - s * cy) * SUB);
    } else {
        c.iz = 0.0;
    }
    *c
}

#[inline(always)]
fn shade(m: &Mesh, cam: &Cam, t: usize) -> u32 {
    let (r, n) = (&cam.r, &m.n[3 * t..3 * t + 3]);
    let nx = r[0] * n[0] + r[1] * n[1] + r[2] * n[2];
    let ny = r[3] * n[0] + r[4] * n[1] + r[5] * n[2];
    let nz = r[6] * n[0] + r[7] * n[1] + r[8] * n[2];
    let l = nx * (0.4 / 1.2329) + ny * (0.6 / 1.2329) + nz * (1.0 / 1.2329);
    let d = 0.25 + 0.75 * l.abs();
    let (cr, cg, cb) = if cam.flags & 4 != 0 { (m.c[3 * t] as f64, m.c[3 * t + 1] as f64, m.c[3 * t + 2] as f64) } else { (0.8, 0.8, 0.8) };
    let q = |v: f64| ((v * d * 255.0) as u32).min(255);
    0xff00_0000 | q(cb) << 16 | q(cg) << 8 | q(cr)
}

/// draw order for culling: triangles grouped by B-rep face, then along a Morton curve of their centroids, so each
/// cluster of `CL` consecutive triangles is compact with a narrow normal cone. Fills `perm` (new -> old index);
/// `keys` is scratch of the same length.
pub fn order(v: &[f64], t: &[u32], face: &[u32], keys: &mut [u128], perm: &mut [u32]) {
    let nt = t.len() / 3;
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for i in 0..v.len() / 3 {
        for q in 0..3 { lo[q] = lo[q].min(v[3 * i + q]); hi[q] = hi[q].max(v[3 * i + q]); }
    }
    let spread = |x: u64| {
        // 10 bits -> every third bit
        let mut x = x & 0x3ff;
        x = (x | x << 16) & 0x0300_00ff;
        x = (x | x << 8) & 0x0300_f00f;
        x = (x | x << 4) & 0x030c_30c3;
        (x | x << 2) & 0x0924_9249
    };
    for k in 0..nt {
        let mut q = [0u64; 3];
        for (d, qd) in q.iter_mut().enumerate() {
            let c = (v[3 * t[3 * k] as usize + d] + v[3 * t[3 * k + 1] as usize + d] + v[3 * t[3 * k + 2] as usize + d]) / 3.0;
            let u = if hi[d] > lo[d] { (c - lo[d]) / (hi[d] - lo[d]) } else { 0.0 };
            *qd = ((u * 1023.0) as u64).min(1023);
        }
        let morton = spread(q[0]) | spread(q[1]) << 1 | spread(q[2]) << 2;
        keys[k] = (face[k] as u128) << 64 | (morton as u128) << 32 | k as u128;
    }
    keys[..nt].sort_unstable();
    for k in 0..nt { perm[k] = keys[k] as u32; }
}

/// renumber vertices in order of first use by `t` (rewritten in place), so the per-frame vertex cache is read
/// nearly sequentially; fills `old` (new -> old index, unused vertices last) using `new` (old -> new) as scratch
pub fn renumber(t: &mut [u32], new: &mut [u32], old: &mut [u32]) {
    new.fill(u32::MAX);
    let mut n = 0u32;
    for i in t.iter_mut() {
        let v = *i as usize;
        if new[v] == u32::MAX {
            new[v] = n;
            old[n as usize] = v as u32;
            n += 1;
        }
        *i = new[v];
    }
    for v in 0..new.len() {
        if new[v] == u32::MAX {
            new[v] = n;
            old[n as usize] = v as u32;
            n += 1;
        }
    }
}

/// bounding spheres and normal cones of the clusters of `CL` consecutive triangles (computed once per mesh)
pub fn clusters(m: &Mesh, out: &mut [Cluster]) {
    let nt = m.t.len() / 3;
    for (k, cl) in out.iter_mut().enumerate().take(nt.div_ceil(CL)) {
        let (t0, t1) = (k * CL, ((k + 1) * CL).min(nt));
        let (mut lo, mut hi, mut ax) = ([f64::MAX; 3], [f64::MIN; 3], [0.0; 3]);
        for t in t0..t1 {
            for j in 0..3 {
                let v = 3 * m.t[3 * t + j] as usize;
                for q in 0..3 { lo[q] = lo[q].min(m.v[v + q]); hi[q] = hi[q].max(m.v[v + q]); }
            }
            for q in 0..3 { ax[q] += m.n[3 * t + q]; }
        }
        let c = [(lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0, (lo[2] + hi[2]) / 2.0];
        let mut r2: f64 = 0.0;
        for t in t0..t1 {
            for j in 0..3 {
                let v = 3 * m.t[3 * t + j] as usize;
                let (dx, dy, dz) = (m.v[v] - c[0], m.v[v + 1] - c[1], m.v[v + 2] - c[2]);
                r2 = r2.max(dx * dx + dy * dy + dz * dz);
            }
        }
        // the sphere is grown slightly: it must contain every vertex despite rounding
        let rad = crate::own_sqrt(r2) * (1.0 + 1e-9) + 1e-12;
        let l = crate::own_sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
        let (mut cos, mut axis): (f64, [f64; 3]) = (-1.0, [0.0; 3]);
        if l > 0.0 {
            axis = [ax[0] / l, ax[1] / l, ax[2] / l];
            cos = 1.0;
            for t in t0..t1 {
                let n = &m.n[3 * t..3 * t + 3];
                // degenerate triangles have a zero normal: they never cover a pixel centre, any cone holds them
                if n[0] == 0.0 && n[1] == 0.0 && n[2] == 0.0 { continue; }
                cos = cos.min(n[0] * axis[0] + n[1] * axis[1] + n[2] * axis[2]);
            }
        }
        let sin = if cos > 0.0 { crate::own_sqrt((1.0 - cos * cos).max(0.0)) } else { 1.0 };
        *cl = Cluster { c, rad, axis, cos, sin, pad: 0.0 };
    }
}

/// the cluster cannot cover a pixel centre: outside the view (with a 2-pixel margin) or behind the near plane, or
/// (with back-face culling) every normal faces away from the eye by more than the margin
#[inline(always)]
fn hidden(cam: &Cam, eye: &[f64; 3], kx: [f64; 2], ky: [f64; 2], cl: &Cluster, cull: bool) -> bool {
    let (r, t, c) = (&cam.r, &cam.t, &cl.c);
    let cx = r[0] * c[0] + r[1] * c[1] + r[2] * c[2] + t[0];
    let cy = r[3] * c[0] + r[4] * c[1] + r[5] * c[2] + t[1];
    let cz = -(r[6] * c[0] + r[7] * c[1] + r[8] * c[2] + t[2]);
    let rad = cl.rad;
    if cz + rad <= NEAR { return true; }
    // side planes |x| <= k z, written as (x - k z) / sqrt(1 + k^2) <= 0 with k[1] = 1 / sqrt(1 + k^2)
    if (cx - kx[0] * cz) * kx[1] > rad || (-cx - kx[0] * cz) * kx[1] > rad { return true; }
    if (cy - ky[0] * cz) * ky[1] > rad || (-cy - ky[0] * cz) * ky[1] > rad { return true; }
    if !cull || cl.cos <= 0.0 { return false; }
    // back-facing: n . (p - eye) > 0 for every normal n in the cone and point p in the sphere; with d = c - eye at
    // angle a to the axis and cone half-angle b, min n . d = |d| cos(a + b) and n . (p - c) >= -rad
    let d = [c[0] - eye[0], c[1] - eye[1], c[2] - eye[2]];
    let l2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    let dt = d[0] * cl.axis[0] + d[1] * cl.axis[1] + d[2] * cl.axis[2];
    if dt <= 0.0 { return false; }
    let l = crate::own_sqrt(l2);
    let cr = crate::own_sqrt((l2 - dt * dt).max(0.0));
    dt * cl.cos - cr * cl.sin > rad + CONE_MARGIN * l
}

/// phase 1: transform, clip, cull and bin triangles [t0, t1); `t0` is a multiple of `CL`
pub fn bin(m: &Mesh, cl: &[Cluster], cam: &Cam, c: &mut [Vtx], b: &mut Bins, t0: usize, t1: usize) {
    let cull = cam.flags & 1 != 0;
    let (cxs, cys) = (cam.w as f64 / 2.0, cam.h as f64 / 2.0);
    let (wl, hl) = (cam.w as i32 - 1, cam.h as i32 - 1);
    let (bh, nb) = (cam.bh as i32, cam.nb as usize);
    let (r, tr) = (&cam.r, &cam.t);
    let eye = [-(r[0] * tr[0] + r[3] * tr[1] + r[6] * tr[2]), -(r[1] * tr[0] + r[4] * tr[1] + r[7] * tr[2]), -(r[2] * tr[0] + r[5] * tr[1] + r[8] * tr[2])];
    let kx0 = (cxs + 2.0) / cam.f;
    let ky0 = (cys + 2.0) / cam.f;
    let (kx, ky) = ([kx0, 1.0 / crate::own_sqrt(1.0 + kx0 * kx0)], [ky0, 1.0 / crate::own_sqrt(1.0 + ky0 * ky0)]);
    b.head[..nb].fill(NIL);
    let cap = b.pool.len() / 2;
    let (mut nr, mut ne) = (0usize, 0usize);
    let (mut sx, mut sy, mut iz) = ([0i32; 4], [0i32; 4], [0f64; 4]);
    let mut t = t0;
    while t < t1 {
        let te = ((t / CL + 1) * CL).min(t1);
        if hidden(cam, &eye, kx, ky, &cl[t / CL], cull) {
            t = te;
            continue;
        }
        for t in t..te {
            let (a, bb, d) = (m.t[3 * t] as usize, m.t[3 * t + 1] as usize, m.t[3 * t + 2] as usize);
            let (va, vb, vd) = (vertex(m, cam, c, a), vertex(m, cam, c, bb), vertex(m, cam, c, d));
            let k;
            if va.iz > 0.0 && vb.iz > 0.0 && vd.iz > 0.0 {
                sx = [va.sx, vb.sx, vd.sx, 0];
                sy = [va.sy, vb.sy, vd.sy, 0];
                iz = [va.iz, vb.iz, vd.iz, 0.0];
                k = 3;
            } else {
                if va.iz == 0.0 && vb.iz == 0.0 && vd.iz == 0.0 { continue; }
                // clip against the near plane (Sutherland-Hodgman): at most 4 vertices
                let vs = [va, vb, vd];
                let ix = [a, bb, d];
                let mut n = 0;
                for j in 0..3 {
                    let (p, q) = (vs[j], vs[(j + 1) % 3]);
                    let (pin, qin) = (p.iz > 0.0, q.iz > 0.0);
                    if pin {
                        // the cached snapped position: identical to the unclipped triangles sharing the vertex
                        sx[n] = p.sx;
                        sy[n] = p.sy;
                        iz[n] = p.iz;
                        n += 1;
                    }
                    if pin != qin {
                        // from the lower vertex index: both triangles of a shared edge get the identical point
                        let (u, v) = if ix[j] < ix[(j + 1) % 3] { (ix[j], ix[(j + 1) % 3]) } else { (ix[(j + 1) % 3], ix[j]) };
                        let (ux, uy, uz) = cam_pos(m, cam, u);
                        let (vx, vy, vz) = cam_pos(m, cam, v);
                        let s = (NEAR - uz) / (vz - uz);
                        let sc = cam.f / NEAR;
                        sx[n] = snap((cxs + sc * (ux + s * (vx - ux))) * SUB);
                        sy[n] = snap((cys - sc * (uy + s * (vy - uy))) * SUB);
                        iz[n] = 1.0 / NEAR;
                        n += 1;
                    }
                }
                k = n;
            }
            let mut rgb = 0;
            for j in 1..k.saturating_sub(1) {
                let (x, y) = ([sx[0], sx[j], sx[j + 1]], [sy[0], sy[j], sy[j + 1]]);
                // pixel-centre bounds (centres at 16 i + 8): rows/columns whose centre can be inside
                let minx = (x[0].min(x[1]).min(x[2]) + 7) >> 4;
                let maxx = (x[0].max(x[1]).max(x[2]) - 8) >> 4;
                let miny = ((y[0].min(y[1]).min(y[2]) + 7) >> 4).max(0);
                let maxy = ((y[0].max(y[1]).max(y[2]) - 8) >> 4).min(hl);
                if minx.max(0) > maxx.min(wl) || miny > maxy { continue; }
                let area = (x[1] - x[0]) as i64 * (y[2] - y[0]) as i64 - (x[2] - x[0]) as i64 * (y[1] - y[0]) as i64;
                if area == 0 { continue; }
                // a front face is counter-clockwise seen from the eye: negative area with y down
                if cull && area > 0 { continue; }
                let (mut x, mut y, mut z) = (x.map(|v| v as f64), y.map(|v| v as f64), [iz[0], iz[j], iz[j + 1]]);
                if area < 0 {
                    x.swap(1, 2);
                    y.swap(1, 2);
                    z.swap(1, 2);
                }
                if rgb == 0 { rgb = shade(m, cam, t); }
                b.rec[nr] = Rec { x, y, z, area: area.abs() as f64, tid: t as u32, rgb, minx: minx.max(0), maxx: maxx.min(wl), miny, maxy };
                for band in (miny / bh) as usize..=(maxy / bh) as usize {
                    if ne < cap {
                        b.pool[2 * ne] = nr as u32;
                        b.pool[2 * ne + 1] = NIL;
                        if b.head[band] == NIL { b.head[band] = ne as u32; } else { b.pool[2 * b.tail[band] as usize + 1] = ne as u32; }
                        b.tail[band] = ne as u32;
                    }
                    ne += 1;
                }
                nr += 1;
            }
        }
        t = te;
    }
    *b.need = ne as u32;
}

pub fn clear_rows(fr: &mut Frame, w: usize, y0: usize, y1: usize) {
    let (s, e) = (y0 * w, y1 * w);
    fr.color[s..e].fill(BACKGROUND);
    fr.depth[s..e].fill(0.0);
    fr.id[s..e].fill(-1);
}

/// phase 2: clear band `band` and draw every chunk's records binned there, in chunk order; returns pixel writes
pub fn raster(cam: &Cam, chunks: &[Bins], fr: &mut Frame, band: usize) -> u64 {
    let (w, h) = (cam.w as usize, cam.h as usize);
    let y0 = band * cam.bh as usize;
    let y1 = (y0 + cam.bh as usize).min(h);
    if y0 >= y1 { return 0; }
    clear_rows(fr, w, y0, y1);
    let edges = cam.flags & 2 != 0;
    let mut drawn = 0;
    for b in chunks {
        let mut e = b.head[band];
        while e != NIL {
            drawn += tri(&b.rec[b.pool[2 * e as usize] as usize], edges, w as i64, y0 as i64, y1 as i64 - 1, fr);
            e = b.pool[2 * e as usize + 1];
        }
    }
    drawn
}

#[inline(always)]
fn clamp(z: f64, lo: f64, hi: f64) -> f64 {
    if z < lo { lo } else if z > hi { hi } else { z }
}

/// pixels [s, e) of one row: depth test and write; `ed` gives the edge values at s and their steps (edge lines)
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn span<const EDGES: bool>(fr: &mut Frame, s: usize, e: usize, mut zz: f64, dzdx: f64, zlo: f64, zhi: f64, rgb: u32, tid: i32, ed: &mut [f64; 9]) -> u64 {
    assert!(e <= fr.depth.len() && e <= fr.id.len() && e <= fr.color.len());
    let mut covered = 0;
    let (dep, ids, col) = (fr.depth.as_mut_ptr(), fr.id.as_mut_ptr(), fr.color.as_mut_ptr());
    for o in s..e {
        // a near-degenerate triangle's snapped barycentrics can extrapolate: clamp to its vertex depths
        let zc = clamp(zz, zlo, zhi) as f32;
        // SAFETY: o < e <= every buffer's length (asserted above)
        unsafe {
            if zc > *dep.add(o) {
                *dep.add(o) = zc;
                *ids.add(o) = tid;
                *col.add(o) = if EDGES && (ed[0] * ed[0] < ed[6] || ed[1] * ed[1] < ed[7] || ed[2] * ed[2] < ed[8]) { EDGE } else { rgb };
                covered += 1;
            }
        }
        if EDGES {
            ed[0] += ed[3];
            ed[1] += ed[4];
            ed[2] += ed[5];
        }
        zz += dzdx;
    }
    covered
}

#[inline(always)]
fn tri(r: &Rec, edges: bool, w: i64, ya: i64, yb: i64, fr: &mut Frame) -> u64 {
    if edges { tri_e::<true>(r, w, ya, yb, fr) } else { tri_e::<false>(r, w, ya, yb, fr) }
}

/// one counter-clockwise triangle in 1/16-pixel integer coordinates, rows [ya, yb]: spans from exact edge
/// functions with a top-left rule
#[inline(always)]
fn tri_e<const EDGES: bool>(r: &Rec, w: i64, ya: i64, yb: i64, fr: &mut Frame) -> u64 {
    let ([x0, x1, x2], [y0, y1, y2], [z0, z1, z2]) = (r.x, r.y, r.z);
    let (minx, maxx) = (r.minx as i64, r.maxx as i64);
    let (miny, maxy) = ((r.miny as i64).max(ya), (r.maxy as i64).min(yb));
    if miny > maxy { return 0; }
    // edge k opposite vertex k: E increases by A per pixel in x and by B per row
    let (a0, b0, a1, b1, a2, b2) = ((y1 - y2) * SUB, (x2 - x1) * SUB, (y2 - y0) * SUB, (x0 - x2) * SUB, (y0 - y1) * SUB, (x1 - x0) * SUB);
    // top-left rule: on an edge that is neither top nor left the boundary is outside (bias -1 on exact integers)
    let tl = |dx: f64, dy: f64| if dy < 0.0 || (dy == 0.0 && dx > 0.0) { 0.0 } else { -1.0 };
    let (o0, o1, o2) = (tl(x2 - x1, y2 - y1), tl(x0 - x2, y0 - y2), tl(x1 - x0, y1 - y0));
    let px = minx as f64 * SUB + SUB / 2.0;
    let py = miny as f64 * SUB + SUB / 2.0;
    let mut r0 = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1) + o0;
    let mut r1 = (x0 - x2) * (py - y2) - (y0 - y2) * (px - x2) + o1;
    let mut r2 = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0) + o2;
    let inv = 1.0 / r.area;
    // 0.6 px edge lines: E < 0.6 SUB |edge|, compared squared (E >= 0 inside)
    let el = |dx: f64, dy: f64| 0.36 * SUB * SUB * (dx * dx + dy * dy);
    let (w0, w1, w2) = if EDGES { (el(x2 - x1, y2 - y1), el(x0 - x2, y0 - y2), el(x1 - x0, y1 - y0)) } else { (0.0, 0.0, 0.0) };
    let (z0u, z1u, z2u) = (z0 * inv, z1 * inv, z2 * inv);
    let dzdx = a0 * z0u + a1 * z1u + a2 * z2u;
    let (zlo, zhi) = (min3(z0, z1, z2), max3(z0, z1, z2));
    let (rgb, tid) = (r.rgb, r.tid as i32);
    let wu = w as usize;
    let mut covered = 0u64;
    // depth from the unbiased edge values (the fill-rule bias is not geometry)
    let depth_at = |e0: f64, e1: f64, e2: f64| (e0 - o0) * z0u + (e1 - o1) * z1u + (e2 - o2) * z2u;
    // narrow triangles (most of a fine mesh): the first and last covered pixel of a row by testing the bounding box
    if maxx - minx < NARROW {
        for yy in miny..=maxy {
            let (mut xs, mut xe) = (maxx + 1, minx - 1);
            let (mut e0, mut e1, mut e2) = (r0, r1, r2);
            for xx in minx..=maxx {
                if (e0 >= 0.0) & (e1 >= 0.0) & (e2 >= 0.0) {
                    xs = xs.min(xx);
                    xe = xx;
                }
                e0 += a0;
                e1 += a1;
                e2 += a2;
            }
            if xs <= xe {
                let i = (xs - minx) as f64;
                let mut ed = [r0 + a0 * i, r1 + a1 * i, r2 + a2 * i, a0, a1, a2, w0, w1, w2];
                let row = yy as usize * wu;
                covered += span::<EDGES>(fr, row + xs as usize, row + xe as usize + 1, depth_at(ed[0], ed[1], ed[2]), dzdx, zlo, zhi, rgb, tid, &mut ed);
            }
            r0 += b0;
            r1 += b1;
            r2 += b2;
        }
        return covered;
    }
    // span ends from the reciprocal slopes (approximate; made exact by the widen-and-shrink below)
    let inva = |a: f64| if a != 0.0 { 1.0 / a } else { 0.0 };
    let (i0, i1, i2) = (inva(a0), inva(a1), inva(a2));
    let inside = |r0: f64, r1: f64, r2: f64, i: i64| {
        let i = i as f64;
        (r0 + a0 * i >= 0.0) & (r1 + a1 * i >= 0.0) & (r2 + a2 * i >= 0.0)
    };
    for yy in miny..=maxy {
        let (mut xs, mut xe) = (minx, maxx);
        let mut skip = false;
        for (r, a, ia) in [(r0, a0, i0), (r1, a1, i1), (r2, a2, i2)] {
            // truncation is within one pixel of floor/ceil: the widen-and-shrink below makes the ends exact
            let q = (-r * ia).max(-1e15).min(1e15) as i64;
            if a > 0.0 { xs = xs.max(minx + q); } else if a < 0.0 { xe = xe.min(minx + q); } else if r < 0.0 { skip = true; }
        }
        if !skip {
            // the products may round and truncate either way: widen by one pixel, then shrink with the exact test
            xs = (xs - 1).max(minx);
            xe = (xe + 1).min(maxx);
            while xs <= xe && !inside(r0, r1, r2, xs - minx) { xs += 1; }
            while xe >= xs && !inside(r0, r1, r2, xe - minx) { xe -= 1; }
            if xs <= xe {
                let i = (xs - minx) as f64;
                let mut ed = [r0 + a0 * i, r1 + a1 * i, r2 + a2 * i, a0, a1, a2, w0, w1, w2];
                let row = yy as usize * wu;
                covered += span::<EDGES>(fr, row + xs as usize, row + xe as usize + 1, depth_at(ed[0], ed[1], ed[2]), dzdx, zlo, zhi, rgb, tid, &mut ed);
            }
        }
        r0 += b0;
        r1 += b1;
        r2 += b2;
    }
    covered
}
