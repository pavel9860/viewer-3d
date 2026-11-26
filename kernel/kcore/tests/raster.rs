//! kcore::raster: the frame does not depend on how triangles are split into chunks or rows into bands, and a closed
//! mesh leaves no crack.
use kcore::raster::*;

/// UV sphere of radius r, outward counter-clockwise triangles; positions, triangles, unit normals
fn sphere(r: f64, nu: usize, nv: usize) -> (Vec<f64>, Vec<u32>, Vec<f64>) {
    let mut v = vec![0.0, 0.0, r];
    for i in 1..nv {
        let th = std::f64::consts::PI * i as f64 / nv as f64;
        for j in 0..nu {
            let ph = 2.0 * std::f64::consts::PI * j as f64 / nu as f64;
            v.extend([r * th.sin() * ph.cos(), r * th.sin() * ph.sin(), r * th.cos()]);
        }
    }
    v.extend([0.0, 0.0, -r]);
    let (top, bot) = (0u32, (v.len() / 3 - 1) as u32);
    let ring = |i: usize, j: usize| (1 + (i - 1) * nu + j % nu) as u32;
    let mut t = vec![];
    for j in 0..nu { t.extend([top, ring(1, j), ring(1, j + 1)]); }
    for i in 1..nv - 1 {
        for j in 0..nu {
            t.extend([ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)]);
            t.extend([ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)]);
        }
    }
    for j in 0..nu { t.extend([ring(nv - 1, j), bot, ring(nv - 1, j + 1)]); }
    let mut n = vec![];
    for k in 0..t.len() / 3 {
        let p = |q: usize| [v[3 * t[3 * k + q] as usize], v[3 * t[3 * k + q] as usize + 1], v[3 * t[3 * k + q] as usize + 2]];
        let (a, b, c) = (p(0), p(1), p(2));
        let (u, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let x = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        n.extend(x.map(|q| if l > 0.0 { q / l } else { 0.0 }));
    }
    (v, t, n)
}

struct Out { color: Vec<u32>, depth: Vec<f32>, id: Vec<i32> }

/// render with the triangles split into `nch` chunks and rows into bands of `bh`
fn render(v: &[f64], t: &[u32], n: &[f64], cam: Cam, nch: usize, bh: u32) -> Out {
    let (nv, nt) = (v.len() / 3, t.len() / 3);
    let c = vec![0.6f32; 3 * nt];
    let m = Mesh { v, t, n, c: &c };
    let mut cl = vec![Cluster::default(); nt.div_ceil(CL)];
    clusters(&m, &mut cl);
    let nb = cam.h.div_ceil(bh) as usize;
    let cam = Cam { bh, nb: nb as u32, ..cam };
    let ncl = nt.div_ceil(CL);
    let ranges: Vec<(usize, usize)> = (0..nch).map(|i| ((CL * (ncl * i / nch)).min(nt), (CL * (ncl * (i + 1) / nch)).min(nt))).collect();
    let mut store: Vec<(Vec<Rec>, Vec<u32>, Vec<u32>, Vec<u32>, u32)> =
        ranges.iter().map(|&(a, b)| (vec![Rec::default(); RECS_PER_TRI * (b - a)], vec![0; nb], vec![0; nb], vec![0; 2 * (8 * (b - a) * nb + 64)], 0)).collect();
    let mut chunks = vec![];
    for (s, &(a, b)) in store.iter_mut().zip(&ranges) {
        let mut cache = vec![Vtx::default(); nv];
        let mut bins = Bins { rec: &mut s.0, head: &mut s.1, tail: &mut s.2, pool: &mut s.3, need: &mut s.4 };
        bin(&m, &cl, &cam, &mut cache, &mut bins, a, b);
        assert!(*bins.need as usize <= bins.pool.len() / 2, "pool too small in the test");
        chunks.push(bins);
    }
    let npx = (cam.w * cam.h) as usize;
    let mut o = Out { color: vec![0; npx], depth: vec![0.0; npx], id: vec![0; npx] };
    let mut fr = Frame { color: &mut o.color, depth: &mut o.depth, id: &mut o.id };
    for band in 0..nb { raster(&cam, &chunks, &mut fr, band); }
    o
}

fn cam(ax: f64, ay: f64, dist: f64, w: u32, h: u32, flags: u32) -> Cam {
    let (cx, sx, cy, sy) = (ax.cos(), ax.sin(), ay.cos(), ay.sin());
    Cam { r: [cy, 0.0, sy, sx * sy, cx, -sx * cy, -cx * sy, sx, cx * cy], t: [0.0, 0.0, -dist], f: (h as f64 / 2.0) / 0.4f64.tan(), w, h, flags, frame: 1, bh: 16, nb: 0 }
}

#[test]
fn split_independent() {
    let (v, t, n) = sphere(0.8, 40, 23);
    for (k, &(ax, ay, dist, flags)) in [(0.3, 0.2, 2.5, 7), (1.1, -0.7, 1.0, 7), (0.0, 0.0, 0.5, 6), (2.0, 1.0, 3.0, 3)].iter().enumerate() {
        let c = cam(ax, ay, dist, 173, 129, flags);
        let a = render(&v, &t, &n, c, 1, 16);
        for (nch, bh) in [(3, 16), (4, 5), (7, 1), (2, 200)] {
            let b = render(&v, &t, &n, c, nch, bh);
            assert!(a.color == b.color && a.id == b.id && a.depth == b.depth, "view {k}: {nch} chunks, bands of {bh} rows differ");
        }
    }
}

#[test]
fn closed_sphere_has_no_crack() {
    let (v, t, n) = sphere(0.8, 57, 31);
    let (w, h, dist) = (201u32, 151u32, 2.5);
    for &(ax, ay) in &[(0.0, 0.0), (0.4, 1.3), (-1.2, 0.5), (1.5707963, 0.0)] {
        let c = cam(ax, ay, dist, w, h, 7);
        let o = render(&v, &t, &n, c, 3, 16);
        // the silhouette of the inscribed polyhedron contains the disk of radius r cos(pi / nu) (and of the rings'
        // angular spacing) seen from the eye: every pixel centre well inside must be drawn
        let r_in = 0.8 * (std::f64::consts::PI / 57.0).cos() * (std::f64::consts::PI / 31.0).cos();
        let rad = c.f * r_in / (dist * dist - r_in * r_in).sqrt() - 1.0;
        for y in 0..h {
            for x in 0..w {
                let (dx, dy) = (x as f64 + 0.5 - w as f64 / 2.0, y as f64 + 0.5 - h as f64 / 2.0);
                if dx * dx + dy * dy < rad * rad {
                    assert!(o.id[(y * w + x) as usize] >= 0, "crack at ({x}, {y}) for view ({ax}, {ay})");
                }
            }
        }
    }
}
