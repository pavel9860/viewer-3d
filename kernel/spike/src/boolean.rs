//! Boolean spike on B-rep solids made of exact analytic faces.
//!
//! Pipeline: one SSI per face pair -> curve restricted to both faces, with exact vertices where it crosses or
//! touches either face's edges or seam (vertex registry: one vertex per surface triple) -> per-face planar
//! arrangement in UV (half-edge graph, exact orient2d) -> regions -> classify against the other solid ->
//! select by operation -> verify closure (every directed boundary segment matched) and volume.
use crate::exact::orient2d;
use crate::geom::*;
use crate::ms;
use crate::ssi::{intersect, Kind};
use std::cell::RefCell;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::time::Instant;

/// samples on circular edges (shared by both faces of the edge)
const NC: usize = 256;
/// UV weld tolerance and boundary-contact tolerance
const UV_TOL: f64 = 1e-8;

#[derive(Clone, Copy, Debug)]
enum Nb { Face(usize), SeamU }
/// boundary loop point; `nb` is what the segment starting here borders
#[derive(Clone)]
struct BPt { uv: [f64; 2], p: V, nb: Nb }
/// face: surface + one outer loop (CCW in UV); `per` is the u-period of a periodic face
struct Face { s: Box<dyn Surf>, lp: Vec<BPt>, per: Option<f64> }
/// solid: faces + point-in-solid (negative inside)
struct Solid { faces: Vec<Face>, inside: Box<dyn Fn(V) -> f64> }

// ------------------------------------------------------------------ primitives
fn plane_face(o: V, a: V, b: V, h: f64, nbs: [usize; 4]) -> Face {
    let s = Plane { o, a, b, ext: 1.6 * h };
    let c = [[-h, -h], [h, -h], [h, h], [-h, h]];
    let lp = (0..4).map(|i| BPt { uv: c[i], p: s.ev(c[i][0], c[i][1]).p, nb: Nb::Face(nbs[i]) }).collect();
    Face { s: Box::new(s), lp, per: None }
}
/// axis-aligned cube, centre c, half size h; faces +x -x +y -y +z -z (a x b = outward normal)
fn cube(c: V, h: f64) -> Solid {
    let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let o = |d: V, s: f64| addv(c, mulv(d, s * h));
    let faces = vec![
        plane_face(o(x, 1.0), y, z, h, [5, 2, 4, 3]), plane_face(o(x, -1.0), z, y, h, [3, 4, 2, 5]),
        plane_face(o(y, 1.0), z, x, h, [1, 4, 0, 5]), plane_face(o(y, -1.0), x, z, h, [5, 0, 4, 1]),
        plane_face(o(z, 1.0), x, y, h, [3, 0, 2, 1]), plane_face(o(z, -1.0), y, x, h, [1, 2, 0, 3]),
    ];
    Solid { faces, inside: Box::new(move |p| (0..3).map(|k| (p[k] - c[k]).abs() - h).fold(f64::MIN, f64::max)) }
}
/// capped cylinder: centre c, radius r, half height h, axes (ex, ey, ez); faces lateral, top cap, bottom cap
fn cylinder(c: V, r: f64, h: f64, ex: V, ey: V, ez: V) -> Solid {
    let lat = Moved { s: Box::new(Cyl { r, h, x: ex, y: ey, z: ez }), r: frame([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]), t: c };
    let ang = |i: usize| 2.0 * PI * i as f64 / NC as f64;
    let at = |u: f64, v: f64| lat.ev(u, v).p;
    // lateral loop, CCW in (u,v): bottom u 0->2pi (cap 2), seam up at 2pi, top u 2pi->0 (cap 1), seam down at 0
    let mut lp: Vec<BPt> = (0..NC).map(|i| BPt { uv: [ang(i), -h], p: at(ang(i), -h), nb: Nb::Face(2) }).collect();
    lp.push(BPt { uv: [2.0 * PI, -h], p: at(0.0, -h), nb: Nb::SeamU });
    lp.extend((1..=NC).rev().map(|i| BPt { uv: [ang(i), h], p: at(ang(i), h), nb: Nb::Face(1) }));
    lp.push(BPt { uv: [0.0, h], p: at(0.0, h), nb: Nb::SeamU });
    let mut faces = vec![Face { s: Box::new(lat), lp, per: Some(2.0 * PI) }];
    for (sg, a, b) in [(1.0, ex, ey), (-1.0, ey, ex)] {
        let o = addv(c, mulv(ez, sg * h));
        let pl = Plane { o, a, b, ext: 1.6 * r };
        let mut lp: Vec<BPt> = (0..NC).map(|i| { let p = addv(o, addv(mulv(ex, r * ang(i).cos()), mulv(ey, r * ang(i).sin()))); BPt { uv: pl.uv_of(p), p, nb: Nb::Face(0) } }).collect();
        if sg < 0.0 { lp.reverse(); } // CCW in the swapped (a, b) frame
        faces.push(Face { s: Box::new(pl), lp, per: None });
    }
    Solid { faces, inside: Box::new(move |p| { let q = sub(p, c); let (x, y, z) = (dot(q, ex), dot(q, ey), dot(q, ez)); ((x * x + y * y).sqrt() - r).max(z.abs() - h) }) }
}

// ------------------------------------------------------------------ small dense Gauss-Newton (n <= 6)
fn gn(mut x: Vec<f64>, m: usize, r: &dyn Fn(&[f64], &mut [f64])) -> Vec<f64> {
    let n = x.len();
    let (mut f0, mut f1) = (vec![0.0; m], vec![0.0; m]);
    for _ in 0..60 {
        r(&x, &mut f0);
        let c0: f64 = f0.iter().map(|v| v * v).sum();
        if !(c0 >= 1e-30) { break; }
        let mut j = vec![vec![0.0; n]; m];
        for k in 0..n { let h = 1e-7 * (1.0 + x[k].abs()); let mut xp = x.clone(); xp[k] += h; r(&xp, &mut f1); for i in 0..m { j[i][k] = (f1[i] - f0[i]) / h; } }
        let mut a = vec![vec![0.0; n + 1]; n];
        for p in 0..n { for q in 0..n { a[p][q] = (0..m).map(|i| j[i][p] * j[i][q]).sum::<f64>() + if p == q { 1e-14 } else { 0.0 }; } a[p][n] = -(0..m).map(|i| j[i][p] * f0[i]).sum::<f64>(); }
        for c in 0..n {
            let pv = (c..n).max_by(|&i, &k| a[i][c].abs().total_cmp(&a[k][c].abs())).unwrap();
            a.swap(c, pv);
            if !(a[c][c].abs() > 1e-300) { return x; }
            for rr in 0..n { if rr != c { let f = a[rr][c] / a[c][c]; for k in c..=n { a[rr][k] -= f * a[c][k]; } } }
        }
        let xn: Vec<f64> = (0..n).map(|k| x[k] + a[k][n] / a[k][k]).collect();
        r(&xn, &mut f1);
        if !(f1.iter().map(|v| v * v).sum::<f64>() < c0) { break; }
        x = xn;
    }
    x
}

// ------------------------------------------------------------------ vertices
thread_local! { static VERTS: RefCell<Vec<([usize; 3], V)>> = RefCell::default(); }
fn sid(s: &dyn Surf) -> usize { s as *const dyn Surf as *const u8 as usize }
/// a vertex is one topological entity per surface triple: solved once, reused bit for bit
/// (tangential vertices are ill-conditioned: independent solves differ by ~1e-7)
fn registered(mut key: [usize; 3], q: V) -> V {
    key.sort_unstable();
    VERTS.with(|v| {
        let mut v = v.borrow_mut();
        if let Some((_, p)) = v.iter().find(|(k, p)| *k == key && norm(sub(*p, q)) < 1e-5) { return *p; }
        v.push((key, q));
        q
    })
}
/// exact vertex where the curve F ∩ S meets the boundary segment of F bordering `nb`
fn vertex(f: &Face, solid: &Solid, nb: Nb, s: &dyn Surf, guess: V) -> V {
    match nb {
        Nb::Face(g) => {
            let gs = solid.faces[g].s.as_ref();
            let (a, b, c) = (f.s.uv_of(guess), gs.uv_of(guess), s.uv_of(guess));
            let x = gn(vec![a[0], a[1], b[0], b[1], c[0], c[1]], 6, &|x, r| {
                let (pf, pg, ps) = (f.s.ev(x[0], x[1]).p, gs.ev(x[2], x[3]).p, s.ev(x[4], x[5]).p);
                r[..3].copy_from_slice(&sub(pf, pg)); r[3..].copy_from_slice(&sub(pf, ps));
            });
            registered([sid(f.s.as_ref()), sid(gs), sid(s)], f.s.ev(x[0], x[1]).p)
        }
        Nb::SeamU => vertex_seam(f, s, guess),
    }
}
fn vertex_seam(f: &Face, s: &dyn Surf, guess: V) -> V {
    let per = f.per.unwrap();
    let (a, c) = (f.s.uv_of(guess), s.uv_of(guess));
    let u0 = if a[0] > per / 2.0 { per } else { 0.0 };
    let x = gn(vec![a[1], c[0], c[1]], 3, &|x, r| r.copy_from_slice(&sub(f.s.ev(u0, x[0]).p, s.ev(x[1], x[2]).p)));
    registered([sid(f.s.as_ref()), sid(f.s.as_ref()) + 1, sid(s)], f.s.ev(u0, x[0]).p)
}
fn on_surface(s: &dyn Surf, q: V) -> bool { let uv = s.uv_of(q); norm(sub(s.ev(uv[0], uv[1]).p, q)) < 1e-9 }

// ------------------------------------------------------------------ 2D helpers (exact predicates)
fn o2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> i8 { orient2d(a, b, c).0 }
/// proper crossing of ab with cd; parameter along ab
fn seg_cross(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> Option<f64> {
    let (d1, d2, d3, d4) = (o2(c, d, a), o2(c, d, b), o2(a, b, c), o2(a, b, d));
    let crosses = (d1 * d2 < 0 && d3 * d4 < 0) || (d1 == 0 && d2 != 0 && d3 * d4 <= 0) || (d2 == 0 && d1 != 0 && d3 * d4 <= 0);
    if !crosses { return None; }
    let den = (b[0] - a[0]) * (d[1] - c[1]) - (b[1] - a[1]) * (d[0] - c[0]);
    (den != 0.0).then(|| ((c[0] - a[0]) * (d[1] - c[1]) - (c[1] - a[1]) * (d[0] - c[0])) / den)
}
fn pt_seg(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
    ((a[0] + t * dx - p[0]).powi(2) + (a[1] + t * dy - p[1]).powi(2)).sqrt()
}
fn seg_dist(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> f64 {
    if seg_cross(a, b, c, d).is_some() { return 0.0; }
    pt_seg(a, c, d).min(pt_seg(b, c, d)).min(pt_seg(c, a, b)).min(pt_seg(d, a, b))
}
fn inside_poly(p: [f64; 2], lp: &[[f64; 2]]) -> bool {
    let mut c = false;
    for i in 0..lp.len() { let (a, b) = (lp[i], lp[(i + 1) % lp.len()]); if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) { c = !c; } }
    c
}
fn area(pts: &[[f64; 2]], cyc: &[usize]) -> f64 { (0..cyc.len()).map(|i| { let (p, q) = (pts[cyc[i]], pts[cyc[(i + 1) % cyc.len()]]); p[0] * q[1] - p[1] * q[0] }).sum::<f64>() / 2.0 }

// ------------------------------------------------------------------ shared edges per face pair
/// uv of a polyline on face f, unwrapped across the period and shifted into [0, per]
fn uv_run(f: &Face, pts: &[V]) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = vec![];
    for &p in pts {
        let mut uv = f.s.uv_of(p);
        if let (Some(per), Some(l)) = (f.per, out.last()) { while uv[0] - l[0] > per / 2.0 { uv[0] -= per; } while l[0] - uv[0] > per / 2.0 { uv[0] += per; } }
        out.push(uv);
    }
    if let Some(per) = f.per {
        let k = (out.iter().map(|q| q[0]).sum::<f64>() / out.len() as f64 / per).floor();
        for q in out.iter_mut() { q[0] = (q[0] - k * per).clamp(0.0, per); }
    }
    out
}
fn in_face(f: &Face, mut uv: [f64; 2]) -> bool {
    if let Some(per) = f.per { uv[0] = uv[0].rem_euclid(per); }
    inside_poly(uv, &f.lp.iter().map(|b| b.uv).collect::<Vec<_>>())
}

/// an edge piece shared by two faces: 3D polyline + indices that are seam points of face a / face b
struct Piece { p: Vec<V>, seam: [Vec<usize>; 2] }

/// the SSI curves of fa and fb restricted to both faces, with exact vertices where a curve crosses or touches
/// either face's boundary or seam
fn face_pair_edges(fa: &Face, sa: &Solid, fb: &Face, sb: &Solid) -> Vec<Piece> {
    let mut out = vec![];
    for c in intersect(fa.s.as_ref(), fb.s.as_ref()).0 {
        if c.kind == Kind::Isolated || c.pts.len() < 2 { continue; }
        let mut pts = c.pts;
        let mut seam = [vec![false; pts.len()], vec![false; pts.len()]];
        let insert = |pts: &mut Vec<V>, seam: &mut [Vec<bool>; 2], at: usize, q: V| { pts.insert(at, q); seam[0].insert(at, false); seam[1].insert(at, false); };
        // a) seam crossings of each periodic face (an existing point on the seam is the crossing itself)
        for (fi, (f, other)) in [(fa, fb), (fb, fa)].into_iter().enumerate() {
            let Some(per) = f.per else { continue };
            let on_seam = |u: f64| u < 1e-9 || per - u < 1e-9;
            let mut k = 0;
            while k + 1 < pts.len() {
                let (u0, u1) = (f.s.uv_of(pts[k])[0], f.s.uv_of(pts[k + 1])[0]);
                if (u1 - u0).abs() > per / 2.0 {
                    if on_seam(u0) || on_seam(u1) { seam[fi][if on_seam(u0) { k } else { k + 1 }] = true; k += 1; continue; }
                    let q = vertex_seam(f, other.s.as_ref(), mulv(addv(pts[k], pts[k + 1]), 0.5));
                    insert(&mut pts, &mut seam, k + 1, q);
                    seam[fi][k + 1] = true;
                    k += 2;
                } else { k += 1; }
            }
        }
        // b) proper crossings of either face's boundary; b2) touch points (reaches an edge without crossing)
        for (f, solid, other) in [(fa, sa, fb), (fb, sb, fa)] {
            let n = f.lp.len();
            let mut k = 0;
            while k + 1 < pts.len() {
                let uv = uv_run(f, &pts[k..=k + 1]);
                let mut hit = None;
                for bi in 0..n {
                    let (b0, b1) = (f.lp[bi].uv, f.lp[(bi + 1) % n].uv);
                    if let Some(t) = seg_cross(uv[0], uv[1], b0, b1) {
                        if t > 1e-9 && t < 1.0 - 1e-9 { hit = Some(vertex(f, solid, f.lp[bi].nb, other.s.as_ref(), addv(pts[k], mulv(sub(pts[k + 1], pts[k]), t)))); break; }
                    } else if seg_dist(uv[0], uv[1], b0, b1) <= 1e-4 && pt_seg(uv[0], b0, b1) >= 1e-9 && pt_seg(uv[1], b0, b1) >= 1e-9 {
                        let q = vertex(f, solid, f.lp[bi].nb, other.s.as_ref(), mulv(addv(pts[k], pts[k + 1]), 0.5));
                        let on3 = on_surface(f.s.as_ref(), q) && on_surface(other.s.as_ref(), q) && match f.lp[bi].nb { Nb::Face(g) => on_surface(solid.faces[g].s.as_ref(), q), Nb::SeamU => true };
                        let seg = sub(pts[k + 1], pts[k]);
                        let t = dot(sub(q, pts[k]), seg) / dot(seg, seg);
                        let quv = uv_run(f, &[pts[k], q])[1];
                        if on3 && pt_seg(quv, b0, b1) < 1e-9 && t > 1e-9 && t < 1.0 - 1e-9 { hit = Some(q); break; }
                    }
                }
                if let Some(q) = hit { insert(&mut pts, &mut seam, k + 1, q); k += 1; }
                k += 1;
            }
        }
        // c) maximal runs of segments inside both faces
        let mut run: Vec<usize> = vec![];
        let flush = |run: &mut Vec<usize>, out: &mut Vec<Piece>| {
            if run.len() >= 2 {
                let idx = |fi: usize| run.iter().enumerate().filter(|(_, &i)| seam[fi][i]).map(|(k, _)| k).collect();
                out.push(Piece { p: run.iter().map(|&i| pts[i]).collect(), seam: [idx(0), idx(1)] });
            }
            run.clear();
        };
        for k in 0..pts.len() - 1 {
            let mid = mulv(addv(pts[k], pts[k + 1]), 0.5);
            if in_face(fa, uv_run(fa, &[pts[k], mid])[1]) && in_face(fb, uv_run(fb, &[pts[k], mid])[1]) { if run.is_empty() { run.push(k); } run.push(k + 1); } else { flush(&mut run, &mut out); }
        }
        flush(&mut run, &mut out);
    }
    out
}
/// per-face cut polylines (uv + 3D), split at this face's seam points
fn cuts_for(f: &Face, p: &[V], seam_idx: &[usize]) -> Vec<Vec<([f64; 2], V)>> {
    let mut bounds = seam_idx.to_vec();
    bounds.push(p.len() - 1);
    let mut start = 0;
    let mut cuts = vec![];
    for e in bounds {
        if e > start { cuts.push(uv_run(f, &p[start..=e]).into_iter().zip(p[start..=e].iter().copied()).collect()); }
        start = e;
    }
    cuts
}

// ------------------------------------------------------------------ arrangement
struct Region { outer: Vec<usize>, holes: Vec<Vec<usize>> }
struct Arr { uv: Vec<[f64; 2]>, p: Vec<V>, regions: Vec<Region> }

fn arrange(f: &Face, cuts: Vec<Vec<([f64; 2], V)>>) -> Arr {
    let n = f.lp.len();
    let mut bsegs: Vec<(BPt, BPt)> = (0..n).map(|i| (f.lp[i].clone(), f.lp[(i + 1) % n].clone())).collect();
    // every cut vertex on the boundary splits it; seam vertices are mirrored onto the other side of the seam
    for c in &cuts {
        let mut pts = c.clone();
        if let Some(per) = f.per { for &(q, p) in c { if q[0].abs() < 1e-9 { pts.push(([per, q[1]], p)); } else if (per - q[0]).abs() < 1e-9 { pts.push(([0.0, q[1]], p)); } } }
        for (q, p) in pts {
            let best = bsegs.iter().enumerate().filter_map(|(bi, (a, b))| {
                let (dx, dy) = (b.uv[0] - a.uv[0], b.uv[1] - a.uv[1]);
                let l2 = dx * dx + dy * dy;
                let t = ((q[0] - a.uv[0]) * dx + (q[1] - a.uv[1]) * dy) / l2;
                (l2 > 0.0 && t > 1e-12 && t < 1.0 - 1e-12).then(|| (bi, pt_seg(q, a.uv, b.uv)))
            }).filter(|x| x.1 < UV_TOL).min_by(|x, y| x.1.total_cmp(&y.1));
            if let Some((bi, _)) = best {
                let (a, b) = bsegs[bi].clone();
                let mid = BPt { uv: q, p, nb: a.nb };
                bsegs[bi] = (a, mid.clone());
                bsegs.insert(bi + 1, (mid, b));
            }
        }
    }
    // a cut segment lying along the boundary is that boundary edge itself
    let on_bnd = |q: [f64; 2]| bsegs.iter().any(|(a, b)| pt_seg(q, a.uv, b.uv) < 1e-9);
    let csegs: Vec<_> = cuts.iter().flat_map(|c| c.windows(2).filter(|w| !(on_bnd(w[0].0) && on_bnd(w[1].0) && on_bnd([(w[0].0[0] + w[1].0[0]) / 2.0, (w[0].0[1] + w[1].0[1]) / 2.0]))).map(|w| (w[0], w[1])).collect::<Vec<_>>()).collect();
    // vertices welded in UV (grid hash)
    let (mut uv, mut p3): (Vec<[f64; 2]>, Vec<V>) = (vec![], vec![]);
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    let mut vid = |q: [f64; 2], p: V| -> usize {
        let key = ((q[0] / UV_TOL / 10.0).floor() as i64, (q[1] / UV_TOL / 10.0).floor() as i64);
        for dx in -1..=1 { for dy in -1..=1 { for &i in grid.get(&(key.0 + dx, key.1 + dy)).into_iter().flatten() { if (uv[i][0] - q[0]).abs() < UV_TOL && (uv[i][1] - q[1]).abs() < UV_TOL { return i; } } } }
        uv.push(q); p3.push(p); grid.entry(key).or_default().push(uv.len() - 1); uv.len() - 1
    };
    let mut edges: Vec<(usize, usize)> = vec![];
    for (a, b) in bsegs.iter().map(|(a, b)| ((a.uv, a.p), (b.uv, b.p))).chain(csegs.iter().copied()) { let (i, j) = (vid(a.0, a.1), vid(b.0, b.1)); if i != j { edges.push((i, j)); } }
    edges.sort_unstable(); edges.dedup();
    // half-edges: h = 2e (a->b), 2e+1 (b->a); faces traced with the region on the left
    let he = |h: usize| if h % 2 == 0 { edges[h / 2] } else { (edges[h / 2].1, edges[h / 2].0) };
    let nh = edges.len() * 2;
    let mut out: Vec<Vec<usize>> = vec![vec![]; uv.len()];
    for h in 0..nh { out[he(h).0].push(h); }
    for (v, l) in out.iter_mut().enumerate() { let o = uv[v]; l.sort_by(|&a, &b| { let (pa, pb) = (uv[he(a).1], uv[he(b).1]); (pa[1] - o[1]).atan2(pa[0] - o[0]).total_cmp(&(pb[1] - o[1]).atan2(pb[0] - o[0])) }); }
    let mut pos = vec![0usize; nh];
    for l in &out { for (k, &h) in l.iter().enumerate() { pos[h] = k; } }
    let next = |h: usize| { let l = &out[he(h).1]; l[(pos[h ^ 1] + l.len() - 1) % l.len()] };
    let mut seen = vec![false; nh];
    let (mut regions, mut holes) = (vec![], vec![]);
    for h0 in 0..nh {
        if seen[h0] { continue; }
        let (mut cyc, mut h) = (vec![], h0);
        while !seen[h] { seen[h] = true; cyc.push(he(h).0); h = next(h); }
        let a = area(&uv, &cyc);
        if a > 1e-14 { regions.push(Region { outer: cyc, holes: vec![] }); }
        else if a < -1e-14 {
            // a CW cycle other than the outside of the face boundary is a hole
            let probe = uv[cyc[0]];
            let outside = cyc.len() >= n && f.lp.iter().any(|b| (b.uv[0] - probe[0]).abs() < UV_TOL && (b.uv[1] - probe[1]).abs() < UV_TOL);
            if !outside { holes.push(cyc); }
        }
    }
    for hc in holes {
        let probe = uv[hc[0]];
        let owner = regions.iter().enumerate().filter(|(_, r)| !r.outer.contains(&hc[0]) && inside_poly(probe, &r.outer.iter().map(|&i| uv[i]).collect::<Vec<_>>()))
            .min_by(|a, b| area(&uv, &a.1.outer).total_cmp(&area(&uv, &b.1.outer))).map(|x| x.0);
        if let Some(ri) = owner { regions[ri].holes.push(hc); }
    }
    Arr { uv, p: p3, regions }
}

/// a point strictly inside the region: offset left of a long boundary edge, verified
fn interior_point(arr: &Arr, r: &Region) -> Option<[f64; 2]> {
    let poly: Vec<[f64; 2]> = r.outer.iter().map(|&i| arr.uv[i]).collect();
    let holes: Vec<Vec<[f64; 2]>> = r.holes.iter().map(|h| h.iter().map(|&i| arr.uv[i]).collect()).collect();
    let n = poly.len();
    let mut cand: Vec<(f64, usize)> = (0..n).map(|i| { let (a, b) = (poly[i], poly[(i + 1) % n]); (((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt(), i) }).filter(|c| c.0 > 1e-12).collect();
    cand.sort_by(|x, y| y.0.total_cmp(&x.0));
    for &(l, i) in cand.iter().take(12) {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        for eps in [1e-3, 1e-4, 1e-5] {
            let d = eps * l.min(1.0) / l;
            let m = [(a[0] + b[0]) / 2.0 - (b[1] - a[1]) * d, (a[1] + b[1]) / 2.0 + (b[0] - a[0]) * d];
            if inside_poly(m, &poly) && !holes.iter().any(|h| inside_poly(m, h)) { return Some(m); }
        }
    }
    None
}

/// ear clipping of a polygon with holes (holes bridged to the nearest visible vertex); O(n^2)
fn triangulate(pts: &[[f64; 2]], outer: &[usize], holes: &[Vec<usize>]) -> Vec<[usize; 3]> {
    let mut poly = outer.to_vec();
    let mut hs = holes.to_vec();
    let maxx = |h: &Vec<usize>| h.iter().map(|&i| pts[i][0]).fold(f64::MIN, f64::max);
    hs.sort_by(|a, b| maxx(b).total_cmp(&maxx(a)));
    for h in hs {
        let hi = (0..h.len()).max_by(|&a, &b| pts[h[a]][0].total_cmp(&pts[h[b]][0])).unwrap();
        let hp = pts[h[hi]];
        let clear = |q: [f64; 2], ring: &[usize], skip: usize| (0..ring.len()).all(|k| { let (a, b) = (ring[k], ring[(k + 1) % ring.len()]); a == skip || b == skip || seg_cross(hp, q, pts[a], pts[b]).is_none() });
        let mut order: Vec<usize> = (0..poly.len()).collect();
        let d2 = |j: usize| { let q = pts[poly[j]]; (q[0] - hp[0]).powi(2) + (q[1] - hp[1]).powi(2) };
        order.sort_by(|&a, &b| d2(a).total_cmp(&d2(b)));
        let j = order.iter().copied().find(|&j| clear(pts[poly[j]], &poly, poly[j]) && clear(pts[poly[j]], &h, h[hi])).unwrap_or(order[0]);
        let mut np = poly[..=j].to_vec();
        np.extend((0..=h.len()).map(|k| h[(hi + k) % h.len()]));
        np.extend_from_slice(&poly[j..]);
        poly = np;
    }
    let mut tris = vec![];
    let (mut v, mut i, mut misses) = (poly, 0usize, 0usize);
    while v.len() > 3 {
        let n = v.len();
        i %= n;
        let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
        let turn = o2(pts[a], pts[b], pts[c]);
        if turn == 0 { v.remove(i); misses = 0; i = i.saturating_sub(1); continue; } // zero-area vertex
        let ear = turn > 0 && v.iter().all(|&q| q == a || q == b || q == c || pts[q] == pts[a] || pts[q] == pts[b] || pts[q] == pts[c]
            || !(o2(pts[a], pts[b], pts[q]) >= 0 && o2(pts[b], pts[c], pts[q]) >= 0 && o2(pts[c], pts[a], pts[q]) >= 0));
        if ear { tris.push([a, b, c]); v.remove(i); misses = 0; i = i.saturating_sub(1); } else { i += 1; misses += 1; if misses > n { break; } }
    }
    if v.len() == 3 { tris.push([v[0], v[1], v[2]]); }
    tris
}
/// signed volume contribution of a UV triangle: (1/3) ∫∫ p · (Su x Sv) du dv (divergence theorem), by a
/// degree-5 7-point Gauss rule after splitting until UV edges < 0.25 (the integrand is smooth)
fn integrate(s: &dyn Surf, a: [f64; 2], b: [f64; 2], c: [f64; 2], depth: u32) -> f64 {
    let l = |p: [f64; 2], q: [f64; 2]| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
    if depth < 16 && l(a, b).max(l(b, c)).max(l(c, a)) > 0.25 {
        let m = |p: [f64; 2], q: [f64; 2]| [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0];
        let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
        return integrate(s, a, ab, ca, depth + 1) + integrate(s, ab, b, bc, depth + 1) + integrate(s, ca, bc, c, depth + 1) + integrate(s, ab, bc, ca, depth + 1);
    }
    // Dunavant degree-5 rule: (barycentric a, b, weight)
    const R: [(f64, f64, f64); 7] = [
        (1.0 / 3.0, 1.0 / 3.0, 0.225),
        (0.059715871789770, 0.470142064105115, 0.132394152788506), (0.470142064105115, 0.059715871789770, 0.132394152788506), (0.470142064105115, 0.470142064105115, 0.132394152788506),
        (0.797426985353087, 0.101286507323456, 0.125939180544827), (0.101286507323456, 0.797426985353087, 0.125939180544827), (0.101286507323456, 0.101286507323456, 0.125939180544827),
    ];
    let jac = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.0; // signed UV area
    R.iter().map(|&(l1, l2, w)| {
        let l3 = 1.0 - l1 - l2;
        let d = s.ev(l1 * a[0] + l2 * b[0] + l3 * c[0], l1 * a[1] + l2 * b[1] + l3 * c[1]);
        w * dot(d.p, cross(d.u, d.v)) / 3.0
    }).sum::<f64>() * jac
}

// ------------------------------------------------------------------ boolean
#[derive(Clone, Copy, PartialEq, Debug)]
enum Op { Union, Inter, Diff }
struct Report { volume: f64, kept: usize, unmatched: usize, nonmanifold: usize, bad: usize }

fn boolean(a: &Solid, b: &Solid, op: Op) -> Report {
    VERTS.with(|v| v.borrow_mut().clear());
    let solids = [a, b];
    // shared edges: one SSI per face pair, given identically to both faces
    let mut cuts: [Vec<Vec<Vec<([f64; 2], V)>>>; 2] = [a.faces.iter().map(|_| vec![]).collect(), b.faces.iter().map(|_| vec![]).collect()];
    for (ia, fa) in a.faces.iter().enumerate() {
        for (ib, fb) in b.faces.iter().enumerate() {
            for pc in face_pair_edges(fa, a, fb, b) {
                cuts[0][ia].extend(cuts_for(fa, &pc.p, &pc.seam[0]));
                cuts[1][ib].extend(cuts_for(fb, &pc.p, &pc.seam[1]));
            }
        }
    }
    let mut r = Report { volume: 0.0, kept: 0, unmatched: 0, nonmanifold: 0, bad: 0 };
    let mut segs: Vec<(V, V)> = vec![];
    for si in 0..2 {
        let other = solids[1 - si];
        for (fi, f) in solids[si].faces.iter().enumerate() {
            let arr = arrange(f, std::mem::take(&mut cuts[si][fi]));
            for reg in &arr.regions {
                let Some(m) = interior_point(&arr, reg) else { r.bad += 1; continue };
                let d = (other.inside)(f.s.ev(m[0], m[1]).p);
                if d.abs() < 1e-9 { r.bad += 1; continue; }
                let (keep, flip) = match (op, si, d < 0.0) {
                    (Op::Union, _, false) | (Op::Inter, _, true) | (Op::Diff, 0, false) => (true, false),
                    (Op::Diff, 1, true) => (true, true),
                    _ => (false, false),
                };
                if !keep { continue; }
                r.kept += 1;
                for cyc in std::iter::once(&reg.outer).chain(&reg.holes) {
                    for i in 0..cyc.len() { let (x, y) = (arr.p[cyc[i]], arr.p[cyc[(i + 1) % cyc.len()]]); segs.push(if flip { (y, x) } else { (x, y) }); }
                }
                let tris = triangulate(&arr.uv, &reg.outer, &reg.holes);
                let want = area(&arr.uv, &reg.outer) + reg.holes.iter().map(|h| area(&arr.uv, h)).sum::<f64>();
                let got: f64 = tris.iter().map(|t| area(&arr.uv, t)).sum();
                if (got - want).abs() > 1e-9 * want.abs().max(1.0) { r.bad += 1; }
                let sgn = if flip { -1.0 } else { 1.0 };
                r.volume += sgn * tris.iter().map(|t| integrate(f.s.as_ref(), arr.uv[t[0]], arr.uv[t[1]], arr.uv[t[2]], 0)).sum::<f64>();
            }
        }
    }
    // closure: every directed boundary segment matched by its reverse (vertices welded at 1e-7)
    let key = |p: V| p.map(|x| (x * 1e7).round() as i64);
    let mut bal: HashMap<([i64; 3], [i64; 3]), i64> = HashMap::new();
    let mut und: HashMap<([i64; 3], [i64; 3]), usize> = HashMap::new();
    for (x, y) in &segs {
        let (kx, ky) = (key(*x), key(*y));
        if kx == ky { continue; }
        *bal.entry((kx, ky)).or_default() += 1;
        *bal.entry((ky, kx)).or_default() -= 1;
        *und.entry(if kx < ky { (kx, ky) } else { (ky, kx) }).or_default() += 1;
    }
    r.unmatched = bal.values().filter(|v| **v > 0).map(|v| *v as usize).sum();
    r.nonmanifold = und.values().filter(|v| **v > 2).count();
    r
}

pub fn run() -> usize {
    println!("== boolean spike (closure = unmatched directed boundary segments) ==");
    let (x, y, z, o) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0; 3]);
    let cases: Vec<(&str, Solid, Solid, Op, f64)> = vec![
        ("S1 cube alone (U far cube)", cube(o, 1.0), cube([10.0, 0.0, 0.0], 1.0), Op::Union, 16.0),
        ("S2 cylinder alone (U far cube)", cylinder(o, 0.5, 2.0, x, y, z), cube([10.0, 0.0, 0.0], 1.0), Op::Union, 8.0 + PI),
        ("B1 cube - cylinder r0.5 (through hole)", cube(o, 1.0), cylinder(o, 0.5, 2.0, x, y, z), Op::Diff, 8.0 - PI * 0.5),
        ("B2 cube U cube offset (0.5,0.5,0.5)", cube(o, 1.0), cube([0.5, 0.5, 0.5], 1.0), Op::Union, 16.0 - 1.5f64.powi(3)),
        ("B3 cyl ^ cyl perpendicular (Steinmetz)", cylinder(o, 1.0, 2.0, x, y, z), cylinder(o, 1.0, 2.0, y, z, x), Op::Inter, 16.0 / 3.0),
        ("B4 cube - cylinder r1 (tangent to 4 faces)", cube(o, 1.0), cylinder(o, 1.0, 2.0, x, y, z), Op::Diff, 8.0 - 2.0 * PI),
        ("B5 cube ^ cylinder on face plane x=1", cube(o, 1.0), cylinder([1.0, 0.0, 0.0], 0.8, 2.0, x, y, z), Op::Inter, PI * 0.64),
    ];
    let mut fails = 0;
    for (name, a, b, op, v) in cases {
        let t0 = Instant::now();
        let r = boolean(&a, &b, op);
        let rel = (r.volume - v).abs() / v;
        let ok = r.unmatched == 0 && rel < 1e-3 && r.bad == 0;
        fails += !ok as usize;
        println!("  {name:<44} faces {:3} | volume {:.6} vs {v:.6} (rel {rel:.1e}) | unmatched {} | non-manifold {} | bad {} | {:6.0} ms {}", r.kept, r.volume, r.unmatched, r.nonmanifold, r.bad, ms(t0), if ok { "PASS" } else { "FAIL" });
    }
    fails
}
