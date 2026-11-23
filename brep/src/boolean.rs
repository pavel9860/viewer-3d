//! Boolean operations on bodies: union, intersection, difference.
//!
//! 1. Both bodies are merged into one working body (B's ids offset).
//! 2. Per face pair, in parallel: SSI curves restricted to both faces, with exact points (edge curve ∩ other
//!    surface) where they cross or touch either face's edges.
//! 3. Vertices are welded globally; original edges are split at those points; curve pieces become edges shared
//!    by the two faces (or are identified with an existing edge they run along).
//! 4. Per face, in parallel: arrangement of its split boundary and its cuts in UV (half-edges) -> regions;
//!    each region is classified against the other body and kept by the operation.
//! 5. Kept regions become the faces of the result (B's reversed for a difference).
use crate::classify::{Side, Solid};
use crate::geom::*;
use crate::nurbs::{t_of, tess, Crv, PolyC};
use crate::par;
use crate::ssi::{intersect, solve_n, End, Kind, SAG as SSI_SAG};
use crate::topo::*;
use crate::trim::{area, in_poly, seg_cross, seg_dist, Trim, P2, SAG};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::Arc;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Op { Union, Intersection, Difference }
/// provenance op id of boolean results; `src` = 1 + working-body id of the input entity (0 = new)
pub const OP_BOOLEAN: u32 = 1;

/// UV distance below which a curve segment counts as touching a boundary segment
const TOUCH: f64 = 1e-4;
/// residual accepted for an edge/surface point
const ACC: f64 = 1e-9;
/// UV distance for points on a face boundary, and for arrangement nodes of the same vertex
const UV_ON: f64 = 1e-9;
const NODE: f64 = 1e-6;

#[derive(Default, Debug, Clone)]
pub struct Stats {
    pub pairs: usize,
    pub pieces: usize,
    pub cut_edges: usize,
    pub regions: usize,
    pub kept: usize,
    /// cut ends that meet nothing in the arrangement (must be 0)
    pub dangling: usize,
    /// regions on the other body's boundary (coincident faces; not handled yet)
    pub on: usize,
}

#[derive(Clone)]
struct Ev { e: Id<Edge>, t: f64, p: V }
/// part of an intersection curve inside both faces; `ev[i]` = edges point i lies on
struct Piece { f: [Id<Face>; 2], pts: Vec<V>, ev: Vec<Vec<Ev>>, closed: bool }

/// append `o` to `w`, returning the face id offset
fn append(w: &mut Body, o: &Body) -> u32 {
    let (nv, ne, nc, nl, nf) = (w.vertices.len() as u32, w.edges.len() as u32, w.coedges.len() as u32, w.loops.len() as u32, w.faces.len() as u32);
    for (_, v) in o.vertices.iter() { w.vertices.push(v.clone()); }
    for (_, e) in o.edges.iter() { let mut e = e.clone(); e.v = e.v.map(|v| Id::new((v.0 + nv) as usize)); w.edges.push(e); }
    for (_, c) in o.coedges.iter() { w.coedges.push(Coedge { edge: Id::new((c.edge.0 + ne) as usize), sense: c.sense }); }
    for (_, l) in o.loops.iter() { w.loops.push(Loop { coedges: l.coedges.iter().map(|c| Id::new((c.0 + nc) as usize)).collect() }); }
    for (_, f) in o.faces.iter() { let mut f = f.clone(); f.loops = f.loops.iter().map(|l| Id::new((l.0 + nl) as usize)).collect(); w.faces.push(f); }
    nf
}

/// point where curve c meets surface s, from a guess (damped Gauss-Newton on c(t) - s(u, v)): (t, residual)
fn curve_surface(c: &dyn Crv, s: &dyn Surf, t0: f64, p0: V) -> (f64, f64) {
    let uv = s.uv_of(p0);
    let mut x = [t0, uv[0], uv[1]];
    let res = |x: &[f64; 3]| sub(c.ev(x[0]).0, s.ev(x[1], x[2]).p);
    let mut r = res(&x);
    let mut lam = 1e-12;
    for _ in 0..80 {
        let f = dot(r, r);
        if f < 1e-30 { break; }
        let (d1, e) = (c.ev(x[0]).1, s.ev(x[1], x[2]));
        let j = [d1, mulv(e.u, -1.0), mulv(e.v, -1.0)];
        let mut a = [[0.0; 3]; 3];
        for p in 0..3 { for q in 0..3 { a[p][q] = dot(j[p], j[q]); } }
        let g = [0, 1, 2].map(|p| -dot(j[p], r));
        let tr = (a[0][0] + a[1][1] + a[2][2]).max(1e-300);
        let mut ok = false;
        for _ in 0..12 {
            let mut ad = a;
            for k in 0..3 { ad[k][k] += lam * tr; }
            if let Some(dx) = solve_n(ad, g) {
                let xn = [x[0] + dx[0], x[1] + dx[1], x[2] + dx[2]];
                let rn = res(&xn);
                if dot(rn, rn) < f { x = xn; r = rn; lam = (lam * 0.1).max(1e-15); ok = true; break; }
            }
            lam *= 10.0;
        }
        if !ok { break; }
    }
    (x[0], norm(r))
}

/// t moved by whole periods into the edge's range (within a relative slack), if possible
fn on_edge(e: &Edge, c: &dyn Crv, t: f64) -> Option<f64> {
    let span = e.t[1] - e.t[0];
    let ks: &[f64] = if c.periodic() { &[0.0, 1.0, -1.0, 2.0, -2.0] } else { &[0.0] };
    ks.iter().map(|k| t + k * 2.0 * PI).find(|t| { let l = (t - e.t[0]) / span; (-1e-7..=1.0 + 1e-7).contains(&l) })
}

/// exact points where the polyline `pts` (on the face's surface and on `other`) crosses or touches the face's edges
/// event points within `weld` of a polyline point are that point
fn boundary_events(w: &Body, tr: &Trim, s: &dyn Surf, other: &dyn Surf, weld: f64, pts: &mut Vec<V>, ev: &mut Vec<Vec<Ev>>) {
    let mut k = 0;
    while k + 1 < pts.len() {
        let (p0, p1) = (pts[k], pts[k + 1]);
        let a = tr.uv(s, p0);
        let b = tr.unwrap(Some(a), tr.uv(s, p1));
        let mut hits: Vec<(f64, Ev)> = vec![];
        for sh in tr.shifts() {
            let (a, b) = ([a[0] + sh[0], a[1] + sh[1]], [b[0] + sh[0], b[1] + sh[1]]);
            let bb = [a[0].min(b[0]) - TOUCH, a[0].max(b[0]) + TOUCH, a[1].min(b[1]) - TOUCH, a[1].max(b[1]) + TOUCH];
            for it in tr.loops.iter().flatten() {
                let ed = w.edges.get(it.e);
                let Some(c) = ed.curve.as_deref() else { continue };
                for j in 0..it.uv.len() - 1 {
                    let (c0, c1) = (it.uv[j], it.uv[j + 1]);
                    if c0[0].max(c1[0]) < bb[0] || c0[0].min(c1[0]) > bb[1] || c0[1].max(c1[1]) < bb[2] || c0[1].min(c1[1]) > bb[3] { continue; }
                    let (sa, sb) = match seg_cross(a, b, c0, c1) {
                        Some(x) => x,
                        None if seg_dist(a, b, c0, c1) <= TOUCH => (0.5, 0.5),
                        None => continue,
                    };
                    let g = addv(p0, mulv(sub(p1, p0), sa));
                    let tg = it.t[j] + (it.t[j + 1] - it.t[j]) * sb;
                    let (t, res) = curve_surface(c, other, tg, g);
                    if !(res <= ACC) { continue; }
                    let Some(t) = on_edge(ed, c, t) else { continue };
                    // an edge running on the other surface overlaps it: no isolated point (its ends come from other edges)
                    let h = 1e-3 * (ed.t[1] - ed.t[0]);
                    if [t - h, t + h].iter().all(|&x| dist(other, c.ev(x).0) <= ACC) { continue; }
                    let q = c.ev(t).0;
                    let ab = sub(p1, p0);
                    let along = (dot(sub(q, p0), ab) / dot(ab, ab)).clamp(0.0, 1.0);
                    if norm(sub(q, addv(p0, mulv(ab, along)))) > 10.0 * SSI_SAG { continue; }
                    if !hits.iter().any(|h| h.1.e == it.e && norm(sub(h.1.p, q)) < weld) { hits.push((along, Ev { e: it.e, t, p: q })); }
                }
            }
        }
        hits.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut ins: Vec<(V, Vec<Ev>)> = vec![];
        for (_, h) in hits {
            let at = if norm(sub(h.p, p0)) < weld { Some(k) } else if norm(sub(h.p, p1)) < weld { Some(k + 1) } else { None };
            match at {
                Some(i) => if !ev[i].iter().any(|x| x.e == h.e) { ev[i].push(h); },
                None => match ins.iter_mut().find(|x| norm(sub(x.0, h.p)) < weld) { Some(x) => x.1.push(h), None => ins.push((h.p, vec![h])) },
            }
        }
        let n = ins.len();
        for (j, (p, e)) in ins.into_iter().enumerate() { pts.insert(k + 1 + j, p); ev.insert(k + 1 + j, e); }
        k += 1 + n;
    }
}

/// SSI of two faces restricted to both, with boundary events
fn pair_pieces(w: &Body, trims: &[Trim], weld: f64, fa: Id<Face>, fb: Id<Face>) -> Vec<Piece> {
    let (sa, sb) = (w.faces.get(fa).surface.as_ref(), w.faces.get(fb).surface.as_ref());
    // only the parts under the faces matter: intersect the surfaces over the faces' parameter boxes
    let (wa, wb) = (window(sa, &trims[fa.ix()]), window(sb, &trims[fb.ix()]));
    let t0 = std::time::Instant::now();
    let curves = intersect(&wa, &wb).0;
    if std::env::var("BREP_TIME").is_ok() && t0.elapsed().as_secs_f64() > 0.5 { eprintln!("  ssi {:.1} s, {} curves", t0.elapsed().as_secs_f64(), curves.len()); }
    let out = restrict(w, trims, weld, fa, fb, curves);
    if std::env::var("BREP_WIN").is_ok() {
        let full = restrict(w, trims, weld, fa, fb, intersect(sa, sb).0);
        let len = |ps: &[Piece]| ps.iter().map(|p| p.pts.windows(2).map(|q| norm(sub(q[1], q[0]))).sum::<f64>()).sum::<f64>();
        if (len(&full) - len(&out)).abs() > 1e-6 { eprintln!("WIN faces src {} x {}: window {:?} pieces {} len {:.5} | full pieces {} len {:.5}", w.faces.get(fa).origin.src, w.faces.get(fb).origin.src, wa.dom(), out.len(), len(&out), full.len(), len(&full)); }
    }
    out
}

/// intersection curves of two faces' surfaces restricted to both faces, with boundary events
fn restrict(w: &Body, trims: &[Trim], weld: f64, fa: Id<Face>, fb: Id<Face>, curves: Vec<crate::ssi::Curve>) -> Vec<Piece> {
    let (sa, sb) = (w.faces.get(fa).surface.as_ref(), w.faces.get(fb).surface.as_ref());
    let mut out = vec![];
    for c in curves {
        if c.kind == Kind::Isolated || c.pts.len() < 2 { continue; }
        let mut pts = c.pts;
        let mut ev = vec![vec![]; pts.len()];
        for (f, s, other) in [(fa, sa, sb), (fb, sb, sa)] { boundary_events(w, &trims[f.ix()], s, other, weld, &mut pts, &mut ev); }
        let inside = |p: V| [(fa, sa), (fb, sb)].iter().all(|&(f, s)| { let t = &trims[f.ix()]; t.covers(t.uv(s, p), UV_ON) });
        let mut seg_in: Vec<bool> = pts.windows(2).map(|w| inside(mulv(addv(w[0], w[1]), 0.5))).collect();
        let closed = c.ends.0 == End::Closed && norm(sub(pts[0], pts[pts.len() - 1])) < 1e-12;
        if closed {
            // cyclic: start at an outside segment, else at an event, else keep as a closed piece
            let n = pts.len() - 1;
            let start = seg_in.iter().position(|x| !x).or_else(|| (0..n).find(|&i| !ev[i].is_empty()));
            match start {
                Some(s) if s > 0 => {
                    let mut e0 = std::mem::take(&mut ev[0]);
                    e0.extend(ev[n].drain(..));
                    ev[0] = e0;
                    pts.pop(); ev.pop();
                    pts.rotate_left(s); ev.rotate_left(s); seg_in.rotate_left(s);
                    pts.push(pts[0]); ev.push(ev[0].clone());
                }
                Some(_) => {}
                None => { out.push(Piece { f: [fa, fb], pts, ev, closed: true }); continue; }
            }
        }
        let mut k = 0;
        while k < seg_in.len() {
            if !seg_in[k] { k += 1; continue; }
            let k0 = k;
            while k < seg_in.len() && seg_in[k] { k += 1; }
            out.push(Piece { f: [fa, fb], pts: pts[k0..=k].to_vec(), ev: ev[k0..=k].to_vec(), closed: false });
        }
    }
    out
}

/// the surface over the face's parameter box (raw parameters), padded by 2%
fn window<'a>(s: &'a dyn Surf, tr: &Trim) -> Window<'a> {
    let (u0, u1) = { let a = tr.raw([tr.bb[0], 0.0])[0]; let b = tr.raw([tr.bb[1], 0.0])[0]; (a.min(b), a.max(b)) };
    let (pu, pv) = (0.02 * (u1 - u0) + 1e-9, 0.02 * (tr.bb[3] - tr.bb[2]) + 1e-9);
    Window::new(s, [u0 - pu, u1 + pu, tr.bb[2] - pv, tr.bb[3] + pv])
}

/// vertices welded within a distance (uniform grid)
struct Weld { tol: f64, grid: HashMap<[i64; 3], Vec<Id<Vertex>>> }
impl Weld {
    fn key(&self, p: V) -> [i64; 3] { p.map(|x| (x / (4.0 * self.tol)).floor() as i64) }
    fn add(&mut self, id: Id<Vertex>, p: V) { let k = self.key(p); self.grid.entry(k).or_default().push(id); }
    fn get(&mut self, w: &mut Body, p: V) -> Id<Vertex> {
        let k = self.key(p);
        let mut best: Option<(f64, Id<Vertex>)> = None;
        for dx in -1..=1 { for dy in -1..=1 { for dz in -1..=1 {
            for &v in self.grid.get(&[k[0] + dx, k[1] + dy, k[2] + dz]).into_iter().flatten() {
                let d = norm(sub(w.vertices.get(v).p, p));
                if d <= self.tol && best.map_or(true, |b| d < b.0) { best = Some((d, v)); }
            }
        } } }
        if let Some((d, v)) = best { let vt = &mut w.vertices.get_mut(v).tol; *vt = vt.max(d); return v; }
        let v = w.vertices.push(Vertex { p, tol: w.tol });
        self.add(v, p);
        v
    }
}

/// a region of a face: loops of (edge, coedge sense), outer first
struct Region { loops: Vec<Vec<(Id<Edge>, bool)>>, probe: P2 }

/// arrangement of a face's split boundary and cuts in UV: regions inside the face, and the dangling count
fn arrange(w: &Body, s: &dyn Surf, tr: &Trim, cuts: &[Id<Edge>]) -> (Vec<Region>, usize) {
    // items: (edge, polyline along the edge use, use sense, end vertices)
    let mut items: Vec<(Id<Edge>, Vec<P2>, bool, [Id<Vertex>; 2])> = tr.loops.iter().flatten().map(|it| (it.e, it.uv.clone(), it.sense, it.v)).collect();
    for &e in cuts {
        let ed = w.edges.get(e);
        let c = ed.curve.as_deref().unwrap();
        let mut p: Vec<V> = tess(c, ed.t[0], ed.t[1], SAG).into_iter().map(|t| c.ev(t).0).collect();
        let n = p.len();
        p[0] = w.vertices.get(ed.v[0]).p;
        p[n - 1] = w.vertices.get(ed.v[1]).p;
        items.push((e, tr.polyline(s, &p), true, ed.v));
    }
    let mut nodes: Vec<(Id<Vertex>, P2)> = vec![];
    let mut node = |v: Id<Vertex>, q: P2| -> usize {
        if let Some(i) = nodes.iter().position(|n| n.0 == v && (n.1[0] - q[0]).hypot(n.1[1] - q[1]) < NODE) { return i; }
        nodes.push((v, q));
        nodes.len() - 1
    };
    let ends: Vec<[usize; 2]> = items.iter().map(|it| [node(it.3[0], it.1[0]), node(it.3[1], it.1[it.1.len() - 1])]).collect();
    // half-edge h: item h/2, forward when even; region on the left
    let nh = 2 * items.len();
    let from = |h: usize| ends[h / 2][h % 2];
    let to = |h: usize| ends[h / 2][1 - h % 2];
    let path = |h: usize| -> Vec<P2> { let p = &items[h / 2].1; if h % 2 == 0 { p.clone() } else { p.iter().rev().copied().collect() } };
    let angle = |h: usize| { let p = path(h); let o = p[0]; let q = p.iter().find(|q| **q != o).copied().unwrap_or(o); (q[1] - o[1]).atan2(q[0] - o[0]) };
    let mut out: Vec<Vec<usize>> = vec![vec![]; nodes.len()];
    for h in 0..nh { out[from(h)].push(h); }
    let ang: Vec<f64> = (0..nh).map(angle).collect();
    for l in out.iter_mut() { l.sort_by(|&a, &b| ang[a].total_cmp(&ang[b])); }
    let dangling = out.iter().filter(|l| l.len() == 1).count();
    let mut pos = vec![0usize; nh];
    for l in &out { for (k, &h) in l.iter().enumerate() { pos[h] = k; } }
    let next = |h: usize| { let l = &out[to(h)]; l[(pos[h ^ 1] + l.len() - 1) % l.len()] };
    let mut seen = vec![false; nh];
    let mut outer: Vec<(Vec<usize>, Vec<P2>)> = vec![];
    let mut holes: Vec<(Vec<usize>, Vec<P2>)> = vec![];
    for h0 in 0..nh {
        if seen[h0] { continue; }
        let (mut cyc, mut h) = (vec![], h0);
        while !seen[h] { seen[h] = true; cyc.push(h); h = next(h); }
        let poly: Vec<P2> = cyc.iter().flat_map(|&h| { let mut p = path(h); p.pop(); p }).collect();
        let a = area(&poly);
        if a > 0.0 { outer.push((cyc, poly)); } else if a < 0.0 { holes.push((cyc, poly)); }
    }
    // a point left of the longest segments of `poly`, strictly inside it, accepted by `ok`
    let probe = |poly: &[P2], ok: &dyn Fn(P2) -> bool| -> Option<P2> {
        let n = poly.len();
        let size = area(poly).abs().sqrt();
        let mut segs: Vec<(f64, usize)> = (0..n).map(|i| { let (a, b) = (poly[i], poly[(i + 1) % n]); ((b[0] - a[0]).hypot(b[1] - a[1]), i) }).filter(|x| x.0 > 0.0).collect();
        segs.sort_by(|x, y| y.0.total_cmp(&x.0));
        for f in [0.1, 0.03, 0.01, 1e-3, 1e-4, 1e-5] {
            for &(l, i) in segs.iter().take(8) {
                let (a, b) = (poly[i], poly[(i + 1) % n]);
                let d = f * size / l;
                let m = [(a[0] + b[0]) / 2.0 - (b[1] - a[1]) * d, (a[1] + b[1]) / 2.0 + (b[0] - a[0]) * d];
                if in_poly(m, poly) && ok(m) { return Some(m); }
            }
        }
        None
    };
    // regions: counter-clockwise cycles bounding face material (not the inside of a face hole)
    let mut regs: Vec<(Vec<usize>, Vec<P2>, Vec<usize>)> = outer.into_iter().filter(|(_, p)| probe(p, &|m| tr.inside(m)).is_some()).map(|(c, p)| (c, p, vec![])).collect();
    // holes: clockwise cycles, owned by the smallest region containing a point just left of them
    for (hi, (_, hp)) in holes.iter().enumerate() {
        let n = hp.len();
        let Some(i) = (0..n).max_by(|&i, &j| { let l = |i: usize| (hp[(i + 1) % n][0] - hp[i][0]).hypot(hp[(i + 1) % n][1] - hp[i][1]); l(i).total_cmp(&l(j)) }) else { continue };
        let (a, b) = (hp[i], hp[(i + 1) % n]);
        let d = 1e-7 / (b[0] - a[0]).hypot(b[1] - a[1]).max(1e-300);
        let m = [(a[0] + b[0]) / 2.0 - (b[1] - a[1]) * d, (a[1] + b[1]) / 2.0 + (b[0] - a[0]) * d];
        if let Some(r) = (0..regs.len()).filter(|&r| in_poly(m, &regs[r].1)).min_by(|&x, &y| area(&regs[x].1).total_cmp(&area(&regs[y].1))) { regs[r].2.push(hi); }
    }
    let lp = |cyc: &[usize]| cyc.iter().map(|&h| (items[h / 2].0, items[h / 2].2 == (h % 2 == 0))).collect::<Vec<_>>();
    let regions = regs.iter().filter_map(|(c, p, hs)| {
        let hp: Vec<&Vec<P2>> = hs.iter().map(|&i| &holes[i].1).collect();
        let m = probe(p, &|m| !hp.iter().any(|h| in_poly(m, h)))?;
        Some(Region { loops: std::iter::once(lp(c)).chain(hs.iter().map(|&i| lp(&holes[i].0))).collect(), probe: m })
    }).collect();
    (regions, dangling)
}

/// 3D box of a face: boundary points and a surface grid over its UV box, padded
fn face_box(w: &Body, f: Id<Face>, tr: &Trim) -> [f64; 6] {
    let s = w.faces.get(f).surface.as_ref();
    let mut bb = [f64::MAX, f64::MIN, f64::MAX, f64::MIN, f64::MAX, f64::MIN];
    let mut add = |p: V| for k in 0..3 { bb[2 * k] = bb[2 * k].min(p[k]); bb[2 * k + 1] = bb[2 * k + 1].max(p[k]); };
    for it in tr.loops.iter().flatten() { for &p in &it.p { add(p); } }
    for i in 0..=8 { for j in 0..=8 {
        let q = [tr.bb[0] + (tr.bb[1] - tr.bb[0]) * i as f64 / 8.0, tr.bb[2] + (tr.bb[3] - tr.bb[2]) * j as f64 / 8.0];
        let r = tr.raw(q);
        add(s.ev(r[0], r[1]).p);
    } }
    let pad = 0.05 * (0..3).map(|k| bb[2 * k + 1] - bb[2 * k]).fold(0.0, f64::max) + 1e-6;
    for k in 0..3 { bb[2 * k] -= pad; bb[2 * k + 1] += pad; }
    bb
}
fn loops_of(w: &Body, f: Id<Face>, split: &HashMap<Id<Edge>, Vec<Id<Edge>>>) -> Vec<Vec<(Id<Edge>, bool)>> {
    w.faces.get(f).loops.iter().map(|&l| w.loops.get(l).coedges.iter().flat_map(|&c| {
        let ce = w.coedges.get(c);
        let subs = split.get(&ce.edge).cloned().unwrap_or_else(|| vec![ce.edge]);
        let v: Vec<(Id<Edge>, bool)> = if ce.sense { subs.into_iter().map(|e| (e, true)).collect() } else { subs.into_iter().rev().map(|e| (e, false)).collect() };
        v
    }).collect()).collect()
}

/// boolean of two closed bodies
pub fn boolean(a: &Body, b: &Body, op: Op, workers: usize) -> (Body, Stats) {
    let mut st = Stats::default();
    let clock = std::time::Instant::now();
    let lap = |what: &str| if std::env::var("BREP_TIME").is_ok() { eprintln!("  {what}: {:.0} ms", clock.elapsed().as_secs_f64() * 1e3); };
    let mut w = a.clone();
    w.tol = a.tol.max(b.tol);
    let boff = append(&mut w, b);
    let faces: Vec<Id<Face>> = w.faces.ids().collect();
    let (fa, fb): (Vec<Id<Face>>, Vec<Id<Face>>) = faces.iter().partition(|f| f.0 < boff);
    let mut loops0: Vec<Vec<Vec<(Id<Edge>, bool)>>> = faces.iter().map(|&f| loops_of(&w, f, &HashMap::new())).collect();
    let trims: Vec<Trim> = par::map(faces.len(), workers, |i| { let fc = w.faces.get(faces[i]); Trim::new(&w, fc.surface.as_ref(), fc.sense, &loops0[i]) });
    let boxes: Vec<[f64; 6]> = par::map(faces.len(), workers, |i| face_box(&w, faces[i], &trims[i]));
    let overlap = |p: &[f64; 6], q: &[f64; 6]| (0..3).all(|k| p[2 * k] <= q[2 * k + 1] && q[2 * k] <= p[2 * k + 1]);
    let pairs: Vec<(Id<Face>, Id<Face>)> = fa.iter().flat_map(|&x| fb.iter().map(move |&y| (x, y))).filter(|(x, y)| overlap(&boxes[x.ix()], &boxes[y.ix()])).collect();
    st.pairs = pairs.len();
    lap("trims, boxes, pairs");
    let weld_tol = (10.0 * w.tol).max(1e-7);
    let pieces: Vec<Piece> = par::map(pairs.len(), workers, |i| {
        let t0 = std::time::Instant::now();
        let r = pair_pieces(&w, &trims, weld_tol, pairs[i].0, pairs[i].1);
        let dt = t0.elapsed().as_secs_f64();
        if dt > 0.5 && std::env::var("BREP_TIME").is_ok() { eprintln!("  slow pair: faces src {} x {}: {dt:.1} s, {} pieces", w.faces.get(pairs[i].0).origin.src, w.faces.get(pairs[i].1).origin.src, r.len()); }
        r
    }).into_iter().flatten().collect();
    st.pieces = pieces.len();
    lap("face-pair intersections");

    // vertices: every event point welded (existing vertices first), piece ends without events become vertices
    let mut weld = Weld { tol: weld_tol, grid: HashMap::new() };
    for (id, v) in w.vertices.iter() { weld.add(id, v.p); }
    let mut splits: HashMap<Id<Edge>, Vec<(f64, Id<Vertex>)>> = HashMap::new();
    let mut pverts: Vec<Vec<Option<Id<Vertex>>>> = vec![];
    for pc in &pieces {
        let n = pc.pts.len();
        let mut pv = vec![None; n];
        for i in 0..n {
            if pc.closed && i == n - 1 { pv[i] = pv[0]; continue; }
            if let Some(e0) = pc.ev[i].first() {
                let v = weld.get(&mut w, e0.p);
                for e in &pc.ev[i] {
                    // the vertex tolerance covers every edge end that meets here
                    let d = norm(sub(e.p, w.vertices.get(v).p));
                    let vt = &mut w.vertices.get_mut(v).tol;
                    *vt = vt.max(d);
                    splits.entry(e.e).or_default().push((e.t, v));
                }
                pv[i] = Some(v);
            } else if (!pc.closed && (i == 0 || i == n - 1)) || (pc.closed && i == 0 && pc.ev.iter().all(|e| e.is_empty())) {
                pv[i] = Some(weld.get(&mut w, pc.pts[i]));
            }
        }
        if pc.closed { pv[n - 1] = pv[0]; }
        pverts.push(pv);
    }
    // split edges at their event vertices
    let mut split: HashMap<Id<Edge>, Vec<Id<Edge>>> = HashMap::new();
    let mut ekeys: Vec<Id<Edge>> = splits.keys().copied().collect();
    ekeys.sort_unstable_by_key(|e| e.0);
    for e in ekeys {
        let ed = w.edges.get(e).clone();
        let span = ed.t[1] - ed.t[0];
        let mut sp: Vec<(f64, Id<Vertex>)> = splits[&e].iter().filter(|(_, v)| !ed.v.contains(v)).copied().collect();
        sp.sort_by(|x, y| ((x.0 - ed.t[0]) / span).total_cmp(&((y.0 - ed.t[0]) / span)));
        sp.dedup_by_key(|x| x.1);
        if sp.is_empty() { continue; }
        let mut cuts = vec![(ed.t[0], ed.v[0])];
        cuts.extend(sp);
        cuts.push((ed.t[1], ed.v[1]));
        let subs = cuts.windows(2).map(|c| w.edges.push(Edge { curve: ed.curve.clone(), v: [c[0].1, c[1].1], t: [c[0].0, c[1].0], tol: ed.tol, origin: Origin { op: OP_BOOLEAN, src: e.0 + 1 } })).collect();
        split.insert(e, subs);
    }
    // live edges by end vertices, to identify curve pieces running along existing edges
    let mut by_ends: HashMap<(u32, u32), Vec<Id<Edge>>> = HashMap::new();
    let key = |a: Id<Vertex>, b: Id<Vertex>| (a.0.min(b.0), a.0.max(b.0));
    for (id, e) in w.edges.iter() { if !split.contains_key(&id) { by_ends.entry(key(e.v[0], e.v[1])).or_default().push(id); } }
    let mut cuts: Vec<Vec<Id<Edge>>> = vec![vec![]; faces.len()];
    for (pc, pv) in pieces.iter().zip(&pverts) {
        let idx: Vec<usize> = (0..pc.pts.len()).filter(|&i| pv[i].is_some()).collect();
        for s in idx.windows(2) {
            let (i, j) = (s[0], s[1]);
            let (va, vb) = (pv[i].unwrap(), pv[j].unwrap());
            let mut p = pc.pts[i..=j].to_vec();
            p[0] = w.vertices.get(va).p;
            *p.last_mut().unwrap() = w.vertices.get(vb).p;
            // a piece shorter than the weld tolerance that closes on one vertex is below model resolution
            if va == vb && p.windows(2).map(|q| norm(sub(q[1], q[0]))).sum::<f64>() <= weld_tol { continue; }
            let mid = p[p.len() / 2];
            let same = by_ends.get(&key(va, vb)).into_iter().flatten().copied().find(|&e| {
                let ed = w.edges.get(e);
                ed.curve.as_deref().map_or(false, |c| { let t = t_of(c, mid); on_edge(ed, c, t).is_some() && norm(sub(c.ev(t).0, mid)) < 10.0 * SSI_SAG })
            });
            let e = same.unwrap_or_else(|| {
                let (sa, sb) = (w.faces.get(pc.f[0]).surface.clone(), w.faces.get(pc.f[1]).surface.clone());
                let dev = p.windows(2).flat_map(|q| [q[0], mulv(addv(q[0], q[1]), 0.5)]).map(|q| dist(sa.as_ref(), q).max(dist(sb.as_ref(), q))).fold(0.0, f64::max);
                let c = PolyC::new(p);
                let t = [0.0, c.range().1];
                st.cut_edges += 1;
                let e = w.edges.push(Edge { curve: Some(Arc::new(c)), v: [va, vb], t, tol: (2.0 * dev).max(w.tol), origin: Origin { op: OP_BOOLEAN, src: 0 } });
                by_ends.entry(key(va, vb)).or_default().push(e);
                e
            });
            for f in pc.f { if !cuts[f.ix()].contains(&e) { cuts[f.ix()].push(e); } }
        }
    }

    lap("vertices, splits, cut edges");
    // per face: arrangement, classification against the other body, selection
    for (i, &f) in faces.iter().enumerate() { loops0[i] = loops_of(&w, f, &split); }
    let (sol_a, sol_b) = (Solid::new(&w, &fa), Solid::new(&w, &fb));
    let debug = std::env::var("BREP_DEBUG").is_ok();
    let res: Vec<(Vec<(Region, bool)>, usize, usize, usize)> = par::map(faces.len(), workers, |i| {
        let f = faces[i];
        let fc = w.faces.get(f);
        let s = fc.surface.as_ref();
        let tr = Trim::new(&w, s, fc.sense, &loops0[i]);
        let bnd: Vec<Id<Edge>> = loops0[i].iter().flatten().map(|x| x.0).collect();
        let cs: Vec<Id<Edge>> = cuts[i].iter().copied().filter(|e| !bnd.contains(e)).collect();
        let (regions, dangling) = arrange(&w, s, &tr, &cs);
        let in_b = f.0 >= boff;
        let other = if in_b { &sol_a } else { &sol_b };
        let (n, mut on) = (regions.len(), 0);
        if debug { eprintln!("face {i} (src {}) loops {} cuts {} regions {n} dangling {dangling} trim bb {:?}", fc.origin.src, loops0[i].len(), cs.len(), tr.bb); }
        let kept = regions.into_iter().filter_map(|r| {
            let q = tr.raw(r.probe);
            let side = other.side(s.ev(q[0], q[1]).p, 1e3 * ACC);
            if debug { eprintln!("   region loops {} probe {:?} p {:?} -> {side:?}", r.loops.len(), r.probe, s.ev(q[0], q[1]).p); }
            if side == Side::On { on += 1; return None; }
            let inside = side == Side::In;
            match (op, in_b, inside) {
                (Op::Union, _, false) | (Op::Intersection, _, true) | (Op::Difference, false, false) => Some((r, false)),
                (Op::Difference, true, true) => Some((r, true)),
                _ => None,
            }
        }).collect();
        (kept, n, dangling, on)
    });

    lap("arrangements + classification");
    // result body
    let mut r = Body { tol: w.tol, ..Default::default() };
    let mut vmap: HashMap<Id<Vertex>, Id<Vertex>> = HashMap::new();
    let mut emap: HashMap<Id<Edge>, Id<Edge>> = HashMap::new();
    let mut shell = vec![];
    for (i, (kept, n, dangling, on)) in res.into_iter().enumerate() {
        st.regions += n; st.dangling += dangling; st.on += on; st.kept += kept.len();
        let fc = w.faces.get(faces[i]);
        for (reg, flip) in kept {
            let mut loops = vec![];
            for mut l in reg.loops {
                if flip { l.reverse(); for c in l.iter_mut() { c.1 = !c.1; } }
                let mut coedges = vec![];
                for (e, sense) in l {
                    let ne = *emap.entry(e).or_insert_with(|| {
                        let mut ed = w.edges.get(e).clone();
                        ed.v = ed.v.map(|v| *vmap.entry(v).or_insert_with(|| r.vertices.push(w.vertices.get(v).clone())));
                        r.edges.push(ed)
                    });
                    coedges.push(r.coedges.push(Coedge { edge: ne, sense }));
                }
                loops.push(r.loops.push(Loop { coedges }));
            }
            shell.push(r.faces.push(Face { surface: fc.surface.clone(), sense: fc.sense != flip, loops, tol: fc.tol, origin: Origin { op: OP_BOOLEAN, src: faces[i].0 + 1 } }));
        }
    }
    r.shells.push(Shell { faces: shell, closed: true });
    (r, st)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::make::{block, cylinder};
    const X: V = [1.0, 0.0, 0.0];
    const Y: V = [0.0, 1.0, 0.0];
    const Z: V = [0.0, 0.0, 1.0];
    const O: V = [0.0; 3];

    fn run(name: &str, a: Body, b: Body, op: Op, want: f64, nonmanifold: usize) {
        let t0 = std::time::Instant::now();
        let (r, st) = boolean(&a, &b, op, 4);
        let faces: Vec<Id<Face>> = r.faces.ids().collect();
        let (v, bad) = crate::mesh::volume(&r, &faces);
        let rep = crate::check::check(&r, 5, 4);
        let rel = (v - want).abs() / want;
        eprintln!("{name:<40} {} | volume {v:.6} vs {want:.6} (rel {rel:.1e}) | {st:?} | {rep:?} | {:.0} ms", r.stats(), t0.elapsed().as_secs_f64() * 1e3);
        assert_eq!(st.dangling, 0, "{name}: dangling cut ends");
        assert_eq!(st.on, 0, "{name}: unclassified regions");
        assert_eq!(bad, 0, "{name}: triangulation lost area");
        assert_eq!(rep.violations(), 0, "{name}: checker {rep:?}");
        assert_eq!(rep.nonmanifold, nonmanifold, "{name}: non-manifold edges");
        assert!(rel < 1e-4, "{name}: volume {v} vs {want}");
    }
    #[test] fn s1_disjoint_union() { run("S1 cube U far cube", block(O, [1.0; 3]), block([10.0, 0.0, 0.0], [1.0; 3]), Op::Union, 16.0, 0) }
    #[test] fn s2_cylinder_alone() { run("S2 cylinder U far cube", cylinder(O, 0.5, 2.0, X, Y, Z), block([10.0, 0.0, 0.0], [1.0; 3]), Op::Union, 8.0 + PI, 0) }
    #[test] fn b1_through_hole() { run("B1 cube - cylinder r0.5", block(O, [1.0; 3]), cylinder(O, 0.5, 2.0, X, Y, Z), Op::Difference, 8.0 - PI * 0.5, 0) }
    #[test] fn b2_offset_cubes() { run("B2 cube U offset cube", block(O, [1.0; 3]), block([0.5; 3], [1.0; 3]), Op::Union, 16.0 - 1.5f64.powi(3), 0) }
    #[test] fn b3_steinmetz() { run("B3 cyl ^ cyl (Steinmetz)", cylinder(O, 1.0, 2.0, X, Y, Z), cylinder(O, 1.0, 2.0, Y, Z, X), Op::Intersection, 16.0 / 3.0, 0) }
    #[test] fn b4_tangent_hole() { run("B4 cube - cylinder r1 (tangent)", block(O, [1.0; 3]), cylinder(O, 1.0, 2.0, X, Y, Z), Op::Difference, 8.0 - 2.0 * PI, 4) }
    #[test] fn b5_axis_on_face() { run("B5 cube ^ cylinder on face x=1", block(O, [1.0; 3]), cylinder(X, 0.8, 2.0, X, Y, Z), Op::Intersection, PI * 0.64, 0) }
}
