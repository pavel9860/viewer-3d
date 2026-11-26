//! Detail levels of a triangle mesh by error-bounded edge collapse.
//!
//! The collapse order is by quadric error (Garland-Heckbert, unweighted triangle planes, so the square root of a
//! vertex's quadric bounds its distance to every merged plane). A collapse is refused when it would change the
//! topology (link condition, duplicate triangles), flip or degenerate a triangle, or move a boundary or non-manifold
//! vertex: a closed manifold mesh stays closed and manifold at every level. Collapses run until the next one would
//! exceed a level's tolerance; the mesh is then saved as that level and collapsing continues for the next.
//! The error of each level is measured, not estimated: the larger of the distances from the original vertices and
//! triangle centroids to the level and from the level's vertices to the original.
use crate::geom::*;
use crate::par;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// one level: positions, triangles, the original triangle's tag (face id) per triangle, measured deviation
#[derive(Clone, Debug, Default)]
pub struct Level { pub pos: Vec<V>, pub tri: Vec<[u32; 3]>, pub tag: Vec<u32>, pub err: f64 }

/// symmetric 4x4 quadric: a00 a01 a02 a11 a12 a22 b0 b1 b2 c, cost(v) = v.A.v + 2 b.v + c
#[derive(Clone, Copy, Default)]
struct Q([f64; 10]);
impl Q {
    fn plane(n: V, d: f64) -> Q { Q([n[0] * n[0], n[0] * n[1], n[0] * n[2], n[1] * n[1], n[1] * n[2], n[2] * n[2], n[0] * d, n[1] * d, n[2] * d, d * d]) }
    fn add(&mut self, o: &Q) { for k in 0..10 { self.0[k] += o.0[k]; } }
    fn cost(&self, v: V) -> f64 {
        let q = &self.0;
        let c = v[0] * (q[0] * v[0] + 2.0 * q[1] * v[1] + 2.0 * q[2] * v[2]) + v[1] * (q[3] * v[1] + 2.0 * q[4] * v[2]) + q[5] * v[2] * v[2]
            + 2.0 * (q[6] * v[0] + q[7] * v[1] + q[8] * v[2]) + q[9];
        c.max(0.0)
    }
    /// minimizer, when A is well conditioned
    fn optimum(&self) -> Option<V> {
        let q = &self.0;
        let a = [[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], q[5]]];
        let det = dot(a[0], cross(a[1], a[2]));
        let tr = q[0] + q[3] + q[5];
        if det.abs() <= 1e-9 * tr * tr * tr { return None; }
        let b = [-q[6], -q[7], -q[8]];
        let r = |i: usize| { let mut m = a; for (k, row) in m.iter_mut().enumerate() { row[i] = b[k]; } dot(m[0], cross(m[1], m[2])) / det };
        Some([r(0), r(1), r(2)])
    }
}

struct S<'a> {
    pos: Vec<V>,
    tri: Vec<[u32; 3]>,
    alive: Vec<bool>,
    vt: Vec<Vec<u32>>,
    q: Vec<Q>,
    lock: Vec<bool>,
    stamp: Vec<u32>,
    gone: Vec<bool>,
    tag: &'a [u32],
}

fn normal(a: V, b: V, c: V) -> V { cross(sub(b, a), sub(c, a)) }

impl S<'_> {
    fn ring(&self, v: u32) -> Vec<u32> {
        let mut r: Vec<u32> = self.vt[v as usize].iter().flat_map(|&t| self.tri[t as usize]).filter(|&x| x != v).collect();
        r.sort_unstable();
        r.dedup();
        r
    }
    /// best position for merging a and b and its cost
    fn target(&self, a: u32, b: u32) -> (f64, V) {
        let mut q = self.q[a as usize];
        q.add(&self.q[b as usize]);
        let (pa, pb) = (self.pos[a as usize], self.pos[b as usize]);
        let mut best = (q.cost(pb), pb);
        let mut try_ = |p: V| { let c = q.cost(p); if c < best.0 { best = (c, p); } };
        try_(pa);
        try_(mulv(addv(pa, pb), 0.5));
        // the optimum only near the edge (far optima of nearly flat regions would stretch the mesh)
        if let Some(p) = q.optimum() {
            let l2 = dot(sub(pa, pb), sub(pa, pb));
            if dot(sub(p, pa), sub(p, pa)) <= l2 && dot(sub(p, pb), sub(p, pb)) <= l2 { try_(p); }
        }
        best
    }
    /// can a and b merge at p without changing the topology or folding the surface
    fn valid(&self, a: u32, b: u32, p: V) -> bool {
        if self.lock[a as usize] || self.lock[b as usize] { return false; }
        // bounded valence: keeps the checks local and the triangles well shaped
        if self.vt[a as usize].len() + self.vt[b as usize].len() > MAX_FAN + 2 { return false; }
        let shared: Vec<u32> = self.vt[a as usize].iter().copied().filter(|t| self.tri[*t as usize].contains(&b)).collect();
        if shared.len() != 2 { return false; }
        let opp: Vec<u32> = shared.iter().map(|&t| self.tri[t as usize].into_iter().find(|&x| x != a && x != b).unwrap()).collect();
        let (ra, rb) = (self.ring(a), self.ring(b));
        let common = ra.iter().filter(|x| rb.binary_search(x).is_ok()).count();
        if common != 2 || opp[0] == opp[1] { return false; }
        // every other triangle of either fan keeps its orientation and stays non-degenerate
        for (v, o) in [(a, b), (b, a)] {
            for &t in &self.vt[v as usize] {
                let tr = self.tri[t as usize];
                if tr.contains(&o) { continue; }
                let old = normal(self.pos[tr[0] as usize], self.pos[tr[1] as usize], self.pos[tr[2] as usize]);
                let q = tr.map(|x| if x == v { p } else { self.pos[x as usize] });
                let new = normal(q[0], q[1], q[2]);
                let (lo, ln) = (norm(old), norm(new));
                if ln <= 1e-12 * (dot(sub(q[1], q[0]), sub(q[1], q[0])) + dot(sub(q[2], q[0]), sub(q[2], q[0]))) { return false; }
                if dot(old, new) <= 0.2 * lo * ln { return false; }
            }
        }
        // a's triangles, re-attached to b, must not duplicate one of b's
        for &t in &self.vt[a as usize] {
            let tr = self.tri[t as usize];
            if tr.contains(&b) { continue; }
            let mut k = tr.map(|x| if x == a { b } else { x });
            k.sort_unstable();
            for &u in &self.vt[b as usize] {
                let mut w = self.tri[u as usize];
                w.sort_unstable();
                if w == k { return false; }
            }
        }
        true
    }
    fn collapse(&mut self, a: u32, b: u32, p: V) {
        let fa = std::mem::take(&mut self.vt[a as usize]);
        for t in fa {
            if self.tri[t as usize].contains(&b) {
                self.alive[t as usize] = false;
                for x in self.tri[t as usize] { if x != a { self.vt[x as usize].retain(|&u| u != t); } }
            } else {
                for x in self.tri[t as usize].iter_mut() { if *x == a { *x = b; } }
                self.vt[b as usize].push(t);
            }
        }
        self.pos[b as usize] = p;
        let qa = self.q[a as usize];
        self.q[b as usize].add(&qa);
        self.gone[a as usize] = true;
        // only the merged vertex's quadric and position changed: edges between other vertices keep their cost
        // (their validity is checked when they come up)
        self.stamp[b as usize] += 1;
        self.stamp[a as usize] += 1;
    }
    fn level(&self) -> Level {
        let mut map = vec![u32::MAX; self.pos.len()];
        let mut l = Level::default();
        for (t, tr) in self.tri.iter().enumerate() {
            if !self.alive[t] { continue; }
            l.tri.push(tr.map(|v| {
                if map[v as usize] == u32::MAX { map[v as usize] = l.pos.len() as u32; l.pos.push(self.pos[v as usize]); }
                map[v as usize]
            }));
            l.tag.push(self.tag[t]);
        }
        l
    }
}

/// most triangles around a vertex after a collapse
const MAX_FAN: usize = 24;

type Entry = Reverse<(u64, u32, u32, u32, u32)>;

fn push(s: &S, h: &mut BinaryHeap<Entry>, a: u32, b: u32) {
    let (a, b) = (a.min(b), a.max(b));
    if s.lock[a as usize] || s.lock[b as usize] { return; }
    let (c, _) = s.target(a, b);
    h.push(Reverse((c.to_bits(), a, b, s.stamp[a as usize], s.stamp[b as usize])));
}

/// Levels of increasing tolerance `tols` (model units; a level is kept only when it removed triangles). Every
/// collapse of level k has quadric error at most tols[k]; the reported error is measured.
pub fn levels(pos: &[V], tri: &[[u32; 3]], tag: &[u32], tols: &[f64], workers: usize) -> Vec<Level> {
    let nv = pos.len();
    let mut s = S { pos: pos.to_vec(), tri: tri.to_vec(), alive: vec![true; tri.len()], vt: vec![vec![]; nv], q: vec![Q::default(); nv], lock: vec![false; nv], stamp: vec![0; nv], gone: vec![false; nv], tag };
    let mut edges: std::collections::HashMap<(u32, u32), u32> = Default::default();
    for (t, tr) in tri.iter().enumerate() {
        let n = normal(pos[tr[0] as usize], pos[tr[1] as usize], pos[tr[2] as usize]);
        let l = norm(n);
        for k in 0..3 {
            s.vt[tr[k] as usize].push(t as u32);
            let (a, b) = (tr[k], tr[(k + 1) % 3]);
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
        if l > 0.0 {
            let n = mulv(n, 1.0 / l);
            let pl = Q::plane(n, -dot(n, pos[tr[0] as usize]));
            for v in tr { s.q[*v as usize].add(&pl); }
        }
    }
    // boundary and non-manifold edges keep their vertices
    for (&(a, b), &c) in &edges { if c != 2 { s.lock[a as usize] = true; s.lock[b as usize] = true; } }
    let mut keys: Vec<_> = edges.keys().copied().collect();
    keys.sort_unstable();
    let mut h = BinaryHeap::new();
    for (a, b) in keys { push(&s, &mut h, a, b); }
    let mut out: Vec<Level> = vec![];
    let mut last = tri.len();
    for &tol in tols {
        let lim = tol * tol;
        while let Some(&Reverse((c, a, b, sa, sb))) = h.peek() {
            if f64::from_bits(c) > lim { break; }
            h.pop();
            if s.gone[a as usize] || s.gone[b as usize] || s.stamp[a as usize] != sa || s.stamp[b as usize] != sb { continue; }
            let (_, p) = s.target(a, b);
            if !s.valid(a, b, p) { continue; }
            // the survivor is the lower index (a deterministic choice; the position is p either way)
            let (gone, keep) = (b, a);
            s.collapse(gone, keep, p);
            for n in s.ring(keep) { push(&s, &mut h, keep, n); }
        }
        let n = s.alive.iter().filter(|x| **x).count();
        if n < last {
            last = n;
            let t0 = std::time::Instant::now();
            let mut l = s.level();
            l.err = deviation(pos, tri, &l.pos, &l.tri, workers);
            if std::env::var("BREP_DEBUG").is_ok() { eprintln!("level tol {tol:.1e}: {n} triangles, deviation {:.2e} in {:.0} ms", l.err, t0.elapsed().as_secs_f64() * 1e3); }
            out.push(l);
        }
    }
    out
}

/// symmetric deviation sampled at vertices and centroids: max over the samples of one mesh of the distance to the
/// other mesh
pub fn deviation(pa: &[V], ta: &[[u32; 3]], pb: &[V], tb: &[[u32; 3]], workers: usize) -> f64 {
    let samples = |p: &[V], t: &[[u32; 3]]| -> Vec<V> {
        let mut s = p.to_vec();
        s.extend(t.iter().map(|x| mulv(addv(addv(p[x[0] as usize], p[x[1] as usize]), p[x[2] as usize]), 1.0 / 3.0)));
        s
    };
    let (sa, sb) = (samples(pa, ta), samples(pb, tb));
    let (ga, gb) = (Grid::new(pa, ta), Grid::new(pb, tb));
    let far = |pts: &[V], g: &Grid| -> f64 {
        let n = pts.len().div_ceil(4096);
        par::map(n, workers, |i| pts[4096 * i..(4096 * (i + 1)).min(pts.len())].iter().map(|&q| g.dist(q)).fold(0.0, f64::max)).into_iter().fold(0.0, f64::max)
    };
    far(&sa, &gb).max(far(&sb, &ga))
}

/// bounding-box tree over triangles for nearest-triangle distance queries (branch and bound, nearer child first)
struct Grid<'a> { p: &'a [V], t: &'a [[u32; 3]], nodes: Vec<([f64; 6], u32, u32)>, order: Vec<u32> }

const LEAF: usize = 4;

impl<'a> Grid<'a> {
    fn new(p: &'a [V], t: &'a [[u32; 3]]) -> Grid<'a> {
        let boxes: Vec<[f64; 6]> = t.iter().map(|x| {
            let v = x.map(|i| p[i as usize]);
            let mut b = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
            for q in v { for k in 0..3 { b[k] = b[k].min(q[k]); b[k + 3] = b[k + 3].max(q[k]); } }
            b
        }).collect();
        let mut g = Grid { p, t, nodes: vec![], order: (0..t.len() as u32).collect() };
        // node: box, then (first child, 0) for an inner node or (start, count | 1 << 31) for a leaf
        let mut stack = vec![(0usize, t.len(), 0usize)];
        g.nodes.push(([0.0; 6], 0, 0));
        while let Some((lo, hi, ni)) = stack.pop() {
            let mut bb = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
            for &i in &g.order[lo..hi] { let b = &boxes[i as usize]; for k in 0..3 { bb[k] = bb[k].min(b[k]); bb[k + 3] = bb[k + 3].max(b[k + 3]); } }
            if hi - lo <= LEAF { g.nodes[ni] = (bb, lo as u32, (hi - lo) as u32 | 1 << 31); continue; }
            let ax = (0..3).max_by(|&i, &j| (bb[i + 3] - bb[i]).total_cmp(&(bb[j + 3] - bb[j]))).unwrap();
            let mid = (lo + hi) / 2;
            let c = |i: &u32| boxes[*i as usize][ax] + boxes[*i as usize][ax + 3];
            g.order[lo..hi].select_nth_unstable_by(mid - lo, |a, b| c(a).total_cmp(&c(b)));
            let first = g.nodes.len();
            g.nodes.push(([0.0; 6], 0, 0));
            g.nodes.push(([0.0; 6], 0, 0));
            g.nodes[ni] = (bb, first as u32, 0);
            stack.push((lo, mid, first));
            stack.push((mid, hi, first + 1));
        }
        g
    }
    fn dist(&self, q: V) -> f64 {
        let bd = |b: &[f64; 6]| (0..3).map(|k| { let d = (b[k] - q[k]).max(q[k] - b[k + 3]).max(0.0); d * d }).sum::<f64>();
        let mut best = f64::MAX;
        let mut stack = vec![(0u32, 0.0f64)];
        while let Some((ni, d)) = stack.pop() {
            if d >= best { continue; }
            let (_, a, n) = self.nodes[ni as usize];
            if n & 1 << 31 != 0 {
                for &ti in &self.order[a as usize..(a + (n & !(1 << 31))) as usize] {
                    let t = self.t[ti as usize];
                    best = best.min(point_tri2(q, self.p[t[0] as usize], self.p[t[1] as usize], self.p[t[2] as usize]));
                }
                continue;
            }
            let (d0, d1) = (bd(&self.nodes[a as usize].0), bd(&self.nodes[a as usize + 1].0));
            // the nearer child is popped first
            if d0 < d1 { stack.push((a + 1, d1)); stack.push((a, d0)); } else { stack.push((a, d0)); stack.push((a + 1, d1)); }
        }
        best.sqrt()
    }
}

/// squared distance from p to triangle abc (Ericson, Real-Time Collision Detection 5.1.5)
fn point_tri2(p: V, a: V, b: V, c: V) -> f64 {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let d2 = |x: V| dot(sub(p, x), sub(p, x));
    let (d1, d2_) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2_ <= 0.0 { return d2(a); }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 { return d2(b); }
    let vc = d1 * d4 - d3 * d2_;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 { let v = d1 / (d1 - d3); return d2(addv(a, mulv(ab, v))); }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 { return d2(c); }
    let vb = d5 * d2_ - d1 * d6;
    if vb <= 0.0 && d2_ >= 0.0 && d6 <= 0.0 { let w = d2_ / (d2_ - d6); return d2(addv(a, mulv(ac, w))); }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 { let w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); return d2(addv(b, mulv(sub(c, b), w))); }
    let den = 1.0 / (va + vb + vc);
    let (v, w) = (vb * den, vc * den);
    d2(addv(a, addv(mulv(ab, v), mulv(ac, w))))
}
