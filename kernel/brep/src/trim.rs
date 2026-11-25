//! Faces in parameter space. A face's loops become UV polylines in a normalized frame `uv' = (±u, v)` where
//! the sign makes the face's outer loops counter-clockwise (face normal = surface normal when `sense`).
//! Periodic directions are unwrapped continuously along each loop; points are tested on every periodic sheet.
use crate::exact::orient2d;
use crate::geom::*;
use crate::nurbs::tess;
use crate::topo::*;

/// chordal tolerance of edge tessellation (model units)
pub const SAG: f64 = 1e-5;

pub type P2 = [f64; 2];

/// one edge use in UV: polyline in traversal direction, with curve parameters and end vertices
pub struct Item { pub e: Id<Edge>, pub sense: bool, pub uv: Vec<P2>, pub p: Vec<V>, pub t: Vec<f64>, pub v: [Id<Vertex>; 2] }

pub struct Trim {
    /// +1 or -1: u is mirrored for faces whose normal opposes the surface normal
    pub sg: f64,
    /// periods of u and v (0 = not periodic)
    pub per: [f64; 2],
    pub loops: Vec<Vec<Item>>,
    poly: Vec<Vec<P2>>,
    pub bb: [f64; 4],
}

pub fn area(poly: &[P2]) -> f64 { (0..poly.len()).map(|i| { let (p, q) = (poly[i], poly[(i + 1) % poly.len()]); p[0] * q[1] - p[1] * q[0] }).sum::<f64>() / 2.0 }
fn o2(a: P2, b: P2, c: P2) -> i8 { orient2d(a, b, c).0 }
/// even-odd point in polygon
pub fn in_poly(p: P2, poly: &[P2]) -> bool {
    let mut c = false;
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) { c = !c; }
    }
    c
}
pub fn pt_seg(p: P2, a: P2, b: P2) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
    ((a[0] + t * dx - p[0]).powi(2) + (a[1] + t * dy - p[1]).powi(2)).sqrt()
}
/// crossing or touching of segments ab and cd (exact signs; collinear overlaps excluded): parameters along ab and cd
pub fn seg_cross(a: P2, b: P2, c: P2, d: P2) -> Option<(f64, f64)> {
    let (d1, d2, d3, d4) = (o2(c, d, a), o2(c, d, b), o2(a, b, c), o2(a, b, d));
    if d1 * d2 > 0 || d3 * d4 > 0 || (d1 == 0 && d2 == 0) { return None; }
    let den = (b[0] - a[0]) * (d[1] - c[1]) - (b[1] - a[1]) * (d[0] - c[0]);
    if den == 0.0 { return None; }
    let s = ((c[0] - a[0]) * (d[1] - c[1]) - (c[1] - a[1]) * (d[0] - c[0])) / den;
    let t = ((c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0])) / den;
    Some((s.clamp(0.0, 1.0), t.clamp(0.0, 1.0)))
}
pub fn seg_dist(a: P2, b: P2, c: P2, d: P2) -> f64 {
    if seg_cross(a, b, c, d).is_some() { return 0.0; }
    pt_seg(a, c, d).min(pt_seg(b, c, d)).min(pt_seg(c, a, b)).min(pt_seg(d, a, b))
}

impl Trim {
    pub fn new(b: &Body, s: &dyn Surf, sense: bool, loops: &[Vec<(Id<Edge>, bool)>]) -> Trim { Self::with_sag(b, s, sense, loops, SAG) }
    /// edges tessellated to chordal tolerance `sag`
    pub fn with_sag(b: &Body, s: &dyn Surf, sense: bool, loops: &[Vec<(Id<Edge>, bool)>], sag: f64) -> Trim { Self::with_edge_sag(b, s, sense, loops, &|_| sag) }
    /// each edge tessellated to its own chordal tolerance (shared edges get the same value in every face)
    pub fn with_edge_sag(b: &Body, s: &dyn Surf, sense: bool, loops: &[Vec<(Id<Edge>, bool)>], sag: &dyn Fn(Id<Edge>) -> f64) -> Trim {
        let d = s.dom();
        let per = [if s.per_u() { d[1] - d[0] } else { 0.0 }, if s.per_v() { d[3] - d[2] } else { 0.0 }];
        let mut t = Trim { sg: if sense { 1.0 } else { -1.0 }, per, loops: vec![], poly: vec![], bb: [f64::MAX, f64::MIN, f64::MAX, f64::MIN] };
        for l in loops {
            let mut items: Vec<Item> = vec![];
            let mut prev: Option<P2> = None;
            for &(e, sense) in l {
                let ed = b.edges.get(e);
                let v = if sense { ed.v } else { [ed.v[1], ed.v[0]] };
                let (mut ts, mut p) = match &ed.curve {
                    Some(c) => { let ts = tess(c.as_ref(), ed.t[0], ed.t[1], sag(e)); let p = ts.iter().map(|&x| c.ev(x).0).collect(); (ts, p) }
                    None => (vec![ed.t[0], ed.t[1]], vec![b.vertices.get(ed.v[0]).p, b.vertices.get(ed.v[1]).p]),
                };
                if !sense { ts.reverse(); p.reverse(); }
                let n = p.len();
                p[0] = b.vertices.get(v[0]).p;
                p[n - 1] = b.vertices.get(v[1]).p;
                // the first point is the previous edge's last point (the same vertex): identical parameters
                let first = prev;
                let mut uv: Vec<P2> = p.iter().map(|&q| { let w = t.unwrap(prev, t.uv_near(s, q, prev)); prev = Some(w); w }).collect();
                // a pole has no u of its own: it takes the u of its neighbour along this edge, and the loop runs along
                // the pole line from the previous edge's u to this one
                let pole = |q: P2| s.poles_v() && ((q[1] - d[2]).abs() < 1e-9 || (q[1] - d[3]).abs() < 1e-9);
                let n = uv.len();
                if n > 1 {
                    if pole(uv[0]) { uv[0][0] = uv[1][0]; }
                    if pole(uv[n - 1]) { uv[n - 1][0] = uv[n - 2][0]; prev = Some(uv[n - 1]); }
                }
                if let Some(f) = first { if !pole(uv[0]) { uv[0] = f; } }
                items.push(Item { e, sense, uv, p, t: ts, v });
            }
            // closing vertex: the loop ends where it started (unless it wraps a period, e.g. around a cylinder)
            if let (Some(f), Some(l)) = (items.first().map(|i| i.uv[0]), items.last_mut()) {
                let e = l.uv.last_mut().unwrap();
                if (e[0] - f[0]).abs() + (e[1] - f[1]).abs() < 1e-9 * (1.0 + f[0].abs() + f[1].abs()) { *e = f; }
            }
            t.loops.push(items);
        }
        // a single loop that crosses a pole line: the crossing to the far side (u +- half a period) is a tie; the
        // side after the crossing is the one that orients the loop counter-clockwise (the face on its left)
        if t.loops.len() == 1 && per[0] > 0.0 {
            let l = &mut t.loops[0];
            if let Some(j) = (0..l.len()).find(|&i| l[i].uv[l[i].uv.len() - 1] != l[(i + 1) % l.len()].uv[0]) {
                if j + 1 < l.len() && area(&Self::ring(l)) < 0.0 {
                    for sh in [-per[0], per[0]] {
                        for it in l[j + 1..].iter_mut() { for q in it.uv.iter_mut() { q[0] += sh; } }
                        if area(&Self::ring(l)) > 0.0 { break; }
                        for it in l[j + 1..].iter_mut() { for q in it.uv.iter_mut() { q[0] -= sh; } }
                    }
                }
            }
        }
        // loops on periodic surfaces: shift every loop to the sheet of the largest one
        t.poly = t.loops.iter().map(|l| Self::ring(l)).collect();
        if let Some(r) = (0..t.poly.len()).max_by(|&i, &j| area(&t.poly[i]).abs().total_cmp(&area(&t.poly[j]).abs())) {
            let c0 = centroid(&t.poly[r]);
            for i in 0..t.loops.len() {
                let c = centroid(&t.poly[i]);
                let sh = [0, 1].map(|k| if per[k] > 0.0 { -((c[k] - c0[k]) / per[k]).round() * per[k] } else { 0.0 });
                if sh != [0.0, 0.0] { for it in t.loops[i].iter_mut() { for q in it.uv.iter_mut() { q[0] += sh[0]; q[1] += sh[1]; } } }
            }
        }
        t.poly = t.loops.iter().map(|l| Self::ring(l)).collect();
        for q in t.poly.iter().flatten() { t.bb = [t.bb[0].min(q[0]), t.bb[1].max(q[0]), t.bb[2].min(q[1]), t.bb[3].max(q[1])]; }
        t
    }
    /// the loop polygon: each item's points without its last one, unless the next item starts elsewhere (a pole line)
    fn ring(l: &[Item]) -> Vec<P2> {
        (0..l.len()).flat_map(|i| { let it = &l[i]; let keep = it.uv[it.uv.len() - 1] != l[(i + 1) % l.len()].uv[0]; it.uv[..it.uv.len() - (!keep) as usize].iter().copied().collect::<Vec<_>>() }).collect()
    }
    /// normalized parameters of a point (base sheet)
    pub fn uv(&self, s: &dyn Surf, p: V) -> P2 { let q = s.uv_of(p); [self.sg * q[0], q[1]] }
    /// normalized parameters of a point next to a known one (curve walking; see `uv_from`)
    pub fn uv_near(&self, s: &dyn Surf, p: V, prev: Option<P2>) -> P2 {
        let q = uv_from(s, p, prev.map(|g| self.raw(g)), 1e-4);
        [self.sg * q[0], q[1]]
    }
    /// raw surface parameters of a normalized point
    pub fn raw(&self, q: P2) -> P2 { [self.sg * q[0], q[1]] }
    /// q moved by whole periods next to `prev`
    pub fn unwrap(&self, prev: Option<P2>, mut q: P2) -> P2 {
        if let Some(l) = prev { for k in 0..2 { if self.per[k] > 0.0 { q[k] -= ((q[k] - l[k]) / self.per[k]).round() * self.per[k]; } } }
        q
    }
    /// sheet offsets to try for a base-sheet point
    pub fn shifts(&self) -> Vec<P2> {
        let r = |k: usize| if self.per[k] > 0.0 { vec![-self.per[k], 0.0, self.per[k], 2.0 * self.per[k], -2.0 * self.per[k]] } else { vec![0.0] };
        r(0).into_iter().flat_map(|a| r(1).into_iter().map(move |b| [a, b])).collect()
    }
    fn inside_here(&self, q: P2) -> bool {
        q[0] >= self.bb[0] && q[0] <= self.bb[1] && q[1] >= self.bb[2] && q[1] <= self.bb[3] && self.poly.iter().filter(|p| in_poly(q, p)).count() % 2 == 1
    }
    /// the sheet offset that puts q strictly inside the face, if any
    pub fn shift_into(&self, q: P2) -> Option<P2> { self.shifts().into_iter().find(|s| self.inside_here([q[0] + s[0], q[1] + s[1]])) }
    pub fn inside(&self, q: P2) -> bool { self.shift_into(q).is_some() }
    /// distance from q (any sheet) to the boundary
    pub fn boundary_dist(&self, q: P2) -> f64 {
        self.shifts().into_iter().map(|s| { let q = [q[0] + s[0], q[1] + s[1]]; self.poly.iter().flat_map(|p| (0..p.len()).map(move |i| pt_seg(q, p[i], p[(i + 1) % p.len()]))).fold(f64::MAX, f64::min) }).fold(f64::MAX, f64::min)
    }
    /// inside, or on the boundary within `tol`
    pub fn covers(&self, q: P2, tol: f64) -> bool { self.inside(q) || self.boundary_dist(q) <= tol }
    /// a 3D polyline on the surface in normalized UV, unwrapped and placed on the sheet where its middle is in the face
    pub fn polyline(&self, s: &dyn Surf, pts: &[V]) -> Vec<P2> {
        let mut prev = None;
        let mut uv: Vec<P2> = pts.iter().map(|&p| { let w = self.unwrap(prev, self.uv_near(s, p, prev)); prev = Some(w); w }).collect();
        let m = uv.len() / 2;
        let mid = if uv.len() % 2 == 0 { [(uv[m - 1][0] + uv[m][0]) / 2.0, (uv[m - 1][1] + uv[m][1]) / 2.0] } else { uv[m] };
        let sh = self.shift_into(mid).or_else(|| self.shifts().into_iter().min_by(|a, b| {
            let d = |s: &P2| self.poly.iter().flatten().map(|p| (p[0] - mid[0] - s[0]).hypot(p[1] - mid[1] - s[1])).fold(f64::MAX, f64::min);
            d(a).total_cmp(&d(b))
        })).unwrap_or([0.0, 0.0]);
        for q in uv.iter_mut() { q[0] += sh[0]; q[1] += sh[1]; }
        uv
    }
    pub fn polygons(&self) -> &[Vec<P2>] { &self.poly }
    /// edges with a segment whose UV midpoint maps farther than `tol` from the segment's 3D midpoint (the straight UV
    /// segment does not follow the 3D chord, e.g. near a pole): a face mesh along them cannot meet `tol`
    pub fn bent_edges(&self, s: &dyn Surf, tol: f64) -> Vec<Id<Edge>> {
        let mut out: Vec<Id<Edge>> = self.loops.iter().flatten().filter(|it| it.uv.windows(2).zip(it.p.windows(2)).any(|(w, q)| {
            let r = self.raw([(w[0][0] + w[1][0]) / 2.0, (w[0][1] + w[1][1]) / 2.0]);
            norm(sub(s.ev(r[0], r[1]).p, mulv(addv(q[0], q[1]), 0.5))) > tol
        })).map(|it| it.e).collect();
        out.sort_unstable_by_key(|e| e.0);
        out.dedup();
        out
    }
    /// edges whose tessellated segments cross or touch a non-adjacent segment of any loop (loops closer than the
    /// tessellation tolerance): these need a finer tessellation
    pub fn crossing_edges(&self) -> Vec<Id<Edge>> {
        let segs: Vec<(P2, P2, Id<Edge>)> = self.loops.iter().flatten().flat_map(|it| it.uv.windows(2).map(move |w| (w[0], w[1], it.e))).collect();
        let mut idx: Vec<usize> = (0..segs.len()).collect();
        idx.sort_by(|&i, &j| segs[i].0[0].min(segs[i].1[0]).total_cmp(&segs[j].0[0].min(segs[j].1[0])));
        let mut out = vec![];
        for (k, &i) in idx.iter().enumerate() {
            let (a, b, ea) = segs[i];
            let xmax = a[0].max(b[0]);
            for &j in &idx[k + 1..] {
                let (c, d, eb) = segs[j];
                if c[0].min(d[0]) > xmax { break; }
                // segments sharing an end point meet there by construction
                if a == c || a == d || b == c || b == d { continue; }
                let (d1, d2, d3, d4) = (orient2d(c, d, a).0, orient2d(c, d, b).0, orient2d(a, b, c).0, orient2d(a, b, d).0);
                if d1 * d2 <= 0 && d3 * d4 <= 0 && !(d1 == 0 && d2 == 0) {
                    out.push(ea); out.push(eb);
                }
            }
        }
        out.sort_unstable_by_key(|e| e.0);
        out.dedup();
        out
    }
}
fn centroid(p: &[P2]) -> P2 { let n = p.len().max(1) as f64; [p.iter().map(|q| q[0]).sum::<f64>() / n, p.iter().map(|q| q[1]).sum::<f64>() / n] }
