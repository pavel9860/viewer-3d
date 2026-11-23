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
    pub fn new(b: &Body, s: &dyn Surf, sense: bool, loops: &[Vec<(Id<Edge>, bool)>]) -> Trim {
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
                    Some(c) => { let ts = tess(c.as_ref(), ed.t[0], ed.t[1], SAG); let p = ts.iter().map(|&x| c.ev(x).0).collect(); (ts, p) }
                    None => (vec![ed.t[0], ed.t[1]], vec![b.vertices.get(ed.v[0]).p, b.vertices.get(ed.v[1]).p]),
                };
                if !sense { ts.reverse(); p.reverse(); }
                let n = p.len();
                p[0] = b.vertices.get(v[0]).p;
                p[n - 1] = b.vertices.get(v[1]).p;
                let uv = p.iter().map(|&q| { let w = t.unwrap(prev, t.uv(s, q)); prev = Some(w); w }).collect();
                items.push(Item { e, sense, uv, p, t: ts, v });
            }
            t.loops.push(items);
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
    fn ring(l: &[Item]) -> Vec<P2> { l.iter().flat_map(|it| it.uv[..it.uv.len() - 1].iter().copied()).collect() }
    /// normalized parameters of a point (base sheet)
    pub fn uv(&self, s: &dyn Surf, p: V) -> P2 { let q = s.uv_of(p); [self.sg * q[0], q[1]] }
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
        let mut uv: Vec<P2> = pts.iter().map(|&p| { let w = self.unwrap(prev, self.uv(s, p)); prev = Some(w); w }).collect();
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
}
fn centroid(p: &[P2]) -> P2 { let n = p.len().max(1) as f64; [p.iter().map(|q| q[0]).sum::<f64>() / n, p.iter().map(|q| q[1]).sum::<f64>() / n] }
