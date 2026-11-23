//! NURBS curves/surfaces, extrusion and revolution surfaces, with derivatives up to order 2.
use crate::geom::*;
use std::f64::consts::PI;

pub fn span(n: usize, p: usize, u: f64, k: &[f64]) -> usize {
    if u >= k[n + 1] { return n; }
    if u <= k[p] { return p; }
    let (mut lo, mut hi) = (p, n + 1);
    let mut mid = (lo + hi) / 2;
    while u < k[mid] || u >= k[mid + 1] { if u < k[mid] { hi = mid; } else { lo = mid; } mid = (lo + hi) / 2; }
    mid
}
/// largest supported degree
pub const MAXP: usize = 15;
/// basis functions and their derivatives up to order 2 (The NURBS Book, A2.3), allocation-free; entries 0..=p
pub fn ders(i: usize, u: f64, p: usize, k: &[f64]) -> [[f64; MAXP + 1]; 3] {
    assert!(p <= MAXP, "B-spline degree {p} > {MAXP}");
    let mut ndu = [[0.0; MAXP + 1]; MAXP + 1];
    let (mut left, mut right) = ([0.0; MAXP + 1], [0.0; MAXP + 1]);
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - k[i + 1 - j]; right[j] = k[i + j] - u;
        let mut saved = 0.0;
        for r in 0..j { ndu[j][r] = right[r + 1] + left[j - r]; let t = if ndu[j][r] != 0.0 { ndu[r][j - 1] / ndu[j][r] } else { 0.0 }; ndu[r][j] = saved + right[r + 1] * t; saved = left[j - r] * t; }
        ndu[j][j] = saved;
    }
    let nd = 2.min(p);
    let mut out = [[0.0; MAXP + 1]; 3];
    for j in 0..=p { out[0][j] = ndu[j][p]; }
    let mut a = [[0.0; MAXP + 1]; 2];
    for r in 0..=p {
        let (mut s1, mut s2) = (0, 1);
        a[0][0] = 1.0;
        for kk in 1..=nd {
            let mut d = 0.0;
            let rk = r as isize - kk as isize;
            let pk = p - kk;
            if r >= kk { a[s2][0] = a[s1][0] / ndu[pk + 1][rk as usize]; d = a[s2][0] * ndu[rk as usize][pk]; }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk as isize { kk - 1 } else { p - r };
            for j in j1..=j2 { let idx = (rk + j as isize) as usize; a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx]; d += a[s2][j] * ndu[idx][pk]; }
            if r <= pk { a[s2][kk] = -a[s1][kk - 1] / ndu[pk + 1][r]; d += a[s2][kk] * ndu[r][pk]; }
            out[kk][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut f = p as f64;
    for kk in 1..=nd { for j in 0..=p { out[kk][j] *= f; } f *= (p - kk) as f64; }
    out
}

/// rational B-spline surface, control points homogeneous (wx, wy, wz, w)
pub struct BspS { pub pu: usize, pub pv: usize, pub ku: Vec<f64>, pub kv: Vec<f64>, pub cp: Vec<Vec<[f64; 4]>> }
impl Surf for BspS {
    fn ev(&self, u: f64, v: f64) -> D2 {
        // defined on its knot range only: evaluate at the clamped parameter
        let d = self.dom();
        let (u, v) = (u.clamp(d[0], d[1]), v.clamp(d[2], d[3]));
        let (nu, nv) = (self.cp.len() - 1, self.cp[0].len() - 1);
        let (su, sv) = (span(nu, self.pu, u, &self.ku), span(nv, self.pv, v, &self.kv));
        let (bu, bv) = (ders(su, u, self.pu, &self.ku), ders(sv, v, self.pv, &self.kv));
        let mut a = [[[0.0; 4]; 3]; 3];
        for i in 0..=self.pu { for j in 0..=self.pv {
            let c = self.cp[su - self.pu + i][sv - self.pv + j];
            for du in 0..3 { for dv in 0..3 - du { let w = bu[du][i] * bv[dv][j]; for k in 0..4 { a[du][dv][k] += w * c[k]; } } }
        } }
        let g = |x: [f64; 4]| [x[0], x[1], x[2]];
        let w = a[0][0][3];
        let p = mulv(g(a[0][0]), 1.0 / w);
        let su_ = mulv(sub(g(a[1][0]), mulv(p, a[1][0][3])), 1.0 / w);
        let sv_ = mulv(sub(g(a[0][1]), mulv(p, a[0][1][3])), 1.0 / w);
        let suu = mulv(sub(sub(g(a[2][0]), mulv(su_, 2.0 * a[1][0][3])), mulv(p, a[2][0][3])), 1.0 / w);
        let svv = mulv(sub(sub(g(a[0][2]), mulv(sv_, 2.0 * a[0][1][3])), mulv(p, a[0][2][3])), 1.0 / w);
        let suv = mulv(sub(sub(sub(g(a[1][1]), mulv(su_, a[0][1][3])), mulv(sv_, a[1][0][3])), mulv(p, a[1][1][3])), 1.0 / w);
        D2 { p, u: su_, v: sv_, uu: suu, uv: suv, vv: svv }
    }
    fn dom(&self) -> [f64; 4] { [self.ku[self.pu], self.ku[self.ku.len() - 1 - self.pu], self.kv[self.pv], self.kv[self.kv.len() - 1 - self.pv]] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// curve with first and second derivative
pub trait Crv: Send + Sync {
    fn ev(&self, t: f64) -> (V, V, V);
    fn range(&self) -> (f64, f64);
    fn periodic(&self) -> bool { false }
    /// closed-form parameter of the closest point, where one exists
    fn param_of(&self, _p: V) -> Option<f64> { None }
    /// parameters that must appear in any tessellation of [t0, t1] (polyline breakpoints)
    fn breaks(&self, _t0: f64, _t1: f64) -> Option<Vec<f64>> { None }
}

/// tessellation parameters of [t0, t1] (either direction): the curve's own breakpoints, else bisection of
/// 4 spans until each chord is within `sag` of the curve at its midpoint
pub fn tess(c: &dyn Crv, t0: f64, t1: f64, sag: f64) -> Vec<f64> {
    if let Some(b) = c.breaks(t0, t1) { return b; }
    fn rec(c: &dyn Crv, a: (f64, V), b: (f64, V), sag: f64, depth: u32, out: &mut Vec<f64>) {
        let tm = 0.5 * (a.0 + b.0);
        let pm = c.ev(tm).0;
        if depth < 12 && norm(sub(pm, mulv(addv(a.1, b.1), 0.5))) > sag {
            rec(c, a, (tm, pm), sag, depth + 1, out);
            rec(c, (tm, pm), b, sag, depth + 1, out);
        } else { out.push(b.0); }
    }
    let mut out = vec![t0];
    for k in 0..4 {
        let (a, b) = (t0 + (t1 - t0) * k as f64 / 4.0, t0 + (t1 - t0) * (k + 1) as f64 / 4.0);
        rec(c, (a, c.ev(a).0), (b, c.ev(b).0), sag, 0, &mut out);
    }
    out
}

/// polyline parameterized by arc length: intersection curves, exact to their edge tolerance
pub struct PolyC { pub p: Vec<V>, pub s: Vec<f64> }
impl PolyC {
    pub fn new(pts: Vec<V>) -> Self {
        let mut p: Vec<V> = Vec::with_capacity(pts.len());
        for q in pts { if p.last().map_or(true, |l| norm(sub(q, *l)) > 0.0) { p.push(q); } }
        let mut s = vec![0.0];
        for w in p.windows(2) { s.push(s[s.len() - 1] + norm(sub(w[1], w[0]))); }
        PolyC { p, s }
    }
    fn seg(&self, t: f64) -> usize { self.s.partition_point(|&x| x <= t).clamp(1, self.p.len() - 1) - 1 }
}
impl Crv for PolyC {
    fn ev(&self, t: f64) -> (V, V, V) {
        let i = self.seg(t);
        let d = mulv(sub(self.p[i + 1], self.p[i]), 1.0 / (self.s[i + 1] - self.s[i]));
        (addv(self.p[i], mulv(d, t - self.s[i])), d, [0.0; 3])
    }
    fn range(&self) -> (f64, f64) { (0.0, self.s[self.s.len() - 1]) }
    fn param_of(&self, q: V) -> Option<f64> {
        (0..self.p.len() - 1).map(|i| {
            let ab = sub(self.p[i + 1], self.p[i]);
            let u = (dot(sub(q, self.p[i]), ab) / dot(ab, ab)).clamp(0.0, 1.0);
            (norm(sub(q, addv(self.p[i], mulv(ab, u)))), self.s[i] + u * (self.s[i + 1] - self.s[i]))
        }).min_by(|a, b| a.0.total_cmp(&b.0)).map(|x| x.1)
    }
    fn breaks(&self, t0: f64, t1: f64) -> Option<Vec<f64>> {
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let mut b: Vec<f64> = std::iter::once(lo).chain(self.s.iter().copied().filter(|&x| x > lo && x < hi)).chain(std::iter::once(hi)).collect();
        if t0 > t1 { b.reverse(); }
        Some(b)
    }
}
pub struct LineC { pub o: V, pub d: V }
impl Crv for LineC {
    fn ev(&self, t: f64) -> (V, V, V) { (addv(self.o, mulv(self.d, t)), self.d, [0.0; 3]) }
    fn range(&self) -> (f64, f64) { (-1e300, 1e300) }
    fn param_of(&self, p: V) -> Option<f64> { Some(dot(sub(p, self.o), self.d) / dot(self.d, self.d)) }
}
/// circle (a = b) or ellipse
pub struct ConicC { pub c: V, pub x: V, pub y: V, pub a: f64, pub b: f64 }
impl Crv for ConicC {
    fn ev(&self, t: f64) -> (V, V, V) { let (s, c) = t.sin_cos(); (addv(self.c, addv(mulv(self.x, self.a * c), mulv(self.y, self.b * s))), addv(mulv(self.x, -self.a * s), mulv(self.y, self.b * c)), addv(mulv(self.x, -self.a * c), mulv(self.y, -self.b * s))) }
    fn range(&self) -> (f64, f64) { (0.0, 2.0 * PI) }
    fn periodic(&self) -> bool { true }
    fn param_of(&self, p: V) -> Option<f64> { let q = sub(p, self.c); Some((dot(q, self.y) / self.b).atan2(dot(q, self.x) / self.a).rem_euclid(2.0 * PI)) }
}
pub struct BspC { pub p: usize, pub k: Vec<f64>, pub cp: Vec<[f64; 4]> }
impl Crv for BspC {
    fn ev(&self, t: f64) -> (V, V, V) {
        let (r0, r1) = self.range();
        let t = t.clamp(r0, r1);
        let n = self.cp.len() - 1;
        let s = span(n, self.p, t, &self.k);
        let b = ders(s, t, self.p, &self.k);
        let mut a = [[0.0; 4]; 3];
        for i in 0..=self.p { let c = self.cp[s - self.p + i]; for d in 0..3 { for k in 0..4 { a[d][k] += b[d][i] * c[k]; } } }
        let g = |x: [f64; 4]| [x[0], x[1], x[2]];
        let w = a[0][3];
        let p = mulv(g(a[0]), 1.0 / w);
        let d1 = mulv(sub(g(a[1]), mulv(p, a[1][3])), 1.0 / w);
        let d2 = mulv(sub(sub(g(a[2]), mulv(d1, 2.0 * a[1][3])), mulv(p, a[2][3])), 1.0 / w);
        (p, d1, d2)
    }
    fn range(&self) -> (f64, f64) { (self.k[self.p], self.k[self.k.len() - 1 - self.p]) }
}
/// closest curve parameter on the curve's range (closed form where the curve provides one)
pub fn t_of(c: &dyn Crv, p: V) -> f64 {
    if let Some(t) = c.param_of(p) { return t; }
    let (a, b) = c.range();
    let n = 200;
    let mut best = (a, f64::MAX);
    for i in 0..=n { let t = a + (b - a) * i as f64 / n as f64; let d = norm(sub(c.ev(t).0, p)); if d < best.1 { best = (t, d); } }
    let mut t = best.0;
    for _ in 0..50 { let (q, d1, d2) = c.ev(t); let r = sub(q, p); let f = dot(r, d1); let fp = dot(d1, d1) + dot(r, d2); if !(fp.abs() > 1e-300) { break; } let nt = (t - f / fp).clamp(a, b); if (nt - t).abs() < 1e-16 { break; } t = nt; }
    t
}

/// surface of linear extrusion: c(u) + v d, u bounded by the face
pub struct ExtrS { pub c: Box<dyn Crv>, pub d: V, pub u0: f64, pub u1: f64, pub v0: f64, pub v1: f64 }
impl Surf for ExtrS {
    fn ev(&self, u: f64, v: f64) -> D2 { let (p, d1, d2) = self.c.ev(u); D2 { p: addv(p, mulv(self.d, v)), u: d1, v: self.d, uu: d2, uv: [0.0; 3], vv: [0.0; 3] } }
    fn dom(&self) -> [f64; 4] { [self.u0, self.u1, self.v0, self.v1] }
    fn per_u(&self) -> bool { self.c.periodic() }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}
/// surface of revolution of profile c about the axis (a, z), v = angle
pub struct RevS { pub c: Box<dyn Crv>, pub a: V, pub z: V }
impl Surf for RevS {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (p, d1, d2) = self.c.ev(u);
        let (s, c) = v.sin_cos();
        let z = self.z;
        let rot = |q: V| addv(addv(mulv(q, c), mulv(cross(z, q), s)), mulv(z, dot(z, q) * (1.0 - c)));
        let drot = |q: V| addv(mulv(q, -s), addv(mulv(cross(z, q), c), mulv(z, dot(z, q) * s)));
        let ddrot = |q: V| addv(mulv(q, -c), addv(mulv(cross(z, q), -s), mulv(z, dot(z, q) * c)));
        let q = sub(p, self.a);
        D2 { p: addv(self.a, rot(q)), u: rot(d1), v: drot(q), uu: rot(d2), uv: drot(d1), vv: ddrot(q) }
    }
    fn dom(&self) -> [f64; 4] { let (a, b) = self.c.range(); [a, b, 0.0, 2.0 * PI] }
    fn per_v(&self) -> bool { true }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}
