//! Vectors, the parametric surface interface and the analytic surfaces used by the spikes.
use std::f64::consts::PI;

pub type V = [f64; 3];
pub fn sub(a: V, b: V) -> V { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
pub fn addv(a: V, b: V) -> V { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
pub fn mulv(a: V, s: f64) -> V { [a[0] * s, a[1] * s, a[2] * s] }
pub fn dot(a: V, b: V) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
pub fn cross(a: V, b: V) -> V { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
pub fn norm(a: V) -> f64 { dot(a, a).sqrt() }
pub fn unit(a: V) -> V { mulv(a, 1.0 / norm(a)) }

/// Point with first and second partial derivatives.
pub struct D2 { pub p: V, pub u: V, pub v: V, pub uu: V, pub uv: V, pub vv: V }

pub trait Surf: Send + Sync {
    fn ev(&self, u: f64, v: f64) -> D2;
    /// [u0, u1, v0, v1]
    fn dom(&self) -> [f64; 4];
    fn per_u(&self) -> bool { false }
    fn per_v(&self) -> bool { false }
    /// v-range ends are parametric poles (degenerate points), not edges: crossing one continues on the other side
    fn poles_v(&self) -> bool { false }
    /// unsigned distance from the implicit equation (analytic surfaces only)
    fn implicit(&self, _p: V) -> Option<f64> { None }
    /// parameters of a point on (or near) the surface
    fn uv_of(&self, p: V) -> [f64; 2];
}

pub struct Plane { pub o: V, pub a: V, pub b: V, pub ext: f64 }
impl Surf for Plane {
    fn ev(&self, u: f64, v: f64) -> D2 { D2 { p: addv(self.o, addv(mulv(self.a, u), mulv(self.b, v))), u: self.a, v: self.b, uu: [0.0; 3], uv: [0.0; 3], vv: [0.0; 3] } }
    fn dom(&self) -> [f64; 4] { [-self.ext, self.ext, -self.ext, self.ext] }
    fn implicit(&self, p: V) -> Option<f64> { let n = unit(cross(self.a, self.b)); Some(dot(sub(p, self.o), n).abs()) }
    fn uv_of(&self, p: V) -> [f64; 2] {
        let d = sub(p, self.o);
        let (aa, ab, bb) = (dot(self.a, self.a), dot(self.a, self.b), dot(self.b, self.b));
        let (r1, r2) = (dot(d, self.a), dot(d, self.b));
        let det = aa * bb - ab * ab;
        [(bb * r1 - ab * r2) / det, (aa * r2 - ab * r1) / det]
    }
}

/// cylinder about the origin: p = r cos u x + r sin u y + v z
pub struct Cyl { pub r: f64, pub h: f64, pub x: V, pub y: V, pub z: V }
impl Surf for Cyl {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (s, c) = u.sin_cos();
        let radial = addv(mulv(self.x, self.r * c), mulv(self.y, self.r * s));
        D2 { p: addv(radial, mulv(self.z, v)), u: addv(mulv(self.x, -self.r * s), mulv(self.y, self.r * c)), v: self.z, uu: mulv(radial, -1.0), uv: [0.0; 3], vv: [0.0; 3] }
    }
    fn dom(&self) -> [f64; 4] { [0.0, 2.0 * PI, -self.h, self.h] }
    fn per_u(&self) -> bool { true }
    fn implicit(&self, p: V) -> Option<f64> { let q = [dot(p, self.x), dot(p, self.y)]; Some(((q[0] * q[0] + q[1] * q[1]).sqrt() - self.r).abs()) }
    fn uv_of(&self, p: V) -> [f64; 2] { [dot(p, self.y).atan2(dot(p, self.x)).rem_euclid(2.0 * PI), dot(p, self.z)] }
}

pub struct Sphere { pub r: f64 }
impl Surf for Sphere {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (su, cu, sv, cv, r) = (u.sin(), u.cos(), v.sin(), v.cos(), self.r);
        D2 { p: [r * cv * cu, r * cv * su, r * sv], u: [-r * cv * su, r * cv * cu, 0.0], v: [-r * sv * cu, -r * sv * su, r * cv],
             uu: [-r * cv * cu, -r * cv * su, 0.0], uv: [r * sv * su, -r * sv * cu, 0.0], vv: [-r * cv * cu, -r * cv * su, -r * sv] }
    }
    fn dom(&self) -> [f64; 4] { [0.0, 2.0 * PI, -PI / 2.0, PI / 2.0] }
    fn per_u(&self) -> bool { true }
    fn poles_v(&self) -> bool { true }
    fn implicit(&self, p: V) -> Option<f64> { Some((norm(p) - self.r).abs()) }
    fn uv_of(&self, p: V) -> [f64; 2] { [p[1].atan2(p[0]).rem_euclid(2.0 * PI), (p[2] / norm(p).max(1e-300)).clamp(-1.0, 1.0).asin()] }
}

pub struct Torus { pub big: f64, pub r: f64 }
impl Surf for Torus {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (su, cu, sv, cv) = (u.sin(), u.cos(), v.sin(), v.cos());
        let w = self.big + self.r * cv;
        D2 { p: [w * cu, w * su, self.r * sv], u: [-w * su, w * cu, 0.0], v: [-self.r * sv * cu, -self.r * sv * su, self.r * cv],
             uu: [-w * cu, -w * su, 0.0], uv: [self.r * sv * su, -self.r * sv * cu, 0.0], vv: [-self.r * cv * cu, -self.r * cv * su, -self.r * sv] }
    }
    fn dom(&self) -> [f64; 4] { [0.0, 2.0 * PI, 0.0, 2.0 * PI] }
    fn per_u(&self) -> bool { true }
    fn per_v(&self) -> bool { true }
    fn implicit(&self, p: V) -> Option<f64> { let w = (p[0] * p[0] + p[1] * p[1]).sqrt() - self.big; Some(((w * w + p[2] * p[2]).sqrt() - self.r).abs()) }
    fn uv_of(&self, p: V) -> [f64; 2] { let w = (p[0] * p[0] + p[1] * p[1]).sqrt() - self.big; [p[1].atan2(p[0]).rem_euclid(2.0 * PI), p[2].atan2(w).rem_euclid(2.0 * PI)] }
}

/// cone with apex at the origin: p = v (t cos u, t sin u, 1), v in [v0, v1]
pub struct Cone { pub t: f64, pub v0: f64, pub v1: f64 }
impl Surf for Cone {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (su, cu, t) = (u.sin(), u.cos(), self.t);
        D2 { p: [v * t * cu, v * t * su, v], u: [-v * t * su, v * t * cu, 0.0], v: [t * cu, t * su, 1.0], uu: [-v * t * cu, -v * t * su, 0.0], uv: [-t * su, t * cu, 0.0], vv: [0.0; 3] }
    }
    fn dom(&self) -> [f64; 4] { [0.0, 2.0 * PI, self.v0, self.v1] }
    fn per_u(&self) -> bool { true }
    fn implicit(&self, p: V) -> Option<f64> { let a = self.t.atan(); Some(((p[0] * p[0] + p[1] * p[1]).sqrt() * a.cos() - p[2] * a.sin()).abs()) }
    fn uv_of(&self, p: V) -> [f64; 2] { [p[1].atan2(p[0]).rem_euclid(2.0 * PI), p[2]] }
}

/// bicubic Bezier patch on [0,1]^2
pub struct Bez { pub cp: [[f64; 3]; 16] }
pub fn bern(t: f64) -> [f64; 4] { let s = 1.0 - t; [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t] }
pub fn dbern(t: f64) -> [f64; 4] { let s = 1.0 - t; [-3.0 * s * s, 3.0 * s * s - 6.0 * s * t, 6.0 * s * t - 3.0 * t * t, 3.0 * t * t] }
fn ddbern(t: f64) -> [f64; 4] { let s = 1.0 - t; [6.0 * s, -12.0 * s + 6.0 * t, 6.0 * s - 12.0 * t, 6.0 * t] }
impl Surf for Bez {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (b0u, b1u, b2u, b0v, b1v, b2v) = (bern(u), dbern(u), ddbern(u), bern(v), dbern(v), ddbern(v));
        let mut d = D2 { p: [0.0; 3], u: [0.0; 3], v: [0.0; 3], uu: [0.0; 3], uv: [0.0; 3], vv: [0.0; 3] };
        for i in 0..4 { for j in 0..4 {
            let c = self.cp[i * 4 + j];
            for k in 0..3 {
                d.p[k] += b0u[i] * b0v[j] * c[k]; d.u[k] += b1u[i] * b0v[j] * c[k]; d.v[k] += b0u[i] * b1v[j] * c[k];
                d.uu[k] += b2u[i] * b0v[j] * c[k]; d.uv[k] += b1u[i] * b1v[j] * c[k]; d.vv[k] += b0u[i] * b2v[j] * c[k];
            }
        } }
        d
    }
    fn dom(&self) -> [f64; 4] { [0.0, 1.0, 0.0, 1.0] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// rigid placement of any surface: p' = R p + t (R given by rows)
pub struct Moved { pub s: Box<dyn Surf>, pub r: [V; 3], pub t: V }
impl Moved {
    fn rot(&self, a: V) -> V { [dot(self.r[0], a), dot(self.r[1], a), dot(self.r[2], a)] }
    fn inv(&self, p: V) -> V { let q = sub(p, self.t); [self.r[0][0] * q[0] + self.r[1][0] * q[1] + self.r[2][0] * q[2], self.r[0][1] * q[0] + self.r[1][1] * q[1] + self.r[2][1] * q[2], self.r[0][2] * q[0] + self.r[1][2] * q[1] + self.r[2][2] * q[2]] }
}
impl Surf for Moved {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let d = self.s.ev(u, v);
        D2 { p: addv(self.rot(d.p), self.t), u: self.rot(d.u), v: self.rot(d.v), uu: self.rot(d.uu), uv: self.rot(d.uv), vv: self.rot(d.vv) }
    }
    fn dom(&self) -> [f64; 4] { self.s.dom() }
    fn per_u(&self) -> bool { self.s.per_u() }
    fn per_v(&self) -> bool { self.s.per_v() }
    fn poles_v(&self) -> bool { self.s.poles_v() }
    fn implicit(&self, p: V) -> Option<f64> { self.s.implicit(self.inv(p)) }
    fn uv_of(&self, p: V) -> [f64; 2] { self.s.uv_of(self.inv(p)) }
}
/// rows of the rotation whose columns are the given axes
pub fn frame(x: V, y: V, z: V) -> [V; 3] { [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]] }

/// closest-point parameters for surfaces without a closed-form inverse: grid search + Newton
pub fn uv_newton(s: &dyn Surf, p: V) -> [f64; 2] {
    let d = s.dom();
    let n = 24;
    let mut best = (d[0], d[2], f64::MAX);
    for i in 0..=n { for j in 0..=n {
        let (u, v) = (d[0] + (d[1] - d[0]) * i as f64 / n as f64, d[2] + (d[3] - d[2]) * j as f64 / n as f64);
        let q = norm(sub(s.ev(u, v).p, p));
        if q < best.2 { best = (u, v, q); }
    } }
    uv_refine(s, p, [best.0, best.1])
}
/// Newton on the closest-point conditions from a starting guess, clamped to the domain
pub fn uv_refine(s: &dyn Surf, p: V, g: [f64; 2]) -> [f64; 2] {
    let d = s.dom();
    let (mut u, mut v) = (g[0], g[1]);
    for _ in 0..60 {
        let e = s.ev(u, v);
        let r = sub(e.p, p);
        let (f1, f2) = (dot(r, e.u), dot(r, e.v));
        let (a11, a12, a22) = (dot(e.u, e.u) + dot(r, e.uu), dot(e.u, e.v) + dot(r, e.uv), dot(e.v, e.v) + dot(r, e.vv));
        let det = a11 * a22 - a12 * a12;
        if !(det.abs() > 1e-300) { break; }
        let (du, dv) = ((a22 * f1 - a12 * f2) / det, (a11 * f2 - a12 * f1) / det);
        // periodic parameters wrap, bounded ones clamp to the domain
        let lim = |x: f64, lo: f64, hi: f64, per: bool| if per { lo + (x - lo).rem_euclid(hi - lo) } else { x.clamp(lo, hi) };
        let (nu, nv) = (lim(u - du, d[0], d[1], s.per_u()), lim(v - dv, d[2], d[3], s.per_v()));
        let done = (nu - u).abs() + (nv - v).abs() < 1e-16;
        u = nu; v = nv;
        if done { break; }
    }
    [u, v]
}
pub fn dist(s: &dyn Surf, p: V) -> f64 { s.implicit(p).unwrap_or_else(|| { let uv = s.uv_of(p); norm(sub(s.ev(uv[0], uv[1]).p, p)) }) }
