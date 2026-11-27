//! CATIA CGM exact B-rep, read from the flat IR written by cadcore (`cadcli cgm PART OUT.cgmir`,
//! format in cadcore/src/formats/catia/cgm_model.hpp). Surfaces and p-curves keep CGM's own
//! parametrizations; faces are bounded by p-curves, so each coedge is a parameter range on one.
use crate::geom::*;
use crate::nurbs::BspS;
use std::sync::Arc;

/// C2 piecewise quintic Hermite (C1 cubic without second derivatives) of `dim` values per knot, row-major
pub struct Hermite { pub k: Vec<f64>, pub dim: usize, pub p: Vec<f64>, pub d1: Vec<f64>, pub d2: Option<Vec<f64>> }
impl Hermite {
    pub fn eval(&self, t: f64, out: &mut [f64]) {
        let (k, n, dim) = (&self.k, self.k.len(), self.dim);
        let mut i = 0;
        while i + 2 < n && t > k[i + 1] { i += 1; }
        let h = k[i + 1] - k[i];
        let s = if h != 0.0 { (t - k[i]) / h } else { 0.0 };
        let (s2, s3) = (s * s, s * s * s);
        let (s4, s5) = (s3 * s, s3 * s2);
        let (a, b) = (i * dim, (i + 1) * dim);
        for c in 0..dim {
            out[c] = match &self.d2 {
                Some(d2) => {
                    let h0 = 1.0 - 10.0 * s3 + 15.0 * s4 - 6.0 * s5;
                    let h1 = s - 6.0 * s3 + 8.0 * s4 - 3.0 * s5;
                    let h2 = 0.5 * s2 - 1.5 * s3 + 1.5 * s4 - 0.5 * s5;
                    let h5 = 10.0 * s3 - 15.0 * s4 + 6.0 * s5;
                    let h4 = -4.0 * s3 + 7.0 * s4 - 3.0 * s5;
                    let h3 = 0.5 * s3 - s4 + 0.5 * s5;
                    h0 * self.p[a + c] + h * h1 * self.d1[a + c] + h * h * h2 * d2[a + c] + h5 * self.p[b + c] + h * h4 * self.d1[b + c] + h * h * h3 * d2[b + c]
                }
                None => {
                    let (a0, a1, b0, b1) = (2.0 * s3 - 3.0 * s2 + 1.0, s3 - 2.0 * s2 + s, -2.0 * s3 + 3.0 * s2, s3 - s2);
                    a0 * self.p[a + c] + h * a1 * self.d1[a + c] + b0 * self.p[b + c] + h * b1 * self.d1[b + c]
                }
            };
        }
    }
    pub fn range(&self) -> (f64, f64) { (self.k[0], self.k[self.k.len() - 1]) }
    /// `n | knots[n] | P[n][dim] | D1[n][dim] | D2[n][dim]`; returns it and the values used
    fn read(p: &[f64], dim: usize) -> Option<(Hermite, usize)> {
        let n = *p.first()? as usize;
        let end = 1 + n + 3 * n * dim;
        if n < 2 || p.len() < end { return None; }
        let a = |j: usize| p[1 + n + j * n * dim..1 + n + (j + 1) * n * dim].to_vec();
        Some((Hermite { k: p[1..1 + n].to_vec(), dim, p: a(0), d1: a(1), d2: Some(a(2)) }, end))
    }
}

/// B-spline of `dim` coordinates (expanded knots), poles homogeneous (dim + 1 values, weight last); evaluated past
/// its ends by extending the end spans
pub struct Bspline { pub deg: usize, pub k: Vec<f64>, pub dim: usize, pub cp: Vec<f64> }
impl Bspline {
    pub fn eval(&self, t: f64, out: &mut [f64]) {
        let w = self.dim + 1;
        let n = self.cp.len() / w - 1;
        // the end spans' polynomials continue past the ends
        let i = crate::nurbs::span(n, self.deg, t, &self.k);
        let b = crate::nurbs::ders(i, t, self.deg, &self.k)[0];
        let mut h = [0.0; 4];
        for j in 0..=self.deg { for c in 0..w { h[c] += b[j] * self.cp[(i - self.deg + j) * w + c]; } }
        for c in 0..self.dim { out[c] = h[c] / h[self.dim]; }
    }
    pub fn range(&self) -> (f64, f64) { (self.k[self.deg], self.k[self.k.len() - 1 - self.deg]) }
    /// `degree nk | knots[nk] | np rational | poles[np][dim] | weights[np if rational]`; None also for decreasing
    /// knots or weights that are not positive
    fn read(p: &[f64], dim: usize) -> Option<(Bspline, usize)> {
        let (deg, nk) = (*p.first()? as usize, *p.get(1)? as usize);
        let (np, rat) = (*p.get(2 + nk)? as usize, *p.get(3 + nk)? != 0.0);
        let poles = 4 + nk;
        let end = poles + np * dim + if rat { np } else { 0 };
        if deg == 0 || nk != np + deg + 1 || p.len() < end || knots_ok(&p[2..2 + nk]).is_err() { return None; }
        if rat && p[poles + np * dim..end].iter().any(|&w| w <= 0.0) { return None; }
        let cp = (0..np).flat_map(|i| {
            let w = if rat { p[poles + np * dim + i] } else { 1.0 };
            (0..dim).map(move |c| p[poles + i * dim + c] * w).chain([w])
        }).collect();
        Some((Bspline { deg, k: p[2..2 + nk].to_vec(), dim, cp }, end))
    }
}

/// expression of the curve parameter (CGM class 155)
pub enum Expr {
    /// v0 + v1 t + v2 t^2
    Poly([f64; 3]),
    /// v3 + v0 cos(v1 - v2 t)
    Cos([f64; 4]),
    Quot(Box<Expr>, Box<Expr>),
    /// sum v_i k_i + v_n
    Lin(Vec<f64>, Vec<Expr>),
    /// v4 + v0 atan(v1 cos(v2 t + v3))
    Atan([f64; 5]),
}
impl Expr {
    pub fn at(&self, t: f64) -> f64 {
        match self {
            Expr::Poly(v) => v[0] + t * (v[1] + t * v[2]),
            Expr::Cos(v) => v[3] + v[0] * (v[1] - v[2] * t).cos(),
            Expr::Quot(a, b) => a.at(t) / b.at(t),
            Expr::Lin(v, k) => k.iter().zip(v).map(|(e, c)| c * e.at(t)).sum::<f64>() + v[k.len()],
            Expr::Atan(v) => v[4] + v[0] * (v[1] * (v[2] * t + v[3]).cos()).atan(),
        }
    }
    /// `code nvals vals[nvals] nkids node[nkids]` at p[*i]
    fn read(p: &[f64], i: &mut usize) -> Option<Expr> {
        let mut next = || { let x = *p.get(*i)?; *i += 1; Some(x) };
        let (code, nv) = (next()? as u32, next()? as usize);
        let vals: Vec<f64> = (0..nv).map(|_| next()).collect::<Option<_>>()?;
        let nk = next()? as usize;
        let kids: Vec<Expr> = (0..nk).map(|_| Expr::read(p, i)).collect::<Option<_>>()?;
        let arr = |n: usize| if nv == n && nk == 0 { Some(vals.clone()) } else { None };
        Some(match code {
            1 => Expr::Poly(arr(3)?.try_into().ok()?),
            2 => Expr::Cos(arr(4)?.try_into().ok()?),
            3 if nv == 0 && nk == 2 => { let mut k = kids.into_iter(); Expr::Quot(Box::new(k.next()?), Box::new(k.next()?)) }
            4 if nv == nk + 1 => Expr::Lin(vals, kids),
            7 => Expr::Atan(arr(5)?.try_into().ok()?),
            _ => return None,
        })
    }
}

/// p-curve in a surface's (u, v)
pub enum PCurve {
    Line { o: [f64; 2], d: [f64; 2] },
    /// angle = t / (r s)
    Circle { c: [f64; 2], r: f64, s: f64 },
    Hermite(Hermite),
    Bspline(Bspline),
    Expr(Expr, Expr),
}
impl PCurve {
    pub fn uv(&self, t: f64) -> [f64; 2] {
        let mut o = [0.0; 2];
        match self {
            PCurve::Line { o, d } => return [o[0] + t * d[0], o[1] + t * d[1]],
            PCurve::Circle { c, r, s } => { let a = t / (r * s); return [c[0] + r * a.cos(), c[1] + r * a.sin()]; }
            PCurve::Hermite(h) => h.eval(t, &mut o),
            PCurve::Bspline(b) => b.eval(t, &mut o),
            PCurve::Expr(u, v) => return [u.at(t), v.at(t)],
        }
        o
    }
}

/// a 3D curve: a p-curve on a surface, a line or a B-spline
pub enum Curve3 { On(Arc<PCurve>, Arc<dyn Surf>), Line(V, V), Bspline(Bspline) }
impl Curve3 {
    /// parameter range of a curve with finite data
    pub fn range(&self) -> Option<(f64, f64)> { if let Curve3::Bspline(b) = self { Some(b.range()) } else { None } }
    /// surface domain with this curve's range in u
    fn dom_u(&self) -> [f64; 4] { let (a, b) = self.range().unwrap_or((-1e9, 1e9)); [a, b, -1e9, 1e9] }
    pub fn at(&self, t: f64) -> V {
        match self {
            Curve3::On(c, s) => { let uv = c.uv(t); s.ev(uv[0], uv[1]).p }
            Curve3::Line(o, d) => addv(*o, mulv(*d, t)),
            Curve3::Bspline(b) => { let mut o = [0.0; 3]; b.eval(t, &mut o); o }
        }
    }
}

/// surface derivatives by central differences of `f`
fn diff(f: impl Fn(f64, f64) -> V, u: f64, v: f64) -> D2 {
    let (hu, hv) = (1e-6 * (1.0 + u.abs()), 1e-6 * (1.0 + v.abs()));
    D2 { p: f(u, v), u: mulv(sub(f(u + hu, v), f(u - hu, v)), 0.5 / hu), v: mulv(sub(f(u, v + hv), f(u, v - hv)), 0.5 / hv), uu: ZERO, uv: ZERO, vv: ZERO }
}

/// p = o + u U + v V
pub struct CgmPlane { pub o: V, pub u: V, pub v: V }
/// p = o + r (cos a X + sin a Y) + v Z, a = u / su, Z = X x Y
pub struct CgmCyl { pub o: V, pub x: V, pub y: V, pub z: V, pub r: f64, pub su: f64 }
/// p = o + cos b (cos a A + sin a B) + sin b C, a = u / k, b = v / k (A, B, C scaled by the radius)
pub struct CgmSphere { pub o: V, pub a: V, pub b: V, pub c: V, pub k: f64 }
/// p = o + (r0 + v cos a)(cos t X + sin t Y) + v sin a Z, t = u / su
pub struct CgmCone { pub o: V, pub x: V, pub y: V, pub z: V, pub a: f64, pub r0: f64, pub su: f64 }
/// p = o + (R + r cos b)(cos a X + sin a Y) + r sin b Z, a = u / R, b = v / r
pub struct CgmTorus { pub o: V, pub x: V, pub y: V, pub z: V, pub big: f64, pub r: f64 }

const ZERO: V = [0.0; 3];
fn lin(o: V, terms: &[(V, f64)]) -> V { terms.iter().fold(o, |acc, (d, s)| addv(acc, mulv(*d, *s))) }

impl Surf for CgmPlane {
    fn ev(&self, u: f64, v: f64) -> D2 { D2 { p: lin(self.o, &[(self.u, u), (self.v, v)]), u: self.u, v: self.v, uu: ZERO, uv: ZERO, vv: ZERO } }
    fn dom(&self) -> [f64; 4] { [-1e9, 1e9, -1e9, 1e9] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}
impl Surf for CgmCyl {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let a = u / self.su;
        let (s, c) = a.sin_cos();
        let k = self.r / self.su;
        D2 { p: lin(self.o, &[(self.x, self.r * c), (self.y, self.r * s), (self.z, v)]), u: lin(ZERO, &[(self.x, -k * s), (self.y, k * c)]), v: self.z, uu: ZERO, uv: ZERO, vv: ZERO }
    }
    fn dom(&self) -> [f64; 4] { [-1e9, 1e9, -1e9, 1e9] }
    fn periods(&self) -> [f64; 2] { [std::f64::consts::TAU * self.su, 0.0] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}
impl Surf for CgmSphere {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (a, b) = (u / self.k, v / self.k);
        let ((sa, ca), (sb, cb)) = (a.sin_cos(), b.sin_cos());
        let rim = lin(ZERO, &[(self.a, ca), (self.b, sa)]);
        D2 { p: addv(self.o, addv(mulv(rim, cb), mulv(self.c, sb))), u: mulv(lin(ZERO, &[(self.a, -sa), (self.b, ca)]), cb / self.k), v: mulv(addv(mulv(rim, -sb), mulv(self.c, cb)), 1.0 / self.k), uu: ZERO, uv: ZERO, vv: ZERO }
    }
    fn dom(&self) -> [f64; 4] { [-1e9, 1e9, -1e9, 1e9] }
    fn periods(&self) -> [f64; 2] { [std::f64::consts::TAU * self.k, 0.0] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}
impl Surf for CgmCone {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let t = u / self.su;
        let (st, ct) = t.sin_cos();
        let (sa, ca) = self.a.sin_cos();
        let dir = lin(ZERO, &[(self.x, ct), (self.y, st)]);
        let w = self.r0 + v * ca;
        D2 { p: addv(self.o, addv(mulv(dir, w), mulv(self.z, v * sa))), u: mulv(lin(ZERO, &[(self.x, -st), (self.y, ct)]), w / self.su), v: addv(mulv(dir, ca), mulv(self.z, sa)), uu: ZERO, uv: ZERO, vv: ZERO }
    }
    fn dom(&self) -> [f64; 4] { [-1e9, 1e9, -1e9, 1e9] }
    fn periods(&self) -> [f64; 2] { [std::f64::consts::TAU * self.su, 0.0] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}
impl Surf for CgmTorus {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let (a, b) = (u / self.big, v / self.r);
        let ((sa, ca), (sb, cb)) = (a.sin_cos(), b.sin_cos());
        let dir = lin(ZERO, &[(self.x, ca), (self.y, sa)]);
        let w = self.big + self.r * cb;
        D2 { p: addv(self.o, addv(mulv(dir, w), mulv(self.z, self.r * sb))), u: mulv(lin(ZERO, &[(self.x, -sa), (self.y, ca)]), w / self.big), v: addv(mulv(dir, -sb), mulv(self.z, cb)), uu: ZERO, uv: ZERO, vv: ZERO }
    }
    fn dom(&self) -> [f64; 4] { [-1e9, 1e9, -1e9, 1e9] }
    fn periods(&self) -> [f64; 2] { [std::f64::consts::TAU * self.big, std::f64::consts::TAU * self.r] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// p = C(u) + v d
pub struct CgmExtrusion { pub c: Curve3, pub d: V }
impl Surf for CgmExtrusion {
    fn ev(&self, u: f64, v: f64) -> D2 { let mut e = diff(|u, v| addv(self.c.at(u), mulv(self.d, v)), u, v); e.v = self.d; e }
    fn dom(&self) -> [f64; 4] { self.c.dom_u() }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// p = o + Rot_Z(v / k) (C(u) - o), turning X toward Y
pub struct CgmRevolution { pub o: V, pub x: V, pub y: V, pub z: V, pub k: f64, pub c: Curve3 }
impl Surf for CgmRevolution {
    fn ev(&self, u: f64, v: f64) -> D2 {
        diff(|u, v| {
            let w = sub(self.c.at(u), self.o);
            let (x, y, z) = (dot(w, self.x), dot(w, self.y), dot(w, self.z));
            let (s, c) = (v / self.k).sin_cos();
            lin(self.o, &[(self.x, x * c - y * s), (self.y, x * s + y * c), (self.z, z)])
        }, u, v)
    }
    fn dom(&self) -> [f64; 4] { self.c.dom_u() }
    fn periods(&self) -> [f64; 2] { [0.0, std::f64::consts::TAU * self.k] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// rolling-ball fillet: section(u) = contact points p1, p2, centre c, arc angle th;
/// p = c + (sin(th - f) (p1 - c) + sin f (p2 - c)) / sin th, f = th (v - v0) / (v1 - v0)
pub struct CgmFillet { pub h: Hermite, pub v0: f64, pub v1: f64 }
impl Surf for CgmFillet {
    fn ev(&self, u: f64, v: f64) -> D2 {
        diff(|u, v| {
            let mut s = [0.0; 10];
            self.h.eval(u, &mut s);
            let (p1, p2, c, th) = (v3(&s, 0), v3(&s, 3), v3(&s, 6), s[9]);
            let f = th * (v - self.v0) / (self.v1 - self.v0);
            lin(c, &[(sub(p1, c), (th - f).sin() / th.sin()), (sub(p2, c), f.sin() / th.sin())])
        }, u, v)
    }
    fn dom(&self) -> [f64; 4] { let (a, b) = self.h.range(); [a, b, self.v0.min(self.v1), self.v0.max(self.v1)] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// p = A(u) + v (B(u) - A(u))
pub struct CgmRuled { pub h: Hermite }
impl Surf for CgmRuled {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let ab = |u: f64| { let mut s = [0.0; 6]; self.h.eval(u, &mut s); (v3(&s, 0), v3(&s, 3)) };
        let mut e = diff(|u, v| { let (a, b) = ab(u); addv(a, mulv(sub(b, a), v)) }, u, v);
        let (a, b) = ab(u);
        e.v = sub(b, a);
        e
    }
    fn dom(&self) -> [f64; 4] { let (a, b) = self.h.range(); [a, b, -1e9, 1e9] }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

/// p = base(u, v) + d unit(Su x Sv)
pub struct CgmOffset { pub base: Arc<dyn Surf>, pub d: f64 }
impl Surf for CgmOffset {
    fn ev(&self, u: f64, v: f64) -> D2 {
        let n = |u: f64, v: f64| { let e = self.base.ev(u, v); (e.p, e.u, e.v, unit(cross(e.u, e.v))) };
        let (p, su, sv, nn) = n(u, v);
        let h = 1e-6 * (1.0 + u.abs() + v.abs());
        let (nu, nv) = (sub(n(u + h, v).3, n(u - h, v).3), sub(n(u, v + h).3, n(u, v - h).3));
        D2 { p: addv(p, mulv(nn, self.d)), u: addv(su, mulv(nu, self.d / (2.0 * h))), v: addv(sv, mulv(nv, self.d / (2.0 * h))), uu: ZERO, uv: ZERO, vv: ZERO }
    }
    fn dom(&self) -> [f64; 4] { self.base.dom() }
    fn periods(&self) -> [f64; 2] { self.base.periods() }
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

pub struct Coedge { pub pcurve: usize, pub edge: usize, pub sense: i32, pub t: [f64; 2] }
pub struct Face { pub surface: usize, pub sense: i32, pub loops: Vec<std::ops::Range<usize>>, pub id: u32 }
/// `solid`: bounds a volume; otherwise a sheet (open surfaces, reference planes)
pub struct Shell { pub id: u32, pub solid: bool, pub faces: Vec<(usize, i32)> }
/// a CGM body (versions of a body share unchanged shells); `flags` bit 0: named by the stream
/// header's roots, bit 1 + j: listed in list j of the part's CATSeeBodyCont
pub struct Body { pub id: u32, pub flags: u32, pub shells: Vec<usize> }
pub struct Skip { pub id: u32, pub kind: u32, pub why: String }

/// the B-rep of one part
pub struct Model {
    pub surfaces: Vec<Arc<dyn Surf>>,
    /// IR kind per surface and per p-curve
    pub surface_kinds: Vec<u32>,
    pub pcurve_kinds: Vec<u32>,
    /// (surface, curve)
    pub pcurves: Vec<(usize, Arc<PCurve>)>,
    pub vertices: usize,
    /// vertex indices (start, end)
    pub edges: Vec<[usize; 2]>,
    pub coedges: Vec<Coedge>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub bodies: Vec<Body>,
    pub skipped: Vec<Skip>,
    /// why a surface or p-curve record could not be used (its faces fail with the reason); None when valid
    pub bad_surface: Vec<Option<String>>,
    pub bad_pcurve: Vec<Option<String>>,
}

struct Rd<'a> { b: &'a [u8], k: usize }
impl<'a> Rd<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        if self.k + n > self.b.len() { return Err(format!("IR ends early at byte {}", self.k)); }
        let s = &self.b[self.k..self.k + n];
        self.k += n;
        Ok(s)
    }
    fn u32(&mut self) -> Result<u32, String> { Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap())) }
    fn i32(&mut self) -> Result<i32, String> { Ok(self.u32()? as i32) }
    fn f64(&mut self) -> Result<f64, String> { Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap())) }
    fn n(&mut self) -> Result<usize, String> { Ok(self.u32()? as usize) }
    fn f64s(&mut self, n: usize) -> Result<Vec<f64>, String> { (0..n).map(|_| self.f64()).collect() }
}

/// knots non-decreasing, with a non-empty range
fn knots_ok(k: &[f64]) -> Result<(), String> {
    if k.windows(2).all(|w| w[0] <= w[1]) && k.first() < k.last() { Ok(()) } else { Err("knots not increasing".into()) }
}

/// every parameter finite
fn valid(p: &[f64]) -> Result<(), String> { if p.iter().all(|x| x.is_finite()) { Ok(()) } else { Err("non-finite parameter".into()) } }

fn v3(p: &[f64], i: usize) -> V { [p[i], p[i + 1], p[i + 2]] }

fn surface(kind: u32, p: &[f64]) -> Result<Arc<dyn Surf>, String> {
    let need = |n: usize| if p.len() < n { Err(format!("surface kind {kind}: {} parameters", p.len())) } else { Ok(()) };
    Ok(match kind {
        1 => { need(9)?; Arc::new(CgmPlane { o: v3(p, 0), u: v3(p, 3), v: v3(p, 6) }) }
        2 => { need(11)?; let (x, y) = (v3(p, 3), v3(p, 6)); Arc::new(CgmCyl { o: v3(p, 0), x, y, z: cross(x, y), r: p[9], su: p[10] }) }
        3 => { need(13)?; Arc::new(CgmSphere { o: v3(p, 0), a: v3(p, 3), b: v3(p, 6), c: v3(p, 9), k: p[12] }) }
        4 => { need(14)?; Arc::new(CgmTorus { o: v3(p, 0), x: v3(p, 3), y: v3(p, 6), z: v3(p, 9), big: p[12], r: p[13] }) }
        6 => { need(15)?; Arc::new(CgmCone { o: v3(p, 0), x: v3(p, 3), y: v3(p, 6), z: v3(p, 9), a: p[12], r0: p[13], su: p[14] }) }
        5 => {
            need(5)?;
            let (du, dv, nu, nv, rat) = (p[0] as usize, p[1] as usize, p[2] as usize, p[3] as usize, p[4] != 0.0);
            let (ku, kv) = (nu + du + 1, nv + dv + 1);
            need(5 + ku + kv + 3 * nu * nv + if rat { nu * nv } else { 0 })?;
            let k0 = 5;
            let poles = k0 + ku + kv;
            let wts = poles + 3 * nu * nv;
            knots_ok(&p[k0..k0 + ku])?;
            knots_ok(&p[k0 + ku..k0 + ku + kv])?;
            if rat && p[wts..wts + nu * nv].iter().any(|&w| w <= 0.0) { return Err("NURBS weight not positive".into()); }
            let cp = (0..nu).map(|i| (0..nv).map(|j| {
                let q = v3(p, poles + 3 * (i * nv + j));
                let w = if rat { p[wts + i * nv + j] } else { 1.0 };
                [q[0] * w, q[1] * w, q[2] * w, w]
            }).collect()).collect();
            Arc::new(BspS { pu: du, pv: dv, ku: p[k0..k0 + ku].to_vec(), kv: p[k0 + ku..k0 + ku + kv].to_vec(), cp })
        }
        _ => return Err(format!("surface kind {kind}")),
    })
}

fn pcurve(kind: u32, p: &[f64]) -> Result<PCurve, String> {
    let bad = || format!("p-curve kind {kind}: {} parameters", p.len());
    Ok(match kind {
        1 if p.len() == 4 => PCurve::Line { o: [p[0], p[1]], d: [p[2], p[3]] },
        2 if p.len() == 4 => PCurve::Circle { c: [p[0], p[1]], r: p[2], s: p[3] },
        // knots[n], u[n], v[n], du[n], dv[n], has_d2, (d2u[n], d2v[n])
        3 => {
            let n = *p.first().ok_or_else(bad)? as usize;
            if n < 2 || p.len() < 2 + 5 * n { return Err(bad()); }
            let has2 = p[1 + 5 * n] != 0.0;
            if has2 && p.len() < 2 + 7 * n { return Err(bad()); }
            let rows = |a: usize, b: usize| (0..n).flat_map(|i| [p[a + i], p[b + i]]).collect::<Vec<f64>>();
            let d2 = has2.then(|| rows(2 + 5 * n, 2 + 6 * n));
            PCurve::Hermite(Hermite { k: p[1..1 + n].to_vec(), dim: 2, p: rows(1 + n, 1 + 2 * n), d1: rows(1 + 3 * n, 1 + 4 * n), d2 })
        }
        4 => match Bspline::read(p, 2) { Some((b, end)) if end == p.len() => PCurve::Bspline(b), _ => return Err(bad()) },
        // t0 t1 a b, then the u and v expressions
        5 => {
            let mut i = 4;
            let (u, v) = (Expr::read(p, &mut i).ok_or_else(bad)?, Expr::read(p, &mut i).ok_or_else(bad)?);
            if i != p.len() { return Err(bad()); }
            PCurve::Expr(u, v)
        }
        _ => return Err(bad()),
    })
}

impl Model {
    pub fn read(b: &[u8]) -> Result<Model, String> {
        if b.len() < 8 || &b[..8] != b"CGMIR\0\0\x02" { return Err("not a CGM IR (version 2)".into()); }
        let mut r = Rd { b, k: 8 };
        let c: Vec<usize> = (0..10).map(|_| r.n()).collect::<Result<_, _>>()?;
        let mut m = Model { surfaces: vec![], surface_kinds: vec![], pcurve_kinds: vec![], pcurves: vec![], vertices: c[2], edges: vec![], coedges: vec![], faces: vec![], shells: vec![], bodies: vec![], skipped: vec![], bad_surface: vec![], bad_pcurve: vec![] };
        let raw: Vec<(u32, Vec<f64>)> = (0..c[0]).map(|_| { let (kind, _id, n) = (r.u32()?, r.u32()?, r.n()?); Ok((kind, r.f64s(n)?)) }).collect::<Result<_, String>>()?;
        for _ in 0..c[1] {
            let (kind, s, _id, n) = (r.u32()?, r.n()?, r.u32()?, r.n()?);
            let p = r.f64s(n)?;
            if s >= raw.len() { return Err("p-curve surface out of range".into()); }
            let pc = valid(&p).and_then(|_| pcurve(kind, &p));
            m.bad_pcurve.push(pc.as_ref().err().cloned());
            m.pcurves.push((s, Arc::new(pc.unwrap_or(PCurve::Line { o: [0.0; 2], d: [0.0; 2] }))));
            m.pcurve_kinds.push(kind);
        }
        m.surface_kinds = raw.iter().map(|x| x.0).collect();
        // surfaces in order; an extrusion's profile lies on an earlier surface
        for (kind, p) in &raw {
            let s = valid(p).and_then(|_| Ok(match kind {
                7 => {
                    if p.len() != 4 { return Err("extrusion parameters".into()); }
                    Arc::new(CgmExtrusion { c: m.curve3(&[1.0, p[0]])?.0, d: v3(p, 1) }) as Arc<dyn Surf>
                }
                8 => {
                    if p.len() != 2 { return Err("offset parameters".into()); }
                    let base = m.surfaces.get(p[0] as usize).ok_or("offset base not before it")?.clone();
                    Arc::new(CgmOffset { base, d: p[1] })
                }
                9 => {
                    if p.len() < 14 { return Err("revolution parameters".into()); }
                    let (c, n) = m.curve3(&p[13..])?;
                    if 13 + n != p.len() { return Err("revolution parameters".into()); }
                    Arc::new(CgmRevolution { o: v3(p, 0), x: v3(p, 3), y: v3(p, 6), z: v3(p, 9), k: p[12], c })
                }
                10 | 11 => {
                    let dim = if *kind == 10 { 10 } else { 6 };
                    let (h, n) = Hermite::read(p, dim).ok_or(format!("surface kind {kind} parameters"))?;
                    if n + 4 != p.len() { return Err(format!("surface kind {kind} parameters")); }
                    if *kind == 10 { Arc::new(CgmFillet { h, v0: p[n + 2], v1: p[n + 3] }) as Arc<dyn Surf> } else { Arc::new(CgmRuled { h }) }
                }
                12 => {
                    if p.len() < 4 { return Err("extrusion parameters".into()); }
                    let (c, n) = m.curve3(&p[3..])?;
                    if 3 + n != p.len() { return Err("extrusion parameters".into()); }
                    Arc::new(CgmExtrusion { c, d: v3(p, 0) })
                }
                _ => surface(*kind, p)?,
            }));
            m.bad_surface.push(s.as_ref().err().cloned());
            m.surfaces.push(s.unwrap_or_else(|_| Arc::new(CgmPlane { o: ZERO, u: ZERO, v: ZERO })));
        }
        for _ in 0..c[3] { let (a, b, _id) = (r.n()?, r.n()?, r.u32()?); if a >= m.vertices || b >= m.vertices { return Err("edge vertex out of range".into()); } m.edges.push([a, b]); }
        for _ in 0..c[4] {
            let (pc, e, sense, _) = (r.n()?, r.n()?, r.i32()?, r.u32()?);
            let t = [r.f64()?, r.f64()?];
            if pc >= m.pcurves.len() || e >= m.edges.len() { return Err("coedge reference out of range".into()); }
            m.coedges.push(Coedge { pcurve: pc, edge: e, sense, t });
        }
        let loops: Vec<std::ops::Range<usize>> = (0..c[5]).map(|_| { let (a, n) = (r.n()?, r.n()?); if a + n > m.coedges.len() { Err("loop out of range".to_string()) } else { Ok(a..a + n) } }).collect::<Result<_, _>>()?;
        for _ in 0..c[6] {
            let (s, sense, l0, nl, id) = (r.n()?, r.i32()?, r.n()?, r.n()?, r.u32()?);
            if s >= m.surfaces.len() || l0 + nl > loops.len() { return Err("face reference out of range".into()); }
            m.faces.push(Face { surface: s, sense, loops: loops[l0..l0 + nl].to_vec(), id });
        }
        for _ in 0..c[7] {
            let (id, solid, n) = (r.u32()?, r.u32()? & 1 != 0, r.n()?);
            let faces = (0..n).map(|_| { let (f, s) = (r.n()?, r.i32()?); if f >= m.faces.len() { Err("shell face out of range".to_string()) } else { Ok((f, s)) } }).collect::<Result<_, _>>()?;
            m.shells.push(Shell { id, solid, faces });
        }
        for _ in 0..c[8] {
            let (id, flags, n) = (r.u32()?, r.u32()?, r.n()?);
            let shells = (0..n).map(|_| { let i = r.n()?; if i >= m.shells.len() { Err("body shell out of range".to_string()) } else { Ok(i) } }).collect::<Result<_, _>>()?;
            m.bodies.push(Body { id, flags, shells });
        }
        for _ in 0..c[9] { let (id, kind, n) = (r.u32()?, r.u32()?, r.n()?); let why = String::from_utf8_lossy(r.take(n)?).into_owned(); m.skipped.push(Skip { id, kind, why }); }
        if r.k != b.len() { return Err(format!("{} bytes after the IR", b.len() - r.k)); }
        Ok(m)
    }
    /// a 3D curve `form | data` (1: p-curve index on an earlier surface, 2: line o d, 3: B-spline); returns it and
    /// the values used
    fn curve3(&self, p: &[f64]) -> Result<(Curve3, usize), String> {
        match p.first().map(|x| *x as u32) {
            Some(1) if p.len() >= 2 => {
                let (on, c) = self.pcurves.get(p[1] as usize).ok_or("profile p-curve out of range")?;
                let s = self.surfaces.get(*on).ok_or("profile surface not before it")?.clone();
                Ok((Curve3::On(c.clone(), s), 2))
            }
            Some(2) if p.len() >= 7 => Ok((Curve3::Line(v3(p, 1), v3(p, 4)), 7)),
            Some(3) => { let (b, n) = Bspline::read(&p[1..], 3).ok_or("profile B-spline parameters")?; Ok((Curve3::Bspline(b), n + 1)) }
            _ => Err("profile form".into()),
        }
    }
    /// 3D point of coedge `c` at parameter t
    pub fn point(&self, c: &Coedge, t: f64) -> V {
        let (s, pc) = &self.pcurves[c.pcurve];
        let uv = pc.uv(t);
        self.surfaces[*s].ev(uv[0], uv[1]).p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: V, b: V) -> bool { norm(sub(a, b)) < 1e-12 }
    fn model() -> Model { Model { surfaces: vec![], surface_kinds: vec![], pcurve_kinds: vec![], pcurves: vec![], vertices: 0, edges: vec![], coedges: vec![], faces: vec![], shells: vec![], bodies: vec![], skipped: vec![], bad_surface: vec![], bad_pcurve: vec![] } }

    #[test]
    fn bspline_pcurve_matches_bezier_and_extends_past_its_ends() {
        // quadratic Bezier on knots [2, 5]: (1-s)^2 P0 + 2 s (1-s) P1 + s^2 P2, s = (t - 2) / 3
        let p = [2.0, 6.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 3.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0];
        let c = pcurve(4, &p).unwrap();
        for t in [1.7, 2.0, 3.1, 5.0, 5.2] {
            let s = (t - 2.0) / 3.0;
            let want = [2.0 * s * (1.0 - s) + 3.0 * s * s, 4.0 * s * (1.0 - s)];
            let got = c.uv(t);
            assert!((got[0] - want[0]).abs() < 1e-12 && (got[1] - want[1]).abs() < 1e-12, "t {t}: {got:?} vs {want:?}");
        }
        // rational: weights of a quarter circle
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let q = pcurve(4, &[2.0, 6.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, w, 1.0]).unwrap();
        for t in [0.0, 0.3, 0.8, 1.0] { let p = q.uv(t); assert!((p[0].hypot(p[1]) - 1.0).abs() < 1e-12); }
    }

    #[test]
    fn expression_pcurve() {
        // u = 1 + 2t + 3t^2; v = (0.5 cos(1 - 2t) + 4) / (atan(2 cos(t + 0.1)) - 2) * 3 + 7 (linear combination)
        let p = [0.0, 1.0, 0.0, 0.0,
            1.0, 3.0, 1.0, 2.0, 3.0, 0.0,
            4.0, 2.0, 3.0, 7.0, 1.0,
            3.0, 0.0, 2.0, 2.0, 4.0, 0.5, 1.0, 2.0, 4.0, 0.0, 7.0, 5.0, 1.0, 2.0, 1.0, 0.1, -2.0, 0.0];
        let c = pcurve(5, &p).unwrap();
        let t: f64 = 0.37;
        let v = 3.0 * ((0.5 * (1.0 - 2.0 * t).cos() + 4.0) / ((2.0 * (t + 0.1).cos()).atan() - 2.0)) + 7.0;
        let got = c.uv(t);
        assert!((got[0] - (1.0 + 2.0 * t + 3.0 * t * t)).abs() < 1e-12 && (got[1] - v).abs() < 1e-12, "{got:?}");
        assert!(pcurve(5, &p[..p.len() - 1]).is_err());
    }

    #[test]
    fn revolution_ruled_fillet_and_curve_extrusion() {
        let mut m = model();
        // revolving the line x = 2 (z along t) about Z through o: a cylinder of radius 2, angle v / k
        let (o, k) = ([1.0, 1.0, 0.0], 3.0);
        let rev = CgmRevolution { o, x: [1.0, 0.0, 0.0], y: [0.0, 1.0, 0.0], z: [0.0, 0.0, 1.0], k, c: m.curve3(&[2.0, 3.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap().0 };
        let e = rev.ev(0.5, 1.2);
        assert!(close(e.p, [1.0 + 2.0 * (0.4f64).cos(), 1.0 + 2.0 * (0.4f64).sin(), 0.5]), "{:?}", e.p);
        assert!(norm(sub(e.v, [-2.0 * (0.4f64).sin() / k, 2.0 * (0.4f64).cos() / k, 0.0])) < 1e-6);
        // ruled between two constant curves: A = 0, B = (1, 2, 3)
        let h = |vals: &[f64], dim: usize| Hermite { k: vec![0.0, 1.0], dim, p: [vals, vals].concat(), d1: vec![0.0; 2 * dim], d2: Some(vec![0.0; 2 * dim]) };
        let r = CgmRuled { h: h(&[0.0, 0.0, 0.0, 1.0, 2.0, 3.0], 6) };
        assert!(close(r.ev(0.3, 2.0).p, [2.0, 4.0, 6.0]));
        // fillet: quarter arc from p1 = (1,0,0) to p2 = (0,1,0) about c = 0; v0 = 0, v1 = 1
        let f = CgmFillet { h: h(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2], 10), v0: 0.0, v1: 1.0 };
        for v in [0.0, 0.25, 1.0] { let a = std::f64::consts::FRAC_PI_2 * v; assert!(close(f.ev(0.5, v).p, [a.cos(), a.sin(), 0.0])); }
        // 3D B-spline profile (a straight degree-1 segment) extruded along z
        let (c, n) = m.curve3(&[3.0, 1.0, 4.0, 0.0, 0.0, 1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0]).unwrap();
        assert_eq!(n, 15);
        let ex = CgmExtrusion { c, d: [0.0, 0.0, 2.0] };
        assert!(close(ex.ev(0.25, 0.5).p, [0.0, 0.5, 1.0]));
        m.surfaces.push(Arc::new(CgmPlane { o: [0.0; 3], u: [1.0, 0.0, 0.0], v: [0.0, 1.0, 0.0] }));
        m.pcurves.push((0, Arc::new(PCurve::Line { o: [1.0, 2.0], d: [1.0, 0.0] })));
        assert!(close(m.curve3(&[1.0, 0.0]).unwrap().0.at(2.0), [3.0, 2.0, 0.0]));
    }
}
