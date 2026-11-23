//! CATIA CGM exact B-rep, read from the flat IR written by cadcore (`cadcli cgm PART OUT.cgmir`,
//! format in cadcore/src/formats/catia/cgm_model.hpp). Surfaces and p-curves keep CGM's own
//! parametrizations; faces are bounded by p-curves, so each coedge is a parameter range on one.
use crate::geom::*;
use crate::nurbs::BspS;
use std::sync::Arc;

/// p-curve in a surface's (u, v)
pub enum PCurve {
    Line { o: [f64; 2], d: [f64; 2] },
    /// angle = t / (r s)
    Circle { c: [f64; 2], r: f64, s: f64 },
    /// C2 piecewise quintic Hermite (cubic when `d2` is empty): values and derivatives at the knots
    Hermite { k: Vec<f64>, p: [Vec<f64>; 2], d1: [Vec<f64>; 2], d2: Option<[Vec<f64>; 2]> },
}
impl PCurve {
    pub fn uv(&self, t: f64) -> [f64; 2] {
        match self {
            PCurve::Line { o, d } => [o[0] + t * d[0], o[1] + t * d[1]],
            PCurve::Circle { c, r, s } => { let a = t / (r * s); [c[0] + r * a.cos(), c[1] + r * a.sin()] }
            PCurve::Hermite { k, p, d1, d2 } => {
                let n = k.len();
                let mut i = 0;
                while i + 2 < n && t > k[i + 1] { i += 1; }
                let h = k[i + 1] - k[i];
                let s = if h != 0.0 { (t - k[i]) / h } else { 0.0 };
                let (s2, s3) = (s * s, s * s * s);
                let (s4, s5) = (s3 * s, s3 * s2);
                [0, 1].map(|c| match d2 {
                    Some(d2) => {
                        let h0 = 1.0 - 10.0 * s3 + 15.0 * s4 - 6.0 * s5;
                        let h1 = s - 6.0 * s3 + 8.0 * s4 - 3.0 * s5;
                        let h2 = 0.5 * s2 - 1.5 * s3 + 1.5 * s4 - 0.5 * s5;
                        let h5 = 10.0 * s3 - 15.0 * s4 + 6.0 * s5;
                        let h4 = -4.0 * s3 + 7.0 * s4 - 3.0 * s5;
                        let h3 = 0.5 * s3 - s4 + 0.5 * s5;
                        h0 * p[c][i] + h * h1 * d1[c][i] + h * h * h2 * d2[c][i] + h5 * p[c][i + 1] + h * h4 * d1[c][i + 1] + h * h * h3 * d2[c][i + 1]
                    }
                    None => {
                        let (a0, a1, b0, b1) = (2.0 * s3 - 3.0 * s2 + 1.0, s3 - 2.0 * s2 + s, -2.0 * s3 + 3.0 * s2, s3 - s2);
                        a0 * p[c][i] + h * a1 * d1[c][i] + b0 * p[c][i + 1] + h * b1 * d1[c][i + 1]
                    }
                })
            }
        }
    }
}

/// p = o + u U + v V
pub struct CgmPlane { pub o: V, pub u: V, pub v: V }
/// p = o + r (cos a X + sin a Y) + v Z, a = u / su, Z = X x Y
pub struct CgmCyl { pub o: V, pub x: V, pub y: V, pub z: V, pub r: f64, pub su: f64 }
/// p = o + cos b (cos a A + sin a B) + sin b C, a = u / k, b = v / k (A, B, C scaled by the radius)
pub struct CgmSphere { pub o: V, pub a: V, pub b: V, pub c: V, pub k: f64 }
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
    fn uv_of(&self, p: V) -> [f64; 2] { uv_newton(self, p) }
}

pub struct Coedge { pub pcurve: usize, pub edge: usize, pub sense: i32, pub t: [f64; 2] }
pub struct Face { pub surface: usize, pub sense: i32, pub loops: Vec<std::ops::Range<usize>>, pub id: u32 }
/// `body`: the owning CGM body record (0 if none); `solid`: bounds a volume; otherwise a sheet (open surfaces, reference planes)
pub struct Shell { pub id: u32, pub body: u32, pub solid: bool, pub faces: Vec<(usize, i32)> }
pub struct Skip { pub id: u32, pub kind: u32, pub why: String }

/// the B-rep of one part
pub struct Model {
    pub surfaces: Vec<Arc<dyn Surf>>,
    /// (surface, curve)
    pub pcurves: Vec<(usize, PCurve)>,
    pub vertices: usize,
    /// vertex indices (start, end)
    pub edges: Vec<[usize; 2]>,
    pub coedges: Vec<Coedge>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub skipped: Vec<Skip>,
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

fn v3(p: &[f64], i: usize) -> V { [p[i], p[i + 1], p[i + 2]] }

fn surface(kind: u32, p: &[f64]) -> Result<Arc<dyn Surf>, String> {
    let need = |n: usize| if p.len() < n { Err(format!("surface kind {kind}: {} parameters", p.len())) } else { Ok(()) };
    Ok(match kind {
        1 => { need(9)?; Arc::new(CgmPlane { o: v3(p, 0), u: v3(p, 3), v: v3(p, 6) }) }
        2 => { need(11)?; let (x, y) = (v3(p, 3), v3(p, 6)); Arc::new(CgmCyl { o: v3(p, 0), x, y, z: cross(x, y), r: p[9], su: p[10] }) }
        3 => { need(13)?; Arc::new(CgmSphere { o: v3(p, 0), a: v3(p, 3), b: v3(p, 6), c: v3(p, 9), k: p[12] }) }
        4 => { need(14)?; Arc::new(CgmTorus { o: v3(p, 0), x: v3(p, 3), y: v3(p, 6), z: v3(p, 9), big: p[12], r: p[13] }) }
        5 => {
            need(5)?;
            let (du, dv, nu, nv, rat) = (p[0] as usize, p[1] as usize, p[2] as usize, p[3] as usize, p[4] != 0.0);
            let (ku, kv) = (nu + du + 1, nv + dv + 1);
            need(5 + ku + kv + 3 * nu * nv + if rat { nu * nv } else { 0 })?;
            let k0 = 5;
            let poles = k0 + ku + kv;
            let wts = poles + 3 * nu * nv;
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
        3 => {
            let n = *p.first().ok_or_else(bad)? as usize;
            if p.len() < 1 + 5 * n + 1 { return Err(bad()); }
            let k = p[1..1 + n].to_vec();
            let arr = |i: usize| p[1 + n + i * n..1 + n + (i + 1) * n].to_vec();
            let has2 = p[1 + 5 * n] != 0.0;
            let d2 = if has2 {
                if p.len() < 2 + 7 * n { return Err(bad()); }
                let b = 2 + 5 * n;
                Some([p[b..b + n].to_vec(), p[b + n..b + 2 * n].to_vec()])
            } else { None };
            PCurve::Hermite { k, p: [arr(0), arr(1)], d1: [arr(2), arr(3)], d2 }
        }
        _ => return Err(bad()),
    })
}

impl Model {
    pub fn read(b: &[u8]) -> Result<Model, String> {
        if b.len() < 8 || &b[..8] != b"CGMIR\0\0\x01" { return Err("not a CGM IR (version 1)".into()); }
        let mut r = Rd { b, k: 8 };
        let c: Vec<usize> = (0..9).map(|_| r.n()).collect::<Result<_, _>>()?;
        let mut m = Model { surfaces: vec![], pcurves: vec![], vertices: c[2], edges: vec![], coedges: vec![], faces: vec![], shells: vec![], skipped: vec![] };
        for _ in 0..c[0] { let (kind, _id, n) = (r.u32()?, r.u32()?, r.n()?); let p = r.f64s(n)?; m.surfaces.push(surface(kind, &p)?); }
        for _ in 0..c[1] { let (kind, s, _id, n) = (r.u32()?, r.n()?, r.u32()?, r.n()?); let p = r.f64s(n)?; if s >= m.surfaces.len() { return Err("p-curve surface out of range".into()); } m.pcurves.push((s, pcurve(kind, &p)?)); }
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
            let (id, body, solid, n) = (r.u32()?, r.u32()?, r.u32()? != 0, r.n()?);
            let faces = (0..n).map(|_| { let (f, s) = (r.n()?, r.i32()?); if f >= m.faces.len() { Err("shell face out of range".to_string()) } else { Ok((f, s)) } }).collect::<Result<_, _>>()?;
            m.shells.push(Shell { id, body, solid, faces });
        }
        for _ in 0..c[8] { let (id, kind, n) = (r.u32()?, r.u32()?, r.n()?); let why = String::from_utf8_lossy(r.take(n)?).into_owned(); m.skipped.push(Skip { id, kind, why }); }
        if r.k != b.len() { return Err(format!("{} bytes after the IR", b.len() - r.k)); }
        Ok(m)
    }
    /// 3D point of coedge `c` at parameter t
    pub fn point(&self, c: &Coedge, t: f64) -> V {
        let (s, pc) = &self.pcurves[c.pcurve];
        let uv = pc.uv(t);
        self.surfaces[*s].ev(uv[0], uv[1]).p
    }
}
