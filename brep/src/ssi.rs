//! Surface-surface intersection.
//! Transversal marching (4x4 Newton) + singular handling from second-order (curvature) analysis:
//! branch points are located and curves split there; tangential curves are traced with a least-squares
//! corrector; curves are clipped exactly at domain boundaries; parametric poles are crossed, not treated
//! as edges. Tolerant mode (model tolerance) treats contacts within tolerance as tangential.
use crate::geom::*;
use std::f64::consts::PI;

/// chordal tolerance of traced polylines (model units)
pub const SAG: f64 = 1e-5;

// ---------------- linear algebra ----------------
pub fn solve_n<const N: usize>(mut a: [[f64; N]; N], mut b: [f64; N]) -> Option<[f64; N]> {
    for c in 0..N {
        let p = (c..N).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs())).unwrap();
        if !a[p][c].is_finite() || a[p][c].abs() < 1e-300 { return None; }
        a.swap(c, p); b.swap(c, p);
        for r in c + 1..N { let f = a[r][c] / a[c][c]; for k in c..N { a[r][k] -= f * a[c][k]; } b[r] -= f * b[c]; }
    }
    let mut x = [0.0; N];
    for c in (0..N).rev() { let mut s = b[c]; for k in c + 1..N { s -= a[c][k] * x[k]; } x[c] = s / a[c][c]; }
    Some(x)
}
/// Damped Gauss-Newton on residual r(x) with a numeric Jacobian.
pub fn gauss_newton(x: [f64; 4], m: usize, r: &dyn Fn(&[f64; 4], &mut [f64]), iters: usize) -> [f64; 4] { gauss_newton_stall(x, m, r, iters, 0.0) }
/// Gauss-Newton that also stops when a step reduces the squared residual by less than the fraction `stall`
/// (a positive minimum: linear convergence on a nonzero residual)
pub fn gauss_newton_stall(mut x: [f64; 4], m: usize, r: &dyn Fn(&[f64; 4], &mut [f64]), iters: usize, stall: f64) -> [f64; 4] {
    let (mut f0, mut f1) = (vec![0.0; m], vec![0.0; m]);
    let mut lam = 1e-12;
    for _ in 0..iters {
        r(&x, &mut f0);
        let c0: f64 = f0.iter().map(|v| v * v).sum();
        if !(c0 >= 1e-32) { break; }
        let mut j = vec![[0.0; 4]; m];
        for k in 0..4 {
            let h = 1e-7 * (1.0 + x[k].abs());
            let mut xp = x; xp[k] += h;
            r(&xp, &mut f1);
            for i in 0..m { j[i][k] = (f1[i] - f0[i]) / h; }
        }
        let mut a = [[0.0; 4]; 4];
        let mut g = [0.0; 4];
        for i in 0..m { for p in 0..4 { g[p] -= j[i][p] * f0[i]; for q in 0..4 { a[p][q] += j[i][p] * j[i][q]; } } }
        let tr = (0..4).map(|k| a[k][k]).sum::<f64>().max(1e-300);
        let mut improved = false;
        for _ in 0..12 {
            let mut ad = a;
            for k in 0..4 { ad[k][k] += lam * tr; }
            let Some(dx) = solve_n(ad, g) else { lam *= 10.0; continue };
            let xn = [x[0] + dx[0], x[1] + dx[1], x[2] + dx[2], x[3] + dx[3]];
            r(&xn, &mut f1);
            let c1: f64 = f1.iter().map(|v| v * v).sum();
            if c1 < c0 { x = xn; lam = (lam * 0.1).max(1e-15); improved = c1 < (1.0 - stall) * c0; break; }
            lam *= 10.0;
        }
        if !improved { break; }
    }
    x
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum End { Closed, Boundary, Branch, Tangent, Fail }
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Kind { Transversal, Tangential, Isolated }
enum Sing { Isolated, Tangential(V), Branch(V, V), Degenerate }
pub struct Curve { pub pts: Vec<V>, pub ends: (End, End), pub kind: Kind }

/// tol: distance tolerance (model units); ang: sine of the normal angle below which surfaces count as tangent
/// sag: chordal tolerance of the traced polylines
struct Ctx<'a> { a: &'a dyn Surf, b: &'a dyn Surf, tol: f64, ang: f64, sag: f64, reg: std::cell::RefCell<Vec<[f64; 4]>> }

fn ls2(su: V, sv: V, d: V) -> (f64, f64) {
    let (p, q, r) = (dot(su, su), dot(su, sv), dot(sv, sv));
    let (r1, r2) = (dot(su, d), dot(sv, d));
    let det = p * r - q * q;
    ((r * r1 - q * r2) / det, (p * r2 - q * r1) / det)
}

impl Ctx<'_> {
    fn rec(&self, x: [f64; 4]) { self.reg.borrow_mut().push(x); }
    fn eval(&self, x: &[f64; 4]) -> (D2, D2) { (self.a.ev(x[0], x[1]), self.b.ev(x[2], x[3])) }
    fn gap(&self, x: &[f64; 4]) -> f64 { let (a, b) = self.eval(x); norm(sub(a.p, b.p)) }
    /// sine of the angle between the normals
    fn sin_n(&self, x: &[f64; 4]) -> f64 { let (a, b) = self.eval(x); let (na, nb) = (cross(a.u, a.v), cross(b.u, b.v)); norm(cross(na, nb)) / (norm(na) * norm(nb)) }
    fn wrap(&self, x: &mut [f64; 4]) {
        for (s, o) in [(self.a, 0), (self.b, 2)] {
            if s.poles_v() { let d = s.dom(); if x[o + 1] > d[3] { x[o + 1] = 2.0 * d[3] - x[o + 1]; x[o] += PI; } else if x[o + 1] < d[2] { x[o + 1] = 2.0 * d[2] - x[o + 1]; x[o] += PI; } }
            if s.per_u() { x[o] = x[o].rem_euclid(2.0 * PI); }
            if s.per_v() { x[o + 1] = x[o + 1].rem_euclid(2.0 * PI); }
        }
    }
    /// first violated bound: (param index, bound value)
    fn violated(&self, x: &[f64; 4]) -> Option<(usize, f64)> {
        for (s, o) in [(self.a, 0), (self.b, 2)] {
            let d = s.dom();
            if !s.per_u() { if x[o] < d[0] { return Some((o, d[0])); } if x[o] > d[1] { return Some((o, d[1])); } }
            if !s.per_v() && !s.poles_v() { if x[o + 1] < d[2] { return Some((o + 1, d[2])); } if x[o + 1] > d[3] { return Some((o + 1, d[3])); } }
        }
        None
    }
    /// exact coverage test: is p on an already traced curve? Project p onto the curve through the nearest
    /// traced points (in the plane through p orthogonal to the local tangent) and compare.
    fn on_traced(&self, p: V) -> bool {
        let reg = self.reg.borrow();
        let mut near: Vec<(f64, [f64; 4])> = reg.iter().map(|x| (norm(sub(self.a.ev(x[0], x[1]).p, p)), *x)).filter(|d| d.0 < 0.2).collect();
        near.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (d, x) in near.into_iter().take(4) {
            if d < 1e-9 { return true; }
            let (a, b) = self.eval(&x);
            let tv = cross(cross(a.u, a.v), cross(b.u, b.v));
            let t = if norm(tv) > 1e-12 { unit(tv) } else { unit(sub(p, a.p)) };
            // p is on the traced part only if it is within one tracer step of x and the curve between them stays in
            // the domain (a curve leaving a bounded window and re-entering it is a different traced part)
            let same_arc = |c: &Self| d <= 0.1 * 1.01 && {
                let m = mulv(addv(a.p, p), 0.5);
                c.correct(x, m, t).map_or(true, |z| c.violated(&z).is_none())
            };
            if let Some(y) = self.correct(x, p, t) { if norm(sub(self.a.ev(y[0], y[1]).p, p)) < 1e-7 && same_arc(self) { return true; } }
            let r = |y: &[f64; 4], f: &mut [f64]| { let (a, b) = self.eval(y); f[..3].copy_from_slice(&sub(a.p, b.p)); f[3] = dot(t, sub(a.p, p)); };
            let y = gauss_newton(x, 4, &r, 40);
            if self.gap(&y) < self.tol && norm(sub(self.a.ev(y[0], y[1]).p, p)) < 1e-5 && same_arc(self) { return true; }
            if self.ang > 1e-6 {
                // tolerant contact: the seed lies somewhere in the contact band; project it onto the ridge with the
                // tracer's own least-squares functional and compare with the traced ridge
                let rr = |y: &[f64; 4], f: &mut [f64]| { let (a, b) = self.eval(y); f[..3].copy_from_slice(&sub(a.p, b.p)); f[3..6].copy_from_slice(&cross(unit(cross(a.u, a.v)), unit(cross(b.u, b.v)))); f[6] = dot(t, sub(a.p, p)); };
                let y = gauss_newton(x, 7, &rr, 60);
                let q = self.a.ev(y[0], y[1]).p;
                if self.gap(&y) < self.tol && norm(sub(q, p)).min(d) < 0.2 && reg.iter().any(|z| norm(sub(self.a.ev(z[0], z[1]).p, q)) < 0.05) && dot(sub(q, p), t).abs() < 1e-6 { return true; }
            }
        }
        false
    }
    /// project onto A ∩ B within the plane through q orthogonal to t
    fn correct(&self, x: [f64; 4], q: V, t: V) -> Option<[f64; 4]> {
        let r = |x: &[f64; 4], f: &mut [f64]| { let (a, b) = self.eval(x); f[..3].copy_from_slice(&sub(a.p, b.p)); f[3] = dot(t, sub(a.p, q)); };
        let mut y = x;
        for _ in 0..12 {
            let (a, b) = self.eval(&y);
            let f = sub(a.p, b.p);
            let g = dot(t, sub(a.p, q));
            if norm(f) < 1e-14 && g.abs() < 1e-14 { return Some(y); }
            let m = [[a.u[0], a.v[0], -b.u[0], -b.v[0]], [a.u[1], a.v[1], -b.u[1], -b.v[1]], [a.u[2], a.v[2], -b.u[2], -b.v[2]], [dot(t, a.u), dot(t, a.v), 0.0, 0.0]];
            let dx = solve_n(m, [-f[0], -f[1], -f[2], -g])?;
            for k in 0..4 { y[k] += dx[k]; }
            if !(dx.iter().map(|v| v.abs()).fold(0.0, f64::max) <= 1.0) { return None; }
        }
        // ill-conditioned (near a singular point): least squares, accepted only on both surfaces
        let y = gauss_newton(x, 4, &r, 40);
        if self.gap(&y) < self.tol { Some(y) } else { None }
    }
    /// nearby point where the surfaces meet with parallel normals
    fn singular_near(&self, x: [f64; 4]) -> Option<[f64; 4]> {
        let r = |x: &[f64; 4], f: &mut [f64]| {
            let (a, b) = self.eval(x);
            f[..3].copy_from_slice(&sub(a.p, b.p));
            f[3..].copy_from_slice(&cross(unit(cross(a.u, a.v)), unit(cross(b.u, b.v))));
        };
        let y = gauss_newton(x, 6, &r, 60);
        if self.gap(&y) < self.tol && self.sin_n(&y) < self.ang { Some(y) } else { None }
    }
    /// difference of the second fundamental forms in the common tangent frame (e1, e2)
    fn qform(&self, x: &[f64; 4]) -> (V, V, f64, f64, f64) {
        let (a, b) = self.eval(x);
        let n = unit(cross(a.u, a.v));
        let e1 = unit(a.u);
        let e2 = cross(n, e1);
        let form = |s: &D2| -> [f64; 3] {
            let mut ns = unit(cross(s.u, s.v));
            if dot(ns, n) < 0.0 { ns = mulv(ns, -1.0); }
            let (l, m, nn) = (dot(s.uu, ns), dot(s.uv, ns), dot(s.vv, ns));
            let (g11, g12, g22) = (dot(s.u, s.u), dot(s.u, s.v), dot(s.v, s.v));
            let det = g11 * g22 - g12 * g12;
            let to_uv = |w: V| { let (r1, r2) = (dot(s.u, w), dot(s.v, w)); ((g22 * r1 - g12 * r2) / det, (g11 * r2 - g12 * r1) / det) };
            let (p, q) = (to_uv(e1), to_uv(e2));
            let ii = |a: (f64, f64), b: (f64, f64)| l * a.0 * b.0 + m * (a.0 * b.1 + a.1 * b.0) + nn * a.1 * b.1;
            [ii(p, p), ii(p, q), ii(q, q)]
        };
        let (qa, qb) = (form(&a), form(&b));
        (e1, e2, qa[0] - qb[0], qa[1] - qb[1], qa[2] - qb[2])
    }
    /// direction along a tangential curve: eigenvector of the smallest-magnitude eigenvalue of the form
    fn tangential_dir(&self, x: &[f64; 4]) -> V {
        let (e1, e2, a, b, c) = self.qform(x);
        let m = (a + c) / 2.0;
        let r = (((a - c) / 2.0).powi(2) + b * b).sqrt();
        let l = if (m - r).abs() < (m + r).abs() { m - r } else { m + r };
        let (vx, vy) = if (a - l).abs() + b.abs() > (c - l).abs() + b.abs() { (-b, a - l) } else { (c - l, -b) };
        let (vx, vy) = if vx == 0.0 && vy == 0.0 { (1.0, 0.0) } else { (vx, vy) };
        unit(addv(mulv(e1, vx), mulv(e2, vy)))
    }
    /// second-order analysis at a point with parallel normals: directions where normal curvatures agree
    fn classify(&self, x: &[f64; 4]) -> Sing {
        let (e1, e2, q11, q12, q22) = self.qform(x);
        let sc = q11.abs().max(q12.abs()).max(q22.abs());
        if !(sc >= 1e-9) { return Sing::Degenerate; }
        let (q11, q12, q22) = (q11 / sc, q12 / sc, q22 / sc);
        let disc = q12 * q12 - q11 * q22;
        let dir = |c: f64, s: f64| unit(addv(mulv(e1, c), mulv(e2, s)));
        if disc < -1e-8 { return Sing::Isolated; }
        if disc <= 1e-8 { return if q11.abs() >= q22.abs() { Sing::Tangential(dir(-q12, q11)) } else { Sing::Tangential(dir(q22, -q12)) }; }
        let r = disc.sqrt();
        if q22.abs() > 1e-12 { Sing::Branch(dir(q22, -q12 + r), dir(q22, -q12 - r)) } else { Sing::Branch(dir(0.0, 1.0), dir(-2.0 * q12, q11)) }
    }
    /// exact point where the curve leaves the domain: fix the violated parameter on its bound, solve A = B
    fn clip(&self, x1: [f64; 4]) -> Option<[f64; 4]> {
        let (k, val) = self.violated(&x1)?;
        let r = |x: &[f64; 4], f: &mut [f64]| { let mut y = *x; y[k] = val; let (a, b) = self.eval(&y); f[..3].copy_from_slice(&sub(a.p, b.p)); f[3] = (x[k] - val) * 1e3; };
        let mut s = x1; s[k] = val;
        let mut y = gauss_newton(s, 4, &r, 40);
        y[k] = val;
        if self.gap(&y) < self.tol && self.violated(&y).is_none() { Some(y) } else { None }
    }
    /// march along a transversal curve; stops at closure, domain boundary or a branch point
    fn march(&self, seed: [f64; 4], t0: V, start: V, branches: &mut Vec<[f64; 4]>) -> (Vec<V>, End) {
        let mut x = seed;
        let mut h: f64 = 1e-3;
        let mut pts = vec![self.a.ev(x[0], x[1]).p]; self.rec(x);
        let mut tp = t0;
        let mut s_prev = self.sin_n(&x);
        for step in 0..200_000 {
            let (a, b) = self.eval(&x);
            let tv = cross(cross(a.u, a.v), cross(b.u, b.v));
            let mut t = if norm(tv) > 1e-300 { unit(tv) } else { tp };
            if dot(t, tp) < 0.0 { t = mulv(t, -1.0); }
            loop {
                if h < 1e-9 { return (pts, End::Fail); }
                let d = mulv(t, h);
                let (dua, dva) = ls2(a.u, a.v, d);
                let (dub, dvb) = ls2(b.u, b.v, d);
                let Some(mut nx) = self.correct([x[0] + dua, x[1] + dva, x[2] + dub, x[3] + dvb], addv(a.p, d), t) else { h *= 0.5; continue };
                let p = self.a.ev(nx[0], nx[1]).p;
                let chord = sub(p, pts[pts.len() - 1]);
                if norm(chord) < 0.5 * h || dot(unit(chord), t) < 0.99 { h *= 0.5; continue; } // jumped branch / turned too fast
                // singular point between x and nx?
                let s_new = self.sin_n(&nx);
                if s_new.min(s_prev) < 0.05 {
                    let from = if s_new < s_prev { nx } else { x };
                    if let Some(bp) = self.singular_near(from) {
                        let q = self.a.ev(bp[0], bp[1]).p;
                        let last = pts[pts.len() - 1];
                        let seg = sub(p, last);
                        let tt = (dot(sub(q, last), seg) / dot(seg, seg)).clamp(0.0, 1.0);
                        if (norm(sub(q, addv(last, mulv(seg, tt)))) < 1e-6 || (norm(sub(q, p)) < 1.5 * h && dot(sub(q, p), t) > 0.0)) && norm(sub(q, start)) > 1e-9 {
                            pts.push(q); self.rec(bp);
                            if !branches.iter().any(|b| norm(sub(self.a.ev(b[0], b[1]).p, q)) < 1e-7) { branches.push(bp); }
                            return (pts, End::Branch);
                        }
                    }
                }
                if self.violated(&nx).is_some() {
                    if let Some(c) = self.clip(nx) { pts.push(self.a.ev(c[0], c[1]).p); self.rec(c); return (pts, End::Boundary); }
                    h *= 0.5; continue;
                }
                self.wrap(&mut nx);
                if step > 3 && norm(sub(p, start)) < h * 1.01 && dot(sub(start, pts[pts.len() - 1]), t) > 0.0 { pts.push(start); return (pts, End::Closed); }
                x = nx;
                pts.push(p); self.rec(nx);
                let used = norm(chord);
                let (na2, nb2) = (self.a.ev(nx[0], nx[1]), self.b.ev(nx[2], nx[3]));
                let tn = cross(cross(na2.u, na2.v), cross(nb2.u, nb2.v));
                let turn = if norm(tn) > 1e-300 { dot(unit(tn), t).abs().min(1.0).acos() } else { 0.0 };
                let kappa = turn / used.max(1e-300);
                tp = t;
                h = (h * 1.5).min(0.1);
                if kappa > 0.0 { h = h.min((8.0 * self.sag / kappa).sqrt()); } // chordal deviation <= sag
                // predicted distance to a zero of sin(normal angle): never step over a singular point
                if s_new < s_prev { let reach = s_new * used / (s_prev - s_new); h = h.min((0.5 * reach).max(1e-5)); }
                s_prev = s_new;
                break;
            }
        }
        (pts, End::Fail)
    }
    /// march along a tangential curve (normals parallel along it, exactly or within tolerance)
    fn march_tangent(&self, seed: [f64; 4], d0: V, start: V) -> (Vec<V>, End) {
        let mut x = seed;
        let mut pts = vec![start]; self.rec(seed);
        let mut t = d0;
        let mut h: f64 = 1e-3;
        for step in 0..20_000 {
            let (a, b) = self.eval(&x);
            let p = a.p;
            let nx = loop {
                if h < 1e-9 { return (pts, End::Fail); }
                let q = addv(p, mulv(t, h));
                let r = |y: &[f64; 4], f: &mut [f64]| {
                    let (a, b) = self.eval(y);
                    f[..3].copy_from_slice(&sub(a.p, b.p));
                    f[3..6].copy_from_slice(&cross(unit(cross(a.u, a.v)), unit(cross(b.u, b.v))));
                    f[6] = dot(t, sub(a.p, q));
                };
                let d = mulv(t, h);
                let (dua, dva) = ls2(a.u, a.v, d);
                let (dub, dvb) = ls2(b.u, b.v, d);
                let nx = gauss_newton([x[0] + dua, x[1] + dva, x[2] + dub, x[3] + dvb], 7, &r, 60);
                let np = self.a.ev(nx[0], nx[1]).p;
                let chord = sub(np, p);
                if !(self.gap(&nx) <= self.tol) || norm(chord) < 0.5 * h || dot(unit(chord), t) < 0.99 { h *= 0.5; continue; }
                break nx;
            };
            let mut nx = nx;
            if let Some((k, val)) = self.violated(&nx) {
                // fraction of the step inside the domain; a step running along the bound clips at its end
                let den = nx[k] - x[k];
                let f = if den != 0.0 { ((val - x[k]) / den).clamp(0.0, 1.0) } else { 1.0 };
                let mut c = [0.0; 4];
                for i in 0..4 { c[i] = x[i] + f * (nx[i] - x[i]); }
                c[k] = val;
                let r = |y: &[f64; 4], ff: &mut [f64]| { let mut z = *y; z[k] = val; let (a, b) = self.eval(&z); ff[..3].copy_from_slice(&sub(a.p, b.p)); ff[3..6].copy_from_slice(&cross(unit(cross(a.u, a.v)), unit(cross(b.u, b.v)))); };
                let mut c2 = gauss_newton(c, 6, &r, 40); c2[k] = val;
                if !(self.gap(&c2) <= self.tol) { c2 = c; }
                pts.push(self.a.ev(c2[0], c2[1]).p); self.rec(c2);
                return (pts, End::Boundary);
            }
            self.wrap(&mut nx);
            let np = self.a.ev(nx[0], nx[1]).p;
            if step > 3 && norm(sub(np, start)) < h * 1.01 && dot(sub(start, p), t) > 0.0 { pts.push(start); return (pts, End::Closed); }
            pts.push(np); self.rec(nx);
            let nd = self.tangential_dir(&nx);
            let nd = if dot(nd, t) < 0.0 { mulv(nd, -1.0) } else { nd };
            let used = norm(sub(np, p));
            // curvature only from steps long enough for the turn to be meaningful (tolerant ridges are noisy)
            let kappa = if used > 1e-5 { dot(nd, t).min(1.0).acos() / used } else { 0.0 };
            x = nx;
            t = nd;
            h = (h * 1.5).min(0.1);
            if kappa > 0.0 { h = h.min((8.0 * self.sag / kappa).sqrt()); }
            h = h.max(1e-5);
            if self.sin_n(&x) > 2.0 * self.ang.max(5e-5) { return (pts, End::Tangent); } // left the contact
        }
        (pts, End::Fail)
    }
}

pub fn intersect(a: &dyn Surf, b: &dyn Surf) -> (Vec<Curve>, usize) { intersect_tol(a, b, 1e-10, 1e-6) }

/// Intersection curves of two bounded surfaces. Tolerant mode (`ang` > 1e-6): contacts within `tol`
/// (distance) and `ang` (sine of the normal angle) count as tangential.
pub fn intersect_tol(a: &dyn Surf, b: &dyn Surf, tol: f64, ang: f64) -> (Vec<Curve>, usize) {
    intersect_sag(a, b, tol, ang, SAG)
}
/// `intersect_tol` with traced polylines within `sag` of the curves (model-dependent accuracy)
pub fn intersect_sag(a: &dyn Surf, b: &dyn Surf, tol: f64, ang: f64, sag: f64) -> (Vec<Curve>, usize) {
    let c = Ctx { a, b, tol, ang, sag, reg: Default::default() };
    let n = 24;
    let samp = |s: &dyn Surf| -> Vec<(f64, f64, V)> {
        let d = s.dom();
        let mut v = vec![];
        for i in 0..n { for j in 0..n {
            let (u, w) = (d[0] + (d[1] - d[0]) * i as f64 / (n - 1) as f64, d[2] + (d[3] - d[2]) * j as f64 / (n - 1) as f64);
            v.push((u, w, s.ev(u, w).p));
        } }
        v
    };
    let (sa, sb) = (samp(a), samp(b));
    // largest grid cell (corner span, doubled for curvature bulge): every surface point is within it of a sample
    let cell = |g: &[(f64, f64, V)]| {
        let mut m: f64 = 0.0;
        for i in 0..n - 1 { for j in 0..n - 1 {
            let q = [g[i * n + j].2, g[(i + 1) * n + j].2, g[i * n + j + 1].2, g[(i + 1) * n + j + 1].2];
            for x in 0..4 { for y in x + 1..4 { m = m.max(norm(sub(q[x], q[y]))); } }
        } }
        2.0 * m
    };
    let reach = cell(&sa) + cell(&sb) + tol;
    let mut curves: Vec<Curve> = vec![];
    let mut branches: Vec<[f64; 4]> = vec![];
    let covered = |curves: &Vec<Curve>, p: V, tol: f64| curves.iter().any(|cv| {
        if cv.pts.len() == 1 { return norm(sub(cv.pts[0], p)) < tol; }
        cv.pts.windows(2).any(|w| { let ab = sub(w[1], w[0]); let t = (dot(sub(p, w[0]), ab) / dot(ab, ab).max(1e-300)).clamp(0.0, 1.0); norm(sub(p, addv(w[0], mulv(ab, t)))) < tol + 4.0 * c.sag })
    });
    let r3 = |x: &[f64; 4], f: &mut [f64]| { let (pa, pb) = (a.ev(x[0], x[1]).p, b.ev(x[2], x[3]).p); f[..3].copy_from_slice(&sub(pa, pb)); };
    let trace_from = |s: [f64; 4], curves: &mut Vec<Curve>, branches: &mut Vec<[f64; 4]>| {
        let p = a.ev(s[0], s[1]).p;
        if covered(curves, p, 1e-4) || c.on_traced(p) { return; }
        if c.sin_n(&s) < c.ang {
            let sing = if c.ang > 1e-6 {
                match c.classify(&s) { Sing::Isolated if c.sin_n(&s) > 1e-9 => Sing::Tangential(c.tangential_dir(&s)), Sing::Branch(..) => Sing::Tangential(c.tangential_dir(&s)), o => o }
            } else { c.classify(&s) };
            match sing {
                Sing::Tangential(d) => {
                    let (f, e1) = c.march_tangent(s, d, p);
                    if e1 == End::Closed { curves.push(Curve { pts: f, ends: (e1, e1), kind: Kind::Tangential }); return; }
                    let (bk, e2) = c.march_tangent(s, mulv(d, -1.0), p);
                    let mut pts: Vec<V> = bk.into_iter().rev().collect();
                    pts.pop(); pts.extend(f);
                    curves.push(Curve { pts, ends: (e2, e1), kind: Kind::Tangential });
                }
                Sing::Isolated => curves.push(Curve { pts: vec![p], ends: (End::Tangent, End::Tangent), kind: Kind::Isolated }),
                Sing::Branch(..) => { if !branches.iter().any(|b| norm(sub(a.ev(b[0], b[1]).p, p)) < 1e-7) { branches.push(s); } }
                Sing::Degenerate => curves.push(Curve { pts: vec![p], ends: (End::Fail, End::Fail), kind: Kind::Isolated }), // coincident
            }
            return;
        }
        let (aa, bb) = c.eval(&s);
        let t = unit(cross(cross(aa.u, aa.v), cross(bb.u, bb.v)));
        let (f, e1) = c.march(s, t, p, branches);
        if e1 == End::Closed { curves.push(Curve { pts: f, ends: (e1, e1), kind: Kind::Transversal }); return; }
        let (bk, e2) = c.march(s, mulv(t, -1.0), p, branches);
        let mut pts: Vec<V> = bk.into_iter().rev().collect();
        pts.pop(); pts.extend(f);
        curves.push(Curve { pts, ends: (e2, e1), kind: Kind::Transversal });
    };
    // seeds from both grids (each sample with the nearest sample of the other surface), so a narrow surface
    // window does not see only the crossings its few samples converge to
    let near = |g: &[(f64, f64, V)], p: V| *g.iter().min_by(|x, y| norm(sub(x.2, p)).total_cmp(&norm(sub(y.2, p)))).unwrap();
    let seeds: Vec<[f64; 4]> = sa.iter().map(|&(ua, va, pa)| { let (ub, vb, pb) = near(&sb, pa); (norm(sub(pb, pa)), [ua, va, ub, vb]) })
        .chain(sb.iter().map(|&(ub, vb, pb)| { let (ua, va, pa) = near(&sa, pb); (norm(sub(pb, pa)), [ua, va, ub, vb]) }))
        // no intersection within one cell of either sample: nothing to find from this seed
        .filter(|x| x.0 <= reach).map(|x| x.1).collect();
    for [ua, va, ub, vb] in seeds {
        let mut s = gauss_newton_stall([ua, va, ub, vb], 3, &r3, 80, 0.01);
        // a seed that converged outside the domain: the curve may still enter the domain, through the bound it
        // violates; pull it onto that bound (intersection with the parameter fixed there)
        for _ in 0..2 {
            let Some((k, val)) = c.violated(&s) else { break };
            let rb = |x: &[f64; 4], f: &mut [f64]| { r3(x, &mut f[..3]); f[3] = x[k] - val; };
            let mut y = s;
            y[k] = val;
            s = gauss_newton(y, 4, &rb, 40);
            s[k] = val;
        }
        // NaN-safe: a NaN gap must never pass a tolerance test
        if !(c.gap(&s) <= c.tol) || s.iter().any(|v| !v.is_finite()) || c.violated(&s).is_some() { continue; }
        c.wrap(&mut s);
        // a seed next to a singular point: snap to it so branch points are found even when seeds miss them
        if c.sin_n(&s) < 1e-3f64.max(c.ang) { if let Some(bp) = c.singular_near(s) { s = bp; } }
        trace_from(s, &mut curves, &mut branches);
    }
    // arcs leaving every branch point along both branch directions
    let mut done = 0;
    while done < branches.len() {
        let bp = branches[done];
        done += 1;
        let Sing::Branch(d1, d2) = c.classify(&bp) else { continue };
        let p0 = a.ev(bp[0], bp[1]).p;
        for d in [d1, mulv(d1, -1.0), d2, mulv(d2, -1.0)] {
            let h0 = 1e-3;
            let q = addv(p0, mulv(d, h0));
            if covered(&curves, q, 1e-5) || c.on_traced(q) { continue; }
            let (aa, bb) = c.eval(&bp);
            let (dua, dva) = ls2(aa.u, aa.v, mulv(d, h0));
            let (dub, dvb) = ls2(bb.u, bb.v, mulv(d, h0));
            let Some(mut x) = c.correct([bp[0] + dua, bp[1] + dva, bp[2] + dub, bp[3] + dvb], q, d) else { continue };
            c.wrap(&mut x);
            let mut br2 = branches.clone();
            let (mut f, e) = c.march(x, d, p0, &mut br2);
            for nb in br2 { if !branches.iter().any(|b| norm(sub(a.ev(b[0], b[1]).p, a.ev(nb[0], nb[1]).p)) < 1e-7) { branches.push(nb); } }
            f.insert(0, p0);
            curves.push(Curve { pts: f, ends: (End::Branch, e), kind: Kind::Transversal });
        }
    }
    // curves seeded before their branch point was known pass through it: duplicates of the arcs
    let nb = branches.len();
    let bpts: Vec<V> = branches.iter().map(|b| a.ev(b[0], b[1]).p).collect();
    let keep = curves.into_iter().filter(|cv| {
        let through = cv.pts.len() > 1 && bpts.iter().any(|q| cv.pts.windows(2).any(|w| { let ab = sub(w[1], w[0]); let t = (dot(sub(*q, w[0]), ab) / dot(ab, ab)).clamp(0.0, 1.0); norm(sub(*q, addv(w[0], mulv(ab, t)))) < 1e-6 }));
        !(through && cv.ends.0 != End::Branch && cv.ends.1 == End::Branch) && !(through && cv.ends.0 == End::Branch && cv.ends.1 != End::Branch && false)
    }).collect();
    (keep, nb)
}
