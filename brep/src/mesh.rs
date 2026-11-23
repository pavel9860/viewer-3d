//! Exact-surface volume integration over faces triangulated in parameter space.
use crate::geom::*;
use crate::topo::*;
use crate::trim::{area, Trim, P2};


/// degree-5 Dunavant rule for (1/3) ∫∫ p · (Su x Sv) over a UV triangle (signed by its UV orientation)
fn rule(s: &dyn Surf, a: P2, b: P2, c: P2) -> f64 {
    const R: [(f64, f64, f64); 7] = [
        (1.0 / 3.0, 1.0 / 3.0, 0.225),
        (0.059715871789770, 0.470142064105115, 0.132394152788506), (0.470142064105115, 0.059715871789770, 0.132394152788506), (0.470142064105115, 0.470142064105115, 0.132394152788506),
        (0.797426985353087, 0.101286507323456, 0.125939180544827), (0.101286507323456, 0.797426985353087, 0.125939180544827), (0.101286507323456, 0.101286507323456, 0.125939180544827),
    ];
    let jac = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.0;
    R.iter().map(|&(l1, l2, w)| {
        let l3 = 1.0 - l1 - l2;
        let d = s.ev(l1 * a[0] + l2 * b[0] + l3 * c[0], l1 * a[1] + l2 * b[1] + l3 * c[1]);
        w * dot(d.p, cross(d.u, d.v)) / 3.0
    }).sum::<f64>() * jac
}
/// adaptive quadrature: a triangle is split into 4 while the children's sum differs from its own rule by more than
/// `tol` (absolute, shared out over the children); independent of the parameterization's scale
fn integrate(s: &dyn Surf, a: P2, b: P2, c: P2, whole: f64, tol: f64, depth: u32) -> f64 {
    let m = |p: P2, q: P2| [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0];
    let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
    let k = [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)];
    let parts = k.map(|t| rule(s, t.0, t.1, t.2));
    let sum: f64 = parts.iter().sum();
    if depth >= 10 || (sum - whole).abs() <= tol { return sum; }
    (0..4).map(|i| integrate(s, k[i].0, k[i].1, k[i].2, parts[i], tol / 4.0, depth + 1)).sum()
}

pub fn face_trim(b: &Body, f: Id<Face>) -> Trim {
    let fc = b.faces.get(f);
    let loops: Vec<Vec<(Id<Edge>, bool)>> = fc.loops.iter().map(|&l| b.loops.get(l).coedges.iter().map(|&c| { let ce = b.coedges.get(c); (ce.edge, ce.sense) }).collect()).collect();
    Trim::new(b, fc.surface.as_ref(), fc.sense, &loops)
}

/// volume enclosed by a face set (divergence theorem); also reports faces whose triangulation lost area
pub fn volume(b: &Body, faces: &[Id<Face>]) -> (f64, usize) {
    let mut bad = 0;
    let v = faces.iter().map(|&f| {
        let t = face_trim(b, f);
        let s = b.faces.get(f).surface.as_ref();
        let polys = t.polygons();
        let (outer, holes): (Vec<&Vec<P2>>, Vec<&Vec<P2>>) = polys.iter().partition(|p| area(p) > 0.0);
        let mut sum = 0.0;
        for o in outer {
            let hs: Vec<Vec<P2>> = holes.iter().filter(|h| crate::trim::in_poly(h[0], o)).map(|h| h.to_vec()).collect();
            let rings: Vec<Vec<P2>> = std::iter::once(o.to_vec()).chain(hs.iter().cloned()).collect();
            let (pts, tris) = crate::tri::triangulate(&rings);
            let want = area(o) + hs.iter().map(|h| area(h)).sum::<f64>();
            let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]])).sum();
            if (got - want).abs() > 1e-9 * want.abs().max(1e-12) {
                bad += 1;
                if let Ok(d) = std::env::var("BREP_DUMP") { let _ = std::fs::write(format!("{d}/face{}.txt", b.faces.get(f).origin.src), std::iter::once(o).chain(hs.iter()).map(|r| r.iter().map(|q| format!("{} {}", q[0], q[1])).collect::<Vec<_>>().join(";")).collect::<Vec<_>>().join("\n")); }
                if std::env::var("BREP_DEBUG").is_ok() { eprintln!("triangulation lost area: face src {} outer {} pts, {} holes, area {want:.3e} got {got:.3e}", b.faces.get(f).origin.src, o.len(), hs.len()); }
            }
            sum += tris.iter().map(|tr| { let (a, b, c) = (t.raw(pts[tr[0]]), t.raw(pts[tr[1]]), t.raw(pts[tr[2]])); integrate(s, a, b, c, rule(s, a, b, c), 1e-12, 0) }).sum::<f64>();
        }
        sum
    }).sum();
    (v, bad)
}
