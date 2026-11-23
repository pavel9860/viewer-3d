//! Face triangulation in parameter space and exact-surface volume integration.
use crate::exact::orient2d;
use crate::geom::*;
use crate::topo::*;
use crate::trim::{area, Trim, P2};

fn o2(a: P2, b: P2, c: P2) -> i8 { orient2d(a, b, c).0 }
fn crosses(a: P2, b: P2, c: P2, d: P2) -> bool { o2(a, b, c) * o2(a, b, d) < 0 && o2(c, d, a) * o2(c, d, b) < 0 }

/// ear clipping of a counter-clockwise polygon with clockwise holes (each hole bridged to a visible vertex)
pub fn triangulate(outer: &[P2], holes: &[Vec<P2>]) -> (Vec<P2>, Vec<[usize; 3]>) {
    let mut pts: Vec<P2> = outer.to_vec();
    let mut poly: Vec<usize> = (0..pts.len()).collect();
    let mut hs: Vec<Vec<usize>> = holes.iter().map(|h| { let s = pts.len(); pts.extend_from_slice(h); (s..pts.len()).collect() }).collect();
    let maxx = |h: &Vec<usize>, pts: &[P2]| h.iter().map(|&i| pts[i][0]).fold(f64::MIN, f64::max);
    hs.sort_by(|a, b| maxx(b, &pts).total_cmp(&maxx(a, &pts)));
    for h in hs {
        let hi = (0..h.len()).max_by(|&a, &b| pts[h[a]][0].total_cmp(&pts[h[b]][0])).unwrap();
        let hp = pts[h[hi]];
        let clear = |q: P2, ring: &[usize]| (0..ring.len()).all(|k| { let (a, b) = (pts[ring[k]], pts[ring[(k + 1) % ring.len()]]); !crosses(hp, q, a, b) });
        let d2 = |j: usize| { let q = pts[poly[j]]; (q[0] - hp[0]).powi(2) + (q[1] - hp[1]).powi(2) };
        let mut order: Vec<usize> = (0..poly.len()).collect();
        order.sort_by(|&a, &b| d2(a).total_cmp(&d2(b)));
        let j = order.iter().copied().find(|&j| clear(pts[poly[j]], &poly) && clear(pts[poly[j]], &h)).unwrap_or(order[0]);
        let mut np = poly[..=j].to_vec();
        np.extend((0..=h.len()).map(|k| h[(hi + k) % h.len()]));
        np.extend_from_slice(&poly[j..]);
        poly = np;
    }
    let mut tris = vec![];
    let (mut v, mut i, mut misses) = (poly, 0usize, 0usize);
    while v.len() > 3 {
        let n = v.len();
        i %= n;
        let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
        let turn = o2(pts[a], pts[b], pts[c]);
        if turn == 0 && pts[a] != pts[c] { v.remove(i); misses = 0; i = i.saturating_sub(1); continue; }
        // only reflex vertices can lie inside a convex corner
        let ear = turn > 0 && (0..n).all(|k| {
            let q = v[k];
            if q == a || q == b || q == c || pts[q] == pts[a] || pts[q] == pts[b] || pts[q] == pts[c] { return true; }
            if o2(pts[v[(k + n - 1) % n]], pts[q], pts[v[(k + 1) % n]]) > 0 { return true; }
            !(o2(pts[a], pts[b], pts[q]) >= 0 && o2(pts[b], pts[c], pts[q]) >= 0 && o2(pts[c], pts[a], pts[q]) >= 0)
        });
        if ear { tris.push([a, b, c]); v.remove(i); misses = 0; i = i.saturating_sub(1); } else { i += 1; misses += 1; if misses > n { break; } }
    }
    if v.len() == 3 { tris.push([v[0], v[1], v[2]]); }
    (pts, tris)
}

/// (1/3) ∫∫ p · (Su x Sv) over a UV triangle (signed by its UV orientation): degree-5 Dunavant rule after
/// splitting until UV edges < 0.25
fn integrate(s: &dyn Surf, a: P2, b: P2, c: P2, depth: u32) -> f64 {
    let l = |p: P2, q: P2| (p[0] - q[0]).hypot(p[1] - q[1]);
    if depth < 16 && l(a, b).max(l(b, c)).max(l(c, a)) > 0.25 {
        let m = |p: P2, q: P2| [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0];
        let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
        return integrate(s, a, ab, ca, depth + 1) + integrate(s, ab, b, bc, depth + 1) + integrate(s, ca, bc, c, depth + 1) + integrate(s, ab, bc, ca, depth + 1);
    }
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
            let (pts, tris) = triangulate(o, &hs);
            let want = area(o) + hs.iter().map(|h| area(h)).sum::<f64>();
            let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]])).sum();
            if (got - want).abs() > 1e-9 * want.abs().max(1e-12) { bad += 1; }
            sum += tris.iter().map(|tr| integrate(s, t.raw(pts[tr[0]]), t.raw(pts[tr[1]]), t.raw(pts[tr[2]]), 0)).sum::<f64>();
        }
        sum
    }).sum();
    (v, bad)
}
