//! Refinement of a face's UV triangulation to a chordal tolerance, shared by the STEP (`tess`) and CGM (`cgm_mesh`)
//! meshers.
//!
//! Each pass: Lawson flips to the Delaunay triangulation in the surface's curvature metric, then every triangle
//! whose interior edge midpoints or centroid lie more than the chord off the surface has its longest interior edge
//! bisected, with Rivara closure (a triangle with a marked edge also splits its own longest interior edge), so
//! triangle shapes stay bounded. Chordal checks alone accept needle triangles whose vertices lie on the surface but
//! whose planes cut across it (Schwarz's lantern: area grows without bound while every vertex is exact); complete
//! Lawson flips in the metric and longest-edge bisection keep the facets tangent to the surface.
use crate::geom::*;
use crate::trim::P2;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug)]
pub struct Opts {
    /// chordal tolerance (model units)
    pub chord: f64,
    /// refinement passes
    pub max_passes: u32,
    /// triangle budget
    pub max_tris: usize,
}

/// result of `refine`
#[derive(Clone, Copy, Debug, Default)]
pub struct Outcome {
    /// the pass or triangle budget ran out before every triangle met the chord
    pub coarse: bool,
    /// largest measured chord deviation (edge midpoints and triangle centroids)
    pub dev: f64,
}

fn ekey(a: usize, c: usize) -> (usize, usize) { (a.min(c), a.max(c)) }

/// Lawson flips to the Delaunay triangulation in the surface's local metric (UV through `scale`'s map at each quad);
/// boundary edges stay. Improves shape without moving points. An edge stack with incremental adjacency: after a
/// flip only the quad's four outer edges are rechecked, so a fan onto a straight row (ear clipping's typical
/// output, whose flips propagate one triangle at a time) still reaches the Delaunay triangulation. The metric
/// varies from quad to quad (no global Delaunay to converge to in general), so the flips also have a budget.
pub fn flip(tri: &mut [[usize; 3]], uv: &[P2], scale: &dyn Fn(P2) -> [f64; 2], boundary: &HashSet<(usize, usize)>) {
    let mut adj: HashMap<(usize, usize), [usize; 2]> = HashMap::with_capacity(tri.len() * 2);
    let mut lone: HashSet<(usize, usize)> = HashSet::new();
    for (ti, t) in tri.iter().enumerate() {
        for k in 0..3 {
            let e = ekey(t[k], t[(k + 1) % 3]);
            match adj.get_mut(&e) { Some(v) if v[1] == usize::MAX => v[1] = ti, Some(_) => { lone.insert(e); } None => { adj.insert(e, [ti, usize::MAX]); } }
        }
    }
    // edges of more than two triangles (non-manifold input) never flip
    let mut stack: Vec<(usize, usize)> = adj.iter().filter(|(e, v)| v[1] != usize::MAX && !boundary.contains(e) && !lone.contains(e)).map(|(e, _)| *e).collect();
    stack.sort_unstable_by(|x, y| y.cmp(x));
    let mut budget = 64 * tri.len() + 1024;
    while let Some(e) = stack.pop() {
        let Some(&[t1, t2]) = adj.get(&e) else { continue };
        if t2 == usize::MAX || boundary.contains(&e) || lone.contains(&e) { continue; }
        // t1 = (a, c, b) with the edge a -> c, t2 = (c, a, d)
        let Some(k1) = (0..3).find(|&k| ekey(tri[t1][k], tri[t1][(k + 1) % 3]) == e) else { continue };
        let (a, c, b) = (tri[t1][k1], tri[t1][(k1 + 1) % 3], tri[t1][(k1 + 2) % 3]);
        let Some(k2) = (0..3).find(|&k| tri[t2][k] == c && tri[t2][(k + 1) % 3] == a) else { continue };
        let d = tri[t2][(k2 + 2) % 3];
        if b == d || adj.contains_key(&ekey(b, d)) { continue; }
        // the metric at the quad's centre: both diagonals are judged in one metric, so a flip is never undone by
        // the reverse flip (which a metric taken at the diagonal's own midpoint allows: flips then cycle)
        let m = [(uv[a][0] + uv[b][0] + uv[c][0] + uv[d][0]) / 4.0, (uv[a][1] + uv[b][1] + uv[c][1] + uv[d][1]) / 4.0];
        let sc = scale(m);
        let q = |i: usize| apply(sc, uv[i]);
        let (pa, pc, pb, pd) = (q(a), q(c), q(b), q(d));
        // the flipped pair must stay counter-clockwise (convex quad)
        if crate::exact::orient2d(pa, pd, pb).0 <= 0 || crate::exact::orient2d(pd, pc, pb).0 <= 0 { continue; }
        // d strictly inside the circumcircle of (a, c, b): a robust predicate (nearly cocircular points are common in
        // a strongly anisotropic metric)
        if crate::exact::incircle(pa, pc, pb, pd).0 <= 0 { continue; }
        if budget == 0 { break; }
        budget -= 1;
        tri[t1] = [a, d, b];
        tri[t2] = [d, c, b];
        adj.remove(&e);
        adj.insert(ekey(b, d), [t1, t2]);
        let mut swap = |x: (usize, usize), from: usize, to: usize| if let Some(v) = adj.get_mut(&x) { for s in v.iter_mut() { if *s == from { *s = to; } } };
        swap(ekey(a, d), t2, t1);
        swap(ekey(c, b), t1, t2);
        stack.extend([ekey(a, d), ekey(d, c), ekey(c, b), ekey(b, a)]);
    }
}

/// Anisotropic metric at `q`: each parameter direction scaled by its speed and the root of its normal curvature (the
/// chordal error of a step h is ~ k h^2 / 8, so equal metric lengths mean equal chord errors), floored so that
/// straight directions allow aspect ratios up to 100. Needs true second derivatives from `ev`: without them every
/// direction reads as straight and the metric collapses. `ev` evaluates the surface at a UV point of the
/// triangulation.
pub fn metric(ev: &dyn Fn(P2) -> D2, q: P2) -> [f64; 2] {
    let d = ev(q);
    let n = unit(cross(d.u, d.v));
    let (lu, lv) = (norm(d.u).max(1e-300), norm(d.v).max(1e-300));
    let (ku, kv) = ((dot(n, d.uu) / (lu * lu)).abs(), (dot(n, d.vv) / (lv * lv)).abs());
    let k0 = 1e-4 * ku.max(kv).max(1e-300);
    [lu * ku.max(k0).sqrt(), lv * kv.max(k0).sqrt()]
}

/// `uv` through the metric `m`
fn apply(m: [f64; 2], uv: P2) -> P2 { [m[0] * uv[0], m[1] * uv[1]] }

/// Refines the triangulation `tri` of the points `uv` (3D positions `pos`, on the surface `ev`) until every triangle
/// is within `o.chord` of the surface (measured between the surface and the facets through the surface points at
/// the vertices' UV). Boundary edges (used by one triangle of the given triangulation) are never
/// split or flipped, so the face keeps the edge points it shares with its neighbours. `same_edge(a, c)`: the two
/// points lie on one B-rep edge, so an interior edge between them would duplicate a boundary edge of the
/// neighbouring face (non-manifold connectivity); such edges are always split. New points are appended to `uv` and
/// `pos` in order; `tri` is replaced.
pub fn refine(ev: &dyn Fn(P2) -> D2, uv: &mut Vec<P2>, pos: &mut Vec<V>, tri: &mut Vec<[usize; 3]>, same_edge: &dyn Fn(usize, usize) -> bool, o: &Opts) -> Outcome {
    let at = |q: P2| ev(q).p;
    let mid = |a: V, c: V| mulv(addv(a, c), 0.5);
    let mut uses: HashMap<(usize, usize), u32> = HashMap::new();
    for t in tri.iter() { for k in 0..3 { *uses.entry(ekey(t[k], t[(k + 1) % 3])).or_default() += 1; } }
    let boundary: HashSet<(usize, usize)> = uses.into_iter().filter(|x| x.1 == 1).map(|x| x.0).collect();
    let mut out = Outcome::default();
    // the surface at each point's UV: errors are the surface's own interpolation error. Boundary points are the
    // shared edge points, which may lie off this surface by up to the tolerance the caller accepted; measuring
    // against them would ask refinement to remove an offset it cannot change (it never converges near the boundary)
    let mut on: Vec<V> = uv.iter().map(|&q| at(q)).collect();
    let scale = |q: P2| metric(ev, q);
    // per-edge chord error and metric length, kept across passes (edges are identified by their end points)
    let mut err: HashMap<(usize, usize), f64> = HashMap::new();
    let mut len: HashMap<(usize, usize), f64> = HashMap::new();
    for pass in 0..=o.max_passes {
        flip(tri, uv, &scale, &boundary);
        for t in tri.iter() {
            for k in 0..3 {
                let (a, c) = (t[k], t[(k + 1) % 3]);
                let e = ekey(a, c);
                if err.contains_key(&e) { continue; }
                let m = [(uv[a][0] + uv[c][0]) / 2.0, (uv[a][1] + uv[c][1]) / 2.0];
                err.insert(e, norm(sub(at(m), mid(on[a], on[c]))));
                // length in the curvature metric (a seam-crossing chord of zero 3D length is still long)
                let sc = scale(m);
                let d = apply(sc, [uv[a][0] - uv[c][0], uv[a][1] - uv[c][1]]);
                len.insert(e, d[0].hypot(d[1]));
            }
        }
        // longest interior edge of each triangle (index k of edge t[k]-t[k+1]); None when all edges are boundary
        let lk: Vec<Option<usize>> = tri.iter().map(|t| (0..3).filter(|&k| !boundary.contains(&ekey(t[k], t[(k + 1) % 3]))).max_by(|&i, &j| {
            let l = |k: usize| len[&ekey(t[k], t[(k + 1) % 3])];
            l(i).total_cmp(&l(j)).then(ekey(t[j], t[(j + 1) % 3]).cmp(&ekey(t[i], t[(i + 1) % 3])))
        })).collect();
        let mut marks: HashMap<(usize, usize), usize> = HashMap::new();
        let mut centroids: HashSet<usize> = Default::default();
        out.dev = 0.0;
        for (ti, t) in tri.iter().enumerate() {
            let mut bad = false;
            for k in 0..3 {
                let (a, c) = (t[k], t[(k + 1) % 3]);
                let d = err[&ekey(a, c)];
                out.dev = out.dev.max(d);
                let interior = !boundary.contains(&ekey(a, c));
                bad |= d > o.chord && interior;
                if interior && same_edge(a, c) { marks.insert(ekey(a, c), usize::MAX); }
            }
            let cu = [(uv[t[0]][0] + uv[t[1]][0] + uv[t[2]][0]) / 3.0, (uv[t[0]][1] + uv[t[1]][1] + uv[t[2]][1]) / 3.0];
            let dc = norm(sub(at(cu), mulv(addv(addv(on[t[0]], on[t[1]]), on[t[2]]), 1.0 / 3.0)));
            out.dev = out.dev.max(dc);
            bad |= dc > o.chord;
            if bad { match lk[ti] { Some(k) => { marks.insert(ekey(t[k], t[(k + 1) % 3]), usize::MAX); } None => { centroids.insert(ti); } } }
        }
        if marks.is_empty() && centroids.is_empty() { break; }
        if pass == o.max_passes || tri.len() > o.max_tris { out.coarse = true; break; }
        // closure: a triangle with a marked edge also splits its longest interior edge (Rivara), keeping shapes bounded
        loop {
            let mut grew = false;
            for (ti, t) in tri.iter().enumerate() {
                if (0..3).any(|k| marks.contains_key(&ekey(t[k], t[(k + 1) % 3]))) {
                    if let Some(k) = lk[ti] { if marks.insert(ekey(t[k], t[(k + 1) % 3]), usize::MAX).is_none() { grew = true; } }
                }
            }
            if !grew { break; }
        }
        let push = |q: P2, uv: &mut Vec<P2>, pos: &mut Vec<V>, on: &mut Vec<V>| { let p = at(q); uv.push(q); pos.push(p); on.push(p); uv.len() - 1 };
        let mut me: Vec<(usize, usize)> = marks.keys().copied().collect();
        me.sort_unstable();
        for e in me { let q = [(uv[e.0][0] + uv[e.1][0]) / 2.0, (uv[e.0][1] + uv[e.1][1]) / 2.0]; let m = push(q, uv, pos, &mut on); marks.insert(e, m); }
        let mut next = Vec::with_capacity(tri.len() * 2);
        for (ti, t) in tri.iter().enumerate() {
            if centroids.contains(&ti) {
                let q = [(uv[t[0]][0] + uv[t[1]][0] + uv[t[2]][0]) / 3.0, (uv[t[0]][1] + uv[t[1]][1] + uv[t[2]][1]) / 3.0];
                let m = push(q, uv, pos, &mut on);
                for k in 0..3 { next.push([t[k], t[(k + 1) % 3], m]); }
                continue;
            }
            let mk = |a: usize, c: usize| marks.get(&ekey(a, c)).copied();
            let Some(k) = lk[ti].filter(|&k| mk(t[k], t[(k + 1) % 3]).is_some()) else { next.push(*t); continue };
            // bisect the longest edge, then the other marked edges inside the two halves
            let (a, c, d) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            let ml = mk(a, c).unwrap();
            match mk(c, d) { Some(m2) => { next.push([ml, c, m2]); next.push([ml, m2, d]); } None => next.push([ml, c, d]) }
            match mk(d, a) { Some(m3) => { next.push([a, ml, m3]); next.push([ml, d, m3]); } None => next.push([a, ml, d]) }
        }
        *tri = next;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    /// a thin half-cylinder strip (radius 40, height 5, u = arc length) sampled along its rims, ear-clipped (a fan onto
    /// the rows) and refined: the facets must lie along the surface (area and normals), not across it (Schwarz's
    /// lantern: before the complete flips this strip meshed to 1.3x its area with facets 66x their UV area)
    #[test]
    fn thin_cylinder_strip_has_no_lantern() {
        let r = 40.0;
        let ev = |q: P2| {
            let (s, c) = (q[0] / r).sin_cos();
            D2 { p: [q[1], r * c, r * s], u: [0.0, -s, c], v: [1.0, 0.0, 0.0], uu: [0.0, -c / r, -s / r], uv: [0.0; 3], vv: [0.0; 3] }
        };
        let (u0, u1, n) = (std::f64::consts::PI * r, 2.0 * std::f64::consts::PI * r, 32);
        let mut ring: Vec<P2> = (0..=n).map(|k| [u0 + (u1 - u0) * k as f64 / n as f64, -2.5]).collect();
        ring.extend([[u1, -1.25], [u1, 0.0], [u1, 1.25]]);
        ring.extend((0..=n).rev().map(|k| [u0 + (u1 - u0) * k as f64 / n as f64, 2.5]));
        ring.extend([[u0, 1.25], [u0, 0.0], [u0, -1.25]]);
        let (mut uv, mut tri) = crate::tri::triangulate(&[ring]);
        let mut pos: Vec<V> = uv.iter().map(|&q| ev(q).p).collect();
        let out = refine(&ev, &mut uv, &mut pos, &mut tri, &|_, _| false, &Opts { chord: 0.11, max_passes: 16, max_tris: 200_000 });
        assert!(!out.coarse && out.dev <= 0.11, "{out:?}");
        let exact = (u1 - u0) * 5.0;
        let mut a = 0.0;
        for t in &tri {
            let [p, q, s] = t.map(|i| pos[i]);
            let n = cross(sub(q, p), sub(s, p));
            a += 0.5 * norm(n);
            let c = [(uv[t[0]][0] + uv[t[1]][0] + uv[t[2]][0]) / 3.0, (uv[t[0]][1] + uv[t[1]][1] + uv[t[2]][1]) / 3.0];
            let d = ev(c);
            let sn = unit(cross(d.u, d.v));
            assert!(dot(unit(n), sn).abs() > 0.99, "facet {t:?} tilted off the surface");
        }
        assert!((a - exact).abs() < 1e-3 * exact, "area {a} of {exact}");
        assert!(tri.len() < 200, "{} triangles for a strip", tri.len());
    }
}
