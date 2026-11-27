//! Watertight tessellation of B-rep faces.
//!
//! Every edge is tessellated once (the same parameters for each face that uses it; the end points are the
//! topological vertices), so neighbouring faces share their boundary points exactly. Each face is triangulated in
//! its normalized UV frame (sweep-line triangulation of the trimmed loops) and refined by conforming red-green
//! subdivision of interior edges until every triangle is within the chordal tolerance of the surface. The face
//! meshes are merged on their shared points into one indexed mesh.
use crate::geom::*;
use crate::par;
use crate::topo::*;
use crate::tri::triangulate_src;
use crate::trim::{Trim, P2};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug)]
pub struct Opts {
    /// chordal tolerance (model units): edges and triangles stay this close to the exact geometry
    pub chord: f64,
    /// refinement passes per face
    pub max_passes: u32,
    /// triangle budget per face
    pub max_tris: usize,
}
impl Default for Opts { fn default() -> Self { Opts { chord: 1e-3, max_passes: 16, max_tris: 200_000 } } }

/// identity of a mesh point: a topological vertex, the k-th tessellation point of an edge, or a face-interior point
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Key { V(u32), E(u32, u32), F(u32, u32) }

/// indexed triangle mesh; `face[i]` is the B-rep face of triangle i (index into the face list given)
#[derive(Default)]
pub struct Mesh { pub pos: Vec<V>, pub tri: Vec<[u32; 3]>, pub face: Vec<u32> }

#[derive(Default, Debug, Clone)]
pub struct Stats {
    pub faces: usize,
    /// faces whose loops gave no triangles
    pub empty: usize,
    /// faces that hit the pass or triangle budget before meeting the chordal tolerance
    pub coarse: usize,
    /// largest measured chord deviation (at edge midpoints and triangle centroids)
    pub max_dev: f64,
}

/// Lawson flips to the Delaunay triangulation in the surface's local metric (UV scaled by |Su|, |Sv| at each quad);
/// boundary edges stay. Improves shape without moving points.
fn flip(tri: &mut [[usize; 3]], uv: &[P2], scale: &dyn Fn(P2) -> [f64; 2], boundary: &std::collections::HashSet<(usize, usize)>) {
    let ekey = |a: usize, c: usize| (a.min(c), a.max(c));
    for _ in 0..32 {
        let mut adj: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
        for (ti, t) in tri.iter().enumerate() { for k in 0..3 { adj.entry(ekey(t[k], t[(k + 1) % 3])).or_default().push((ti, k)); } }
        let mut done = vec![false; tri.len()];
        let mut flips = 0;
        let mut edges: Vec<_> = adj.into_iter().filter(|(e, v)| v.len() == 2 && !boundary.contains(e)).collect();
        edges.sort_unstable_by_key(|x| x.0);
        for (_, v) in edges {
            let ((t1, k1), (t2, k2)) = (v[0], v[1]);
            if done[t1] || done[t2] { continue; }
            let (a, c, b) = (tri[t1][k1], tri[t1][(k1 + 1) % 3], tri[t1][(k1 + 2) % 3]);
            let d = tri[t2][(k2 + 2) % 3];
            if tri[t2][k2] != c || tri[t2][(k2 + 1) % 3] != a { continue; }
            let m = [(uv[a][0] + uv[c][0]) / 2.0, (uv[a][1] + uv[c][1]) / 2.0];
            let sc = scale(m);
            let q = |i: usize| [uv[i][0] * sc[0], uv[i][1] * sc[1]];
            let (pa, pc, pb, pd) = (q(a), q(c), q(b), q(d));
            // the flipped pair must stay counter-clockwise (convex quad)
            if crate::exact::orient2d(pa, pd, pb).0 <= 0 || crate::exact::orient2d(pd, pc, pb).0 <= 0 { continue; }
            // d strictly inside the circumcircle of (a, c, b)
            let r = |p: [f64; 2]| [p[0] - pd[0], p[1] - pd[1], (p[0] - pd[0]).powi(2) + (p[1] - pd[1]).powi(2)];
            let (x, y, z) = (r(pa), r(pc), r(pb));
            let det = x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0]) + x[2] * (y[0] * z[1] - y[1] * z[0]);
            let scale2 = (x[2] + y[2] + z[2]).powi(2);
            if det > 1e-12 * scale2 {
                tri[t1] = [a, d, b];
                tri[t2] = [d, c, b];
                done[t1] = true;
                done[t2] = true;
                flips += 1;
            }
        }
        if flips == 0 { break; }
    }
}

struct FaceMesh { keys: Vec<Key>, pos: Vec<V>, tri: Vec<[usize; 3]>, coarse: bool, dev: f64 }


fn face_mesh(b: &Body, fi: u32, f: Id<Face>, o: &Opts, sag: &HashMap<Id<Edge>, f64>) -> FaceMesh {
    let fc = b.faces.get(f);
    let s = fc.surface.as_ref();
    let tr = Trim::with_edge_sag(b, s, fc.sense, &b.face_loops(f), &|e| sag.get(&e).copied().unwrap_or(o.chord));
    // ring points and their keys (the ring of a loop is its items' points without each item's last point)
    let mut rings: Vec<Vec<P2>> = vec![];
    let mut rkeys: Vec<Vec<(Key, V)>> = vec![];
    // B-rep edges each ring point lies on (a vertex lies on two)
    let mut ron: Vec<Vec<[u32; 2]>> = vec![];
    for l in &tr.loops {
        let (mut r, mut k) = (vec![], vec![]);
        let mut on = vec![];
        for (ii, it) in l.iter().enumerate() {
            let n = it.uv.len();
            // as Trim::ring: the last point stays when the next edge starts elsewhere in UV (a pole line)
            let upto = if it.uv[n - 1] != l[(ii + 1) % l.len()].uv[0] { n } else { n - 1 };
            for j in 0..upto {
                let key = if j == 0 { Key::V(it.v[0].0) } else if j == n - 1 { Key::V(it.v[1].0) } else { Key::E(it.e.0, if it.sense { j } else { n - 1 - j } as u32) };
                r.push(it.uv[j]);
                k.push((key, it.p[j]));
                let other = if j == 0 { l[(ii + l.len() - 1) % l.len()].e.0 } else if j == n - 1 { l[(ii + 1) % l.len()].e.0 } else { it.e.0 };
                on.push([it.e.0, other]);
            }
        }
        rings.push(r);
        rkeys.push(k);
        ron.push(on);
    }
    // a single clockwise loop (not wrapping a period) still bounds the region inside it: triangulate it reversed and
    // orient the triangles along the loop, so shared edges stay opposite to the neighbouring faces
    let cw = rings.len() == 1 && crate::trim::area(&rings[0]) < 0.0;
    if cw { rings[0].reverse(); rkeys[0].reverse(); }
    let (mut uv, src, mut tri) = triangulate_src(&rings);
    let mut keys: Vec<Key> = src.iter().map(|&(r, i)| rkeys[r][i].0).collect();
    let on: Vec<[u32; 2]> = src.iter().map(|&(r, i)| ron[r][i]).collect();
    // an interior edge between two points of the same B-rep edge would also exist in the neighbouring face
    // (non-manifold connectivity): such edges are always split
    let same_edge = |a: usize, c: usize| a < on.len() && c < on.len() && on[a].iter().any(|e| on[c].contains(e));
    let mut pos: Vec<V> = src.iter().map(|&(r, i)| rkeys[r][i].1).collect();
    let at = |q: P2| { let r = tr.raw(q); s.ev(r[0], r[1]).p };
    let mid = |a: V, c: V| mulv(addv(a, c), 0.5);
    let ekey = |a: usize, c: usize| (a.min(c), a.max(c));
    // boundary edges: used by one triangle of the initial triangulation
    let mut uses: HashMap<(usize, usize), u32> = HashMap::new();
    for t in &tri { for k in 0..3 { *uses.entry(ekey(t[k], t[(k + 1) % 3])).or_default() += 1; } }
    let boundary: std::collections::HashSet<(usize, usize)> = uses.into_iter().filter(|x| x.1 == 1).map(|x| x.0).collect();
    let (mut coarse, mut dev) = (false, 0.0f64);
    // anisotropic metric: each parameter direction scaled by its speed and the root of its normal curvature (the
    // chordal error of a step h is ~ k h^2 / 8), floored so that straight directions allow aspect ratios up to 100
    let scale = |q: P2| {
        let r = tr.raw(q);
        let d = s.ev(r[0], r[1]);
        let n = unit(cross(d.u, d.v));
        let (lu, lv) = (norm(d.u).max(1e-300), norm(d.v).max(1e-300));
        let (ku, kv) = ((dot(n, d.uu) / (lu * lu)).abs(), (dot(n, d.vv) / (lv * lv)).abs());
        let k0 = 1e-4 * ku.max(kv).max(1e-300);
        [lu * ku.max(k0).sqrt(), lv * kv.max(k0).sqrt()]
    };
    // per-edge chord error and metric length, kept across passes (edges are identified by their end points)
    let mut err: HashMap<(usize, usize), f64> = HashMap::new();
    let mut len: HashMap<(usize, usize), f64> = HashMap::new();
    for pass in 0..=o.max_passes {
        flip(&mut tri, &uv, &scale, &boundary);
        for t in &tri {
            for k in 0..3 {
                let (a, c) = (t[k], t[(k + 1) % 3]);
                let e = ekey(a, c);
                if err.contains_key(&e) { continue; }
                let m = [(uv[a][0] + uv[c][0]) / 2.0, (uv[a][1] + uv[c][1]) / 2.0];
                err.insert(e, norm(sub(at(m), mid(pos[a], pos[c]))));
                // length in the curvature metric (a seam-crossing chord of zero 3D length is still long)
                let sc = scale(m);
                len.insert(e, ((uv[a][0] - uv[c][0]) * sc[0]).hypot((uv[a][1] - uv[c][1]) * sc[1]));
            }
        }
        let e_of = |a: usize, c: usize, err: &mut HashMap<(usize, usize), f64>| err[&ekey(a, c)];
        // longest interior edge of each triangle (index k of edge t[k]-t[k+1]); None when all edges are boundary
        let lk: Vec<Option<usize>> = tri.iter().map(|t| (0..3).filter(|&k| !boundary.contains(&ekey(t[k], t[(k + 1) % 3]))).max_by(|&i, &j| {
            let l = |k: usize| len[&ekey(t[k], t[(k + 1) % 3])];
            l(i).total_cmp(&l(j)).then(ekey(t[j], t[(j + 1) % 3]).cmp(&ekey(t[i], t[(i + 1) % 3])))
        })).collect();
        let mut marks: HashMap<(usize, usize), usize> = HashMap::new();
        let mut centroids: std::collections::HashSet<usize> = Default::default();
        dev = 0.0;
        for (ti, t) in tri.iter().enumerate() {
            let mut bad = false;
            for k in 0..3 {
                let (a, c) = (t[k], t[(k + 1) % 3]);
                let d = e_of(a, c, &mut err);
                dev = dev.max(d);
                let interior = !boundary.contains(&ekey(a, c));
                bad |= d > o.chord && interior;
                if interior && same_edge(a, c) { marks.insert(ekey(a, c), usize::MAX); }
            }
            let cu = [(uv[t[0]][0] + uv[t[1]][0] + uv[t[2]][0]) / 3.0, (uv[t[0]][1] + uv[t[1]][1] + uv[t[2]][1]) / 3.0];
            let dc = norm(sub(at(cu), mulv(addv(addv(pos[t[0]], pos[t[1]]), pos[t[2]]), 1.0 / 3.0)));
            dev = dev.max(dc);
            bad |= dc > o.chord;
            if bad { match lk[ti] { Some(k) => { marks.insert(ekey(t[k], t[(k + 1) % 3]), usize::MAX); } None => { centroids.insert(ti); } } }
        }
        if marks.is_empty() && centroids.is_empty() { break; }
        if pass == o.max_passes || tri.len() > o.max_tris { coarse = true; break; }
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
        let push = |q: P2, uv: &mut Vec<P2>, pos: &mut Vec<V>, keys: &mut Vec<Key>| { uv.push(q); pos.push(at(q)); keys.push(Key::F(fi, keys.len() as u32)); uv.len() - 1 };
        let mut me: Vec<(usize, usize)> = marks.keys().copied().collect();
        me.sort_unstable();
        for e in me { let q = [(uv[e.0][0] + uv[e.1][0]) / 2.0, (uv[e.0][1] + uv[e.1][1]) / 2.0]; let m = push(q, &mut uv, &mut pos, &mut keys); marks.insert(e, m); }
        let mut out = Vec::with_capacity(tri.len() * 2);
        for (ti, t) in tri.iter().enumerate() {
            if centroids.contains(&ti) {
                let q = [(uv[t[0]][0] + uv[t[1]][0] + uv[t[2]][0]) / 3.0, (uv[t[0]][1] + uv[t[1]][1] + uv[t[2]][1]) / 3.0];
                let m = push(q, &mut uv, &mut pos, &mut keys);
                for k in 0..3 { out.push([t[k], t[(k + 1) % 3], m]); }
                continue;
            }
            let mk = |a: usize, c: usize| marks.get(&ekey(a, c)).copied();
            let Some(k) = lk[ti].filter(|&k| mk(t[k], t[(k + 1) % 3]).is_some()) else { out.push(*t); continue };
            // bisect the longest edge, then the other marked edges inside the two halves
            let (a, c, d) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            let ml = mk(a, c).unwrap();
            match mk(c, d) { Some(m2) => { out.push([ml, c, m2]); out.push([ml, m2, d]); } None => out.push([ml, c, d]) }
            match mk(d, a) { Some(m3) => { out.push([a, ml, m3]); out.push([ml, d, m3]); } None => out.push([a, ml, d]) }
        }
        tri = out;
    }
    if cw { for t in tri.iter_mut() { t.swap(1, 2); } }
    FaceMesh { keys, pos, tri, coarse, dev }
}

/// tessellate the given faces into one mesh, merged on shared vertices and edge points
pub fn tessellate(b: &Body, faces: &[Id<Face>], o: &Opts, workers: usize) -> (Mesh, Stats) {
    // edge tolerances: start at the chord; edges whose tessellation makes loops cross (loops closer than the chord)
    // or whose UV segments leave the chord tolerance are refined, for every face that uses them
    let mut sag: HashMap<Id<Edge>, f64> = HashMap::new();
    for _ in 0..8 {
        let bad: Vec<Vec<Id<Edge>>> = par::map(faces.len(), workers, |i| {
            let fc = b.faces.get(faces[i]);
            let t = Trim::with_edge_sag(b, fc.surface.as_ref(), fc.sense, &b.face_loops(faces[i]), &|e| sag.get(&e).copied().unwrap_or(o.chord));
            let mut v = t.crossing_edges();
            v.extend(t.bent_edges(fc.surface.as_ref(), o.chord));
            v
        });
        if crate::debug("tess") { eprintln!("refine round: {} edges in {} faces", bad.iter().map(|v| v.len()).sum::<usize>(), bad.iter().filter(|v| !v.is_empty()).count()); }
        let mut any = false;
        for e in bad.into_iter().flatten() { let v = sag.entry(e).or_insert(o.chord); *v *= 0.25; any = true; }
        if !any { break; }
    }
    let fms: Vec<FaceMesh> = par::map(faces.len(), workers, |i| face_mesh(b, i as u32, faces[i], o, &sag));
    let mut st = Stats { faces: faces.len(), ..Default::default() };
    let mut m = Mesh::default();
    let mut index: HashMap<Key, u32> = HashMap::new();
    for (fi, fm) in fms.into_iter().enumerate() {
        if fm.tri.is_empty() { st.empty += 1; }
        st.coarse += fm.coarse as usize;
        if fm.coarse && crate::debug("tess") { let f = b.faces.get(faces[fi]); eprintln!("coarse face STEP #{} dom {:?} per {} {} tris {} dev {:.2e}", f.origin.src, f.surface.dom(), f.surface.per_u(), f.surface.per_v(), fm.tri.len(), fm.dev); }
        st.max_dev = st.max_dev.max(fm.dev);
        let ids: Vec<u32> = fm.keys.iter().zip(&fm.pos).map(|(k, p)| *index.entry(*k).or_insert_with(|| { m.pos.push(*p); (m.pos.len() - 1) as u32 })).collect();
        for t in fm.tri {
            let t = [ids[t[0]], ids[t[1]], ids[t[2]]];
            if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] { continue; }
            m.tri.push(t);
            m.face.push(fi as u32);
        }
    }
    (m, st)
}

/// watertightness: directed edges without their reverse (open), and undirected edges used more than twice
pub fn open_edges(m: &Mesh) -> (usize, usize) {
    let mut d: HashMap<(u32, u32), i32> = HashMap::new();
    let mut u: HashMap<(u32, u32), u32> = HashMap::new();
    for t in &m.tri {
        for k in 0..3 {
            let (a, c) = (t[k], t[(k + 1) % 3]);
            *d.entry((a, c)).or_default() += 1;
            *d.entry((c, a)).or_default() -= 1;
            *u.entry((a.min(c), a.max(c))).or_default() += 1;
        }
    }
    (d.values().filter(|v| **v > 0).map(|v| *v as usize).sum(), u.values().filter(|v| **v > 2).count())
}

/// enclosed volume of a closed mesh (divergence theorem)
pub fn mesh_volume(m: &Mesh) -> f64 { m.tri.iter().map(|t| dot(m.pos[t[0] as usize], cross(m.pos[t[1] as usize], m.pos[t[2] as usize]))).sum::<f64>() / 6.0 }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::make::{block, cylinder};
    #[test]
    fn primitives_tessellate_watertight_and_converge() {
        for (b, v) in [(block([0.0; 3], [1.0, 2.0, 3.0]), 48.0), (cylinder([0.0; 3], 0.5, 2.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]), std::f64::consts::PI)] {
            let faces: Vec<Id<Face>> = b.faces.ids().collect();
            let (m, st) = tessellate(&b, &faces, &Opts { chord: 1e-4, ..Default::default() }, 2);
            assert_eq!(open_edges(&m), (0, 0));
            assert_eq!((st.empty, st.coarse), (0, 0));
            assert!(st.max_dev <= 1e-4, "{st:?}");
            assert!((mesh_volume(&m) - v).abs() / v < 1e-3, "volume {}", mesh_volume(&m));
        }
    }
}
