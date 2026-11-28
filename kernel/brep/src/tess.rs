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
    // the surface over the triangulation's UV (u mirrored when the face is reversed: derivatives in u follow)
    let ev = |q: P2| { let r = tr.raw(q); let d = s.ev(r[0], r[1]); D2 { u: mulv(d.u, tr.sg), uv: mulv(d.uv, tr.sg), ..d } };
    let n0 = uv.len();
    let ro = crate::refine::Opts { chord: o.chord, max_passes: o.max_passes, max_tris: o.max_tris };
    let r = crate::refine::refine(&ev, &mut uv, &mut pos, &mut tri, &same_edge, &ro);
    keys.extend((n0..uv.len()).map(|i| Key::F(fi, i as u32)));
    let (coarse, dev) = (r.coarse, r.dev);
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
