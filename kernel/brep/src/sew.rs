//! Sewing: open shells (e.g. faces exported one by one) become shells that share their common edges.
//!
//! Only free edges (one coedge) and their end vertices take part; closed shells stay as they are.
//! 1. End vertices of free edges within the tolerance are merged.
//! 2. A free edge is split at every end vertex of another free edge lying on it (T-junctions).
//! 3. Free edges between the same two vertices that run within the tolerance of each other are paired.
//!    Unambiguous pairs go first; where more edges coincide (duplicated faces) an edge pairs with one whose face
//!    is already joined to its own.
//! 4. Faces joined through shared edges form one shell. Each face's loops are turned to run counter-clockwise
//!    about its normal, then faces are flipped so that each shared edge is used once in each direction (keeping
//!    the majority orientation; a closed shell is turned outward). Shared edges that cannot be made consistent
//!    are counted as conflicts.
use crate::geom::*;
use crate::nurbs::t_of;
use crate::topo::*;
use std::collections::HashMap;

#[derive(Debug, Default, Clone)]
pub struct SewStats {
    pub vertices_merged: usize,
    pub edges_merged: usize,
    pub edges_split: usize,
    pub free_before: usize,
    pub free_after: usize,
    pub shells_before: usize,
    pub shells_after: usize,
    pub closed_after: usize,
    pub flipped: usize,
    pub loops_reversed: usize,
    /// shared edges whose two faces could not be oriented consistently
    pub conflicts: usize,
}

const NONE: u32 = u32::MAX;

/// whether edge `a` runs within `tol` of edge `e` (same end vertices), checked at ¼, ½, ¾; Some(same direction)
fn coincide(a: &Edge, e: &Edge, tol: f64) -> Option<bool> {
    let (ca, ce) = (a.curve.as_ref()?, e.curve.as_ref()?);
    let (lo, hi) = (e.t[0].min(e.t[1]), e.t[0].max(e.t[1]));
    let eps = 1e-9 * (hi - lo).max(1.0);
    for k in 1..4 {
        let p = ca.ev(a.t[0] + (a.t[1] - a.t[0]) * k as f64 / 4.0).0;
        let t = t_of(ce.as_ref(), p);
        if t < lo - eps || t > hi + eps || norm(sub(ce.ev(t).0, p)) > tol { return None; }
    }
    if a.v[0] != a.v[1] { return Some(a.v[0] == e.v[0]); }
    // closed edges: compare tangents at the middle
    let ta = ca.ev((a.t[0] + a.t[1]) / 2.0);
    let te = ce.ev(t_of(ce.as_ref(), ta.0));
    Some(dot(ta.1, te.1) * (a.t[1] - a.t[0]) * (e.t[1] - e.t[0]) > 0.0)
}

/// split every free edge at the end vertices of other free edges lying on it within `tol`; returns edges split
fn split_at_vertices(b: &mut Body, tol: f64) -> usize {
    let u = b.edge_uses();
    let free: Vec<usize> = (0..b.edges.len()).filter(|&e| u[e] == 1 && b.edges.get(Id::new(e)).curve.is_some()).collect();
    if free.is_empty() { return 0; }
    // boxes of 17 samples, grown by tol and a quarter of the longest step (curvature between samples)
    let boxes: Vec<[f64; 6]> = free.iter().map(|&e| {
        let ed = b.edges.get(Id::new(e));
        let c = ed.curve.as_ref().unwrap();
        let pts: Vec<V> = (0..=16).map(|k| c.ev(ed.t[0] + (ed.t[1] - ed.t[0]) * k as f64 / 16.0).0).collect();
        let g = tol + pts.windows(2).map(|w| norm(sub(w[1], w[0]))).fold(0.0, f64::max) / 4.0;
        let mut bx = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
        for p in pts { for k in 0..3 { bx[k] = bx[k].min(p[k] - g); bx[k + 3] = bx[k + 3].max(p[k] + g); } }
        bx
    }).collect();
    // grid of the median box size; boxes spanning too many cells go to one list checked for every vertex
    let mut sizes: Vec<f64> = boxes.iter().map(|x| (0..3).map(|k| x[k + 3] - x[k]).fold(0.0, f64::max)).collect();
    sizes.sort_unstable_by(f64::total_cmp);
    let h = sizes[sizes.len() / 2].max(tol * 4.0);
    let cell = |x: f64| (x / h).floor() as i64;
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    let mut huge = vec![];
    for (i, bx) in boxes.iter().enumerate() {
        let (a, z) = ([0, 1, 2].map(|k| cell(bx[k])), [0, 1, 2].map(|k| cell(bx[k + 3])));
        if (z[0] - a[0] + 1) * (z[1] - a[1] + 1) * (z[2] - a[2] + 1) > 4096 { huge.push(i as u32); continue; }
        for x in a[0]..=z[0] { for y in a[1]..=z[1] { for w in a[2]..=z[2] { grid.entry([x, y, w]).or_default().push(i as u32); } } }
    }
    let mut verts: Vec<u32> = free.iter().flat_map(|&e| b.edges.get(Id::new(e)).v.map(|v| v.0)).collect();
    verts.sort_unstable();
    verts.dedup();
    let mut splits: HashMap<usize, Vec<(f64, Id<Vertex>)>> = HashMap::new();
    for &v in &verts {
        let p = b.vertices.get(Id::new(v as usize)).p;
        for &i in grid.get(&p.map(cell)).into_iter().flatten().chain(&huge) {
            let bx = &boxes[i as usize];
            if (0..3).any(|k| p[k] < bx[k] || p[k] > bx[k + 3]) { continue; }
            let e = free[i as usize];
            let ed = b.edges.get(Id::new(e));
            if ed.v.iter().any(|x| x.0 == v) { continue; }
            let c = ed.curve.as_ref().unwrap();
            let (lo, hi) = (ed.t[0].min(ed.t[1]), ed.t[0].max(ed.t[1]));
            let mut t = t_of(c.as_ref(), p);
            if c.periodic() { while t < lo { t += std::f64::consts::TAU; } while t > hi { t -= std::f64::consts::TAU; } }
            let eps = 1e-9 * (hi - lo);
            if t <= lo + eps || t >= hi - eps || norm(sub(c.ev(t).0, p)) > tol { continue; }
            splits.entry(e).or_default().push((t, Id::new(v as usize)));
        }
    }
    // pieces in edge direction (the first reuses the edge id); coedges of the edge become coedges of the pieces
    let mut pieces: HashMap<usize, Vec<Id<Edge>>> = HashMap::new();
    let mut keys: Vec<usize> = splits.keys().copied().collect();
    keys.sort_unstable();
    for e in keys {
        let mut sp = splits.remove(&e).unwrap();
        let ed = b.edges.get(Id::new(e)).clone();
        let fwd = ed.t[1] >= ed.t[0];
        sp.sort_by(|x, y| if fwd { x.0.total_cmp(&y.0) } else { y.0.total_cmp(&x.0) });
        sp.dedup_by(|x, y| x.1 == y.1);
        let ts: Vec<f64> = std::iter::once(ed.t[0]).chain(sp.iter().map(|s| s.0)).chain([ed.t[1]]).collect();
        let vs: Vec<Id<Vertex>> = std::iter::once(ed.v[0]).chain(sp.iter().map(|s| s.1)).chain([ed.v[1]]).collect();
        let ids = (0..ts.len() - 1).map(|k| {
            let piece = Edge { curve: ed.curve.clone(), v: [vs[k], vs[k + 1]], t: [ts[k], ts[k + 1]], tol: ed.tol, origin: ed.origin };
            if k == 0 { *b.edges.get_mut(Id::new(e)) = piece; Id::new(e) } else { b.edges.push(piece) }
        }).collect();
        pieces.insert(e, ids);
    }
    for l in 0..b.loops.len() {
        let lp = b.loops.get(Id::new(l)).clone();
        if !lp.coedges.iter().any(|&c| pieces.contains_key(&b.coedges.get(c).edge.ix())) { continue; }
        let mut out = vec![];
        for &c in &lp.coedges {
            let ce = b.coedges.get(c).clone();
            let Some(ps) = pieces.get(&ce.edge.ix()) else { out.push(c); continue };
            let order: Vec<Id<Edge>> = if ce.sense { ps.clone() } else { ps.iter().rev().copied().collect() };
            for (k, pe) in order.into_iter().enumerate() {
                if k == 0 { b.coedges.get_mut(c).edge = pe; out.push(c); } else { out.push(b.coedges.push(Coedge { edge: pe, sense: ce.sense })); }
            }
        }
        b.loops.get_mut(Id::new(l)).coedges = out;
    }
    pieces.len()
}

/// reverse every loop of a face (the face normal stays)
fn reverse_loops(b: &mut Body, f: Id<Face>) {
    for &l in &b.faces.get(f).loops.clone() {
        let mut cs = b.loops.get(l).coedges.clone();
        cs.reverse();
        for &c in &cs { let ce = b.coedges.get_mut(c); ce.sense = !ce.sense; }
        b.loops.get_mut(l).coedges = cs;
    }
}

fn flip_face(b: &mut Body, f: Id<Face>) {
    reverse_loops(b, f);
    let s = b.faces.get(f).sense;
    b.faces.get_mut(f).sense = !s;
}

/// orient the faces of one shell so that each shared edge is used once in each direction, keeping the majority
/// orientation; a closed shell is then turned outward. Returns the faces flipped.
fn orient(b: &mut Body, faces: &[Id<Face>], conflicts: &mut usize) -> usize {
    let n = faces.len();
    let mut by_edge: HashMap<u32, Vec<(usize, bool)>> = HashMap::new();
    for (i, &f) in faces.iter().enumerate() { for ce in b.face_coedges(f) { by_edge.entry(ce.edge.0).or_default().push((i, ce.sense)); } }
    // (neighbour, must differ in flip, edge); same senses on a shared edge: exactly one of the faces must flip
    let mut adj: Vec<Vec<(usize, bool, u32)>> = vec![vec![]; n];
    let mut keys: Vec<u32> = by_edge.keys().copied().collect();
    keys.sort_unstable();
    for k in keys {
        let u = &by_edge[&k];
        if u.len() != 2 || u[0].0 == u[1].0 { continue; }
        let differ = u[0].1 == u[1].1;
        adj[u[0].0].push((u[1].0, differ, k));
        adj[u[1].0].push((u[0].0, differ, k));
    }
    let debug = crate::debug("sew");
    let mut flip: Vec<Option<bool>> = vec![None; n];
    for s in 0..n {
        if flip[s].is_some() { continue; }
        flip[s] = Some(false);
        let mut q = std::collections::VecDeque::from([s]);
        while let Some(i) = q.pop_front() {
            for &(j, d, ek) in &adj[i] {
                let want = flip[i].unwrap() ^ d;
                match flip[j] {
                    None => { flip[j] = Some(want); q.push_back(j); }
                    Some(x) if x != want => {
                        *conflicts += 1;
                        if debug {
                            let ed = b.edges.get(Id::new(ek as usize));
                            eprintln!("conflict faces #{} #{} edge #{}", b.faces.get(faces[i]).origin.src, b.faces.get(faces[j]).origin.src, ed.origin.src);
                        }
                    }
                    _ => {}
                }
            }
        }
    }
    let mut fl: Vec<bool> = flip.into_iter().map(Option::unwrap).collect();
    if 2 * fl.iter().filter(|x| **x).count() > n { for x in fl.iter_mut() { *x = !*x; } }
    for (i, &x) in fl.iter().enumerate() { if x { flip_face(b, faces[i]); } }
    let count = fl.iter().filter(|x| **x).count();
    if by_edge.values().all(|u| u.len() == 2) && crate::mesh::volume(b, faces).0 < 0.0 {
        for &f in faces { flip_face(b, f); }
        return n - count;
    }
    count
}

/// sew the open shells of `b` with tolerance `tol` (model units)
pub fn sew(b: &mut Body, tol: f64) -> SewStats {
    let mut st = SewStats { shells_before: b.shells.len(), ..Default::default() };
    // 1. merge free-edge end vertices within tol (grid of cell tol, 27 neighbour cells)
    let nv = b.vertices.len();
    let u = b.edge_uses();
    let mut on_free = vec![false; nv];
    for (e, ed) in b.edges.iter() { if u[e.ix()] == 1 { for v in ed.v { on_free[v.ix()] = true; } } }
    let cell = |p: V| p.map(|x| (x / tol).floor() as i64);
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    for (i, v) in b.vertices.iter() { if on_free[i.ix()] { grid.entry(cell(v.p)).or_default().push(i.0); } }
    let mut vs = Dsu::new(nv);
    for (i, v) in b.vertices.iter() {
        if !on_free[i.ix()] { continue; }
        let c = cell(v.p);
        for d in 0..27 {
            let k = [c[0] + d % 3 - 1, c[1] + d / 3 % 3 - 1, c[2] + d / 9 - 1];
            for &j in grid.get(&k).into_iter().flatten() {
                if j > i.0 && norm(sub(b.vertices.get(Id::new(j as usize)).p, v.p)) <= tol { vs.union(i.0, j); }
            }
        }
    }
    st.vertices_merged = (0..nv as u32).filter(|&i| vs.find(i) != i).count();
    for e in 0..b.edges.len() {
        let ed = b.edges.get_mut(Id::new(e));
        ed.v = ed.v.map(|v| Id::new(vs.find(v.0) as usize));
    }
    // 2. T-junctions
    st.edges_split = split_at_vertices(b, tol);
    // 3. pair coincident free edges, grouped by their end vertices
    let ne = b.edges.len();
    let u = b.edge_uses();
    st.free_before = u.iter().filter(|&&n| n == 1).count();
    let mut groups: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for e in (0..ne).filter(|&e| u[e] == 1) {
        let v = b.edges.get(Id::new(e)).v;
        groups.entry((v[0].0.min(v[1].0), v[0].0.max(v[1].0))).or_default().push(e);
    }
    // faces are joined by their original shells and by each pairing
    let mut face_of = vec![NONE; ne];
    let mut fs = Dsu::new(b.faces.len());
    for (_, sh) in b.shells.iter() {
        for &f in &sh.faces {
            for ce in b.face_coedges(f) { face_of[ce.edge.ix()] = f.0; }
            fs.union(sh.faces[0].0, f.0);
        }
    }
    let mut keys: Vec<_> = groups.keys().copied().collect();
    keys.sort_unstable();
    // edge -> (kept edge, reversed)
    let mut repl: Vec<Option<(usize, bool)>> = vec![None; ne];
    let mut taken = vec![false; ne];
    // `strict`: only edges whose faces are already joined
    let mut pair = |fs: &mut Dsu, i: usize, j: usize, strict: bool| {
        if taken[i] || taken[j] { return; }
        if strict && (face_of[i] == NONE || face_of[j] == NONE || fs.find(face_of[i]) != fs.find(face_of[j])) { return; }
        let Some(same) = coincide(b.edges.get(Id::new(i)), b.edges.get(Id::new(j)), tol) else { return };
        taken[i] = true;
        taken[j] = true;
        repl[j] = Some((i, !same));
        if face_of[i] != NONE && face_of[j] != NONE { fs.union(face_of[i], face_of[j]); }
        st.edges_merged += 1;
    };
    for g in keys.iter().map(|k| &groups[k]).filter(|g| g.len() == 2) { pair(&mut fs, g[0], g[1], false); }
    // ambiguous groups (duplicated faces): three passes within joined faces (joins grow each pass), then any
    for strict in [true, true, true, false] {
        for g in keys.iter().map(|k| &groups[k]).filter(|g| g.len() > 2) {
            for i in 0..g.len() { for j in i + 1..g.len() { pair(&mut fs, g[i], g[j], strict); } }
        }
    }
    for c in 0..b.coedges.len() {
        let ce = b.coedges.get_mut(Id::new(c));
        if let Some((k, rev)) = repl[ce.edge.ix()] { ce.edge = Id::new(k); ce.sense ^= rev; }
    }
    st.free_after = b.edge_uses().iter().filter(|&&n| n == 1).count();
    // 4. shells of faces joined by shared edges; closed shells are kept
    let mut kept: Vec<Shell> = vec![];
    let mut open_faces: Vec<Id<Face>> = vec![];
    for (_, sh) in b.shells.iter() { if sh.closed { kept.push(sh.clone()); } else { open_faces.extend(&sh.faces); } }
    let mut fs = Dsu::new(open_faces.len());
    let mut by_edge: HashMap<u32, u32> = HashMap::new();
    for (i, &f) in open_faces.iter().enumerate() {
        for ce in b.face_coedges(f) { let j = *by_edge.entry(ce.edge.0).or_insert(i as u32); fs.union(i as u32, j); }
    }
    let mut comp: HashMap<u32, usize> = HashMap::new();
    let mut shells: Vec<Shell> = vec![];
    for (i, &f) in open_faces.iter().enumerate() {
        let k = *comp.entry(fs.find(i as u32)).or_insert_with(|| { shells.push(Shell { faces: vec![], closed: false }); shells.len() - 1 });
        shells[k].faces.push(f);
    }
    // loops counter-clockwise about the face normal: positive signed UV area
    for &f in &open_faces {
        if crate::mesh::face_trim(b, f).polygons().iter().map(|p| crate::trim::area(p)).sum::<f64>() < 0.0 { reverse_loops(b, f); st.loops_reversed += 1; }
    }
    for sh in &mut shells {
        st.flipped += orient(b, &sh.faces, &mut st.conflicts);
        let mut n: HashMap<u32, u32> = HashMap::new();
        for &f in &sh.faces { for ce in b.face_coedges(f) { *n.entry(ce.edge.0).or_default() += 1; } }
        sh.closed = n.values().all(|&x| x == 2);
    }
    let mut out = Arena::default();
    for s in kept.into_iter().chain(shells) { out.push(s); }
    b.shells = out;
    st.shells_after = b.shells.len();
    st.closed_after = b.shells.iter().filter(|s| s.1.closed).count();
    st
}
