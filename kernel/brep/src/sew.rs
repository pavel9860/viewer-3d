//! Sewing: faces exported as separate open shells become shells that share their common edges.
//!
//! Only free edges (used by one coedge) and their end vertices take part; closed shells stay as they are.
//! 1. End vertices of free edges closer than the tolerance are merged.
//! 2. A free edge is split at every end vertex of another free edge lying on it within the tolerance
//!    (T-junctions: one face's long edge against its neighbour's shorter ones).
//! 3. Free edges between the same two vertices that run within the tolerance of each other (checked at a
//!    quarter, half and three quarters) are paired; the coedges of the second move to the first, reversed when it
//!    starts at the other vertex. Unambiguous pairs go first; where more edges coincide (duplicated faces) an edge
//!    pairs with one whose face is already joined to its own.
//! 4. Faces joined through shared edges form one shell. Each face's loops are turned to run counter-clockwise
//!    about its normal, then faces are flipped so that each shared edge is used once in each direction (keeping
//!    the orientation most faces had; a closed shell is turned outward). Shared edges that cannot be made
//!    consistent are counted as conflicts.
use crate::geom::*;
use crate::nurbs::t_of;
use crate::topo::*;
use std::collections::HashMap;

#[derive(Debug, Default, Clone)]
pub struct SewStats {
    pub vertices_merged: usize,
    pub edges_merged: usize,
    /// free edges split at a vertex of another free edge lying on them (T-junctions)
    pub edges_split: usize,
    pub free_before: usize,
    pub free_after: usize,
    pub shells_before: usize,
    pub shells_after: usize,
    pub closed_after: usize,
    /// faces reversed so that every shared edge is used in opposite directions by its two faces
    pub flipped: usize,
    /// faces whose loops ran against the face normal (signed UV area of all loops negative), reversed
    pub loops_reversed: usize,
    /// shared edges whose faces could not be made consistent (non-orientable joins)
    pub conflicts: usize,
}

fn find(p: &mut [u32], mut i: u32) -> u32 {
    while p[i as usize] != i { p[i as usize] = p[p[i as usize] as usize]; i = p[i as usize]; }
    i
}
fn union(p: &mut [u32], a: u32, b: u32) {
    let (a, b) = (find(p, a), find(p, b));
    if a != b { let (lo, hi) = (a.min(b), a.max(b)); p[hi as usize] = lo; }
}

/// coedge count per edge
fn uses(b: &Body) -> Vec<u32> {
    let mut u = vec![0u32; b.edges.len()];
    for (_, sh) in b.shells.iter() {
        for &f in &sh.faces { for &l in &b.faces.get(f).loops { for &c in &b.loops.get(l).coedges { u[b.coedges.get(c).edge.ix()] += 1; } } }
    }
    u
}

/// does edge `a` run within `tol` of edge `e` (between the same end points)
fn coincide(a: &Edge, e: &Edge, tol: f64) -> Option<bool> {
    let (ca, ce) = (a.curve.as_ref()?, e.curve.as_ref()?);
    for k in 1..4 {
        let p = ca.ev(a.t[0] + (a.t[1] - a.t[0]) * k as f64 / 4.0).0;
        let (lo, hi) = (e.t[0].min(e.t[1]), e.t[0].max(e.t[1]));
        let t = t_of(ce.as_ref(), p);
        // the closest point must lie on the edge, not elsewhere on its curve
        if t < lo - 1e-9 * (hi - lo).max(1.0) || t > hi + 1e-9 * (hi - lo).max(1.0) { return None; }
        if norm(sub(ce.ev(t).0, p)) > tol { return None; }
    }
    // same direction: the same start vertex; for closed edges (one vertex) the tangents at the middle
    if a.v[0] != a.v[1] { return Some(a.v[0] == e.v[0]); }
    let ta = ca.ev((a.t[0] + a.t[1]) / 2.0);
    let te = ce.ev(t_of(ce.as_ref(), ta.0));
    let (da, de) = (mulv(ta.1, (a.t[1] - a.t[0]).signum()), mulv(te.1, (e.t[1] - e.t[0]).signum()));
    Some(dot(da, de) > 0.0)
}

/// split every free edge at the vertices of other free edges that lie on it within `tol`; returns the number of
/// edges split
fn split_at_vertices(b: &mut Body, tol: f64) -> usize {
    let u = uses(b);
    let free: Vec<usize> = (0..b.edges.len()).filter(|&e| u[e] == 1 && b.edges.get(Id::new(e)).curve.is_some()).collect();
    if free.is_empty() { return 0; }
    // edge boxes from 17 samples, grown by a quarter of the longest sample step (curvature between samples) and tol
    let boxes: Vec<[f64; 6]> = free.iter().map(|&e| {
        let ed = b.edges.get(Id::new(e));
        let c = ed.curve.as_ref().unwrap();
        let pts: Vec<V> = (0..=16).map(|k| c.ev(ed.t[0] + (ed.t[1] - ed.t[0]) * k as f64 / 16.0).0).collect();
        let step = pts.windows(2).map(|w| norm(sub(w[1], w[0]))).fold(0.0, f64::max);
        let g = tol + step / 4.0;
        let mut bx = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
        for p in pts { for k in 0..3 { bx[k] = bx[k].min(p[k] - g); bx[k + 3] = bx[k + 3].max(p[k] + g); } }
        bx
    }).collect();
    let mut sizes: Vec<f64> = boxes.iter().map(|x| (0..3).map(|k| x[k + 3] - x[k]).fold(0.0, f64::max)).collect();
    sizes.sort_unstable_by(f64::total_cmp);
    let h = sizes[sizes.len() / 2].max(tol * 4.0);
    let cell = |x: f64| (x / h).floor() as i64;
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    for (i, bx) in boxes.iter().enumerate() {
        let (a, z) = ([0, 1, 2].map(|k| cell(bx[k])), [0, 1, 2].map(|k| cell(bx[k + 3])));
        // a huge edge would fill too many cells: it goes into every cell of a coarser walk only when small enough
        if (z[0] - a[0] + 1) * (z[1] - a[1] + 1) * (z[2] - a[2] + 1) > 4096 { grid.entry([i64::MIN; 3]).or_default().push(i as u32); continue; }
        for x in a[0]..=z[0] { for y in a[1]..=z[1] { for w in a[2]..=z[2] { grid.entry([x, y, w]).or_default().push(i as u32); } } }
    }
    let mut verts: Vec<u32> = free.iter().flat_map(|&e| b.edges.get(Id::new(e)).v.map(|v| v.0)).collect();
    verts.sort_unstable();
    verts.dedup();
    let mut splits: HashMap<usize, Vec<(f64, Id<Vertex>)>> = HashMap::new();
    let huge = grid.get(&[i64::MIN; 3]).cloned().unwrap_or_default();
    for &v in &verts {
        let p = b.vertices.get(Id::new(v as usize)).p;
        let cand = grid.get(&p.map(cell)).into_iter().flatten().chain(&huge);
        for &i in cand {
            let bx = &boxes[i as usize];
            if (0..3).any(|k| p[k] < bx[k] || p[k] > bx[k + 3]) { continue; }
            let e = free[i as usize];
            let ed = b.edges.get(Id::new(e));
            if ed.v[0].0 == v || ed.v[1].0 == v { continue; }
            let c = ed.curve.as_ref().unwrap();
            let (lo, hi) = (ed.t[0].min(ed.t[1]), ed.t[0].max(ed.t[1]));
            let mut t = t_of(c.as_ref(), p);
            if c.periodic() { while t < lo { t += 2.0 * std::f64::consts::PI; } while t > hi { t -= 2.0 * std::f64::consts::PI; } }
            let eps = 1e-9 * (hi - lo);
            if t <= lo + eps || t >= hi - eps || norm(sub(c.ev(t).0, p)) > tol { continue; }
            splits.entry(e).or_default().push((t, Id::new(v as usize)));
        }
    }
    // pieces in edge direction; the coedge using the edge is replaced by coedges of the pieces
    let mut pieces: HashMap<usize, Vec<Id<Edge>>> = HashMap::new();
    let mut keys: Vec<usize> = splits.keys().copied().collect();
    keys.sort_unstable();
    for e in keys {
        let mut sp = splits.remove(&e).unwrap();
        let ed = b.edges.get(Id::new(e)).clone();
        let fwd = ed.t[1] >= ed.t[0];
        sp.sort_by(|x, y| if fwd { x.0.total_cmp(&y.0) } else { y.0.total_cmp(&x.0) });
        sp.dedup_by(|x, y| x.1 == y.1);
        let mut ts = vec![ed.t[0]];
        let mut vs = vec![ed.v[0]];
        for (t, v) in sp { ts.push(t); vs.push(v); }
        ts.push(ed.t[1]);
        vs.push(ed.v[1]);
        let mut ids = vec![];
        for k in 0..ts.len() - 1 {
            let piece = Edge { curve: ed.curve.clone(), v: [vs[k], vs[k + 1]], t: [ts[k], ts[k + 1]], tol: ed.tol, origin: ed.origin };
            if k == 0 { *b.edges.get_mut(Id::new(e)) = piece; ids.push(Id::new(e)); } else { ids.push(b.edges.push(piece)); }
        }
        pieces.insert(e, ids);
    }
    let n = pieces.len();
    if n == 0 { return 0; }
    for l in 0..b.loops.len() {
        let lp = b.loops.get(Id::new(l)).clone();
        if !lp.coedges.iter().any(|&c| pieces.contains_key(&b.coedges.get(c).edge.ix())) { continue; }
        let mut out = vec![];
        for &c in &lp.coedges {
            let ce = b.coedges.get(c).clone();
            match pieces.get(&ce.edge.ix()) {
                None => out.push(c),
                Some(ps) => {
                    let order: Vec<Id<Edge>> = if ce.sense { ps.clone() } else { ps.iter().rev().copied().collect() };
                    for (k, pe) in order.into_iter().enumerate() {
                        if k == 0 { b.coedges.get_mut(c).edge = pe; out.push(c); } else { out.push(b.coedges.push(Coedge { edge: pe, sense: ce.sense })); }
                    }
                }
            }
        }
        b.loops.get_mut(Id::new(l)).coedges = out;
    }
    n
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

/// the loops of a face bound it counter-clockwise about its normal: their signed areas in the face's UV frame add up
/// to the (positive) face area. Faces whose loops run the other way get them reversed; returns how many.
fn loops_to_normal(b: &mut Body, faces: &[Id<Face>]) -> usize {
    let mut n = 0;
    for &f in faces {
        let t = crate::mesh::face_trim(b, f);
        let a: f64 = t.polygons().iter().map(|p| crate::trim::area(p)).sum();
        if a < 0.0 { reverse_loops(b, f); n += 1; }
    }
    n
}

fn flip_face(b: &mut Body, f: Id<Face>) {
    reverse_loops(b, f);
    let s = b.faces.get(f).sense;
    b.faces.get_mut(f).sense = !s;
}

/// orient the faces of one shell consistently (each shared edge used once in each direction), keeping the
/// orientation most faces had; a closed shell is then turned outward (positive volume). Returns the faces flipped.
fn orient(b: &mut Body, faces: &[Id<Face>], conflicts: &mut usize) -> usize {
    let n = faces.len();
    let mut by_edge: HashMap<u32, Vec<(usize, bool)>> = HashMap::new();
    for (i, &f) in faces.iter().enumerate() {
        for &l in &b.faces.get(f).loops { for &c in &b.loops.get(l).coedges { let ce = b.coedges.get(c); by_edge.entry(ce.edge.0).or_default().push((i, ce.sense)); } }
    }
    let mut adj: Vec<Vec<(usize, bool, u32)>> = vec![vec![]; n]; // (neighbour, must differ in flip, edge)
    let mut keys: Vec<u32> = by_edge.keys().copied().collect();
    keys.sort_unstable();
    for k in keys {
        let u = &by_edge[&k];
        if u.len() != 2 || u[0].0 == u[1].0 { continue; }
        // consistent when the senses differ; same senses: exactly one of the two faces must flip
        let differ = u[0].1 == u[1].1;
        adj[u[0].0].push((u[1].0, differ, k));
        adj[u[1].0].push((u[0].0, differ, k));
    }
    let mut flip: Vec<Option<bool>> = vec![None; n];
    for s in 0..n {
        if flip[s].is_some() { continue; }
        flip[s] = Some(false);
        let mut st = std::collections::VecDeque::from([s]);
        while let Some(i) = st.pop_front() {
            let fi = flip[i].unwrap();
            for &(j, d, ek) in &adj[i] {
                let want = fi ^ d;
                match flip[j] {
                    None => { flip[j] = Some(want); st.push_back(j); }
                    Some(x) if x != want => {
                        *conflicts += 1;
                        if std::env::var("SEW_DEBUG").is_ok() { let ed = b.edges.get(Id::new(ek as usize)); let (p0, p1) = (b.vertices.get(ed.v[0]).p, b.vertices.get(ed.v[1]).p); eprintln!("conflict faces STEP #{} #{} edge #{} {:?} -> {:?} len {:.2e}", b.faces.get(faces[i]).origin.src, b.faces.get(faces[j]).origin.src, ed.origin.src, p0, p1, norm(sub(p0, p1))); }
                    }
                    _ => {}
                }
            }
        }
    }
    let mut fl: Vec<bool> = flip.into_iter().map(|x| x.unwrap()).collect();
    if 2 * fl.iter().filter(|x| **x).count() > n { for x in fl.iter_mut() { *x = !*x; } }
    for (i, &x) in fl.iter().enumerate() { if x { flip_face(b, faces[i]); } }
    let mut count = fl.iter().filter(|x| **x).count();
    let closed = by_edge.values().all(|u| u.len() == 2);
    if closed && crate::mesh::volume(b, faces).0 < 0.0 {
        for &f in faces { flip_face(b, f); }
        count = n - count;
    }
    count
}

/// sew the open shells of `b` with tolerance `tol` (model units of the body)
pub fn sew(b: &mut Body, tol: f64) -> SewStats {
    let mut st = SewStats { shells_before: b.shells.len(), ..Default::default() };
    // vertices within tol: cells of size tol, neighbours in the 27 surrounding cells
    let nv = b.vertices.len();
    let mut par: Vec<u32> = (0..nv as u32).collect();
    let cell = |p: V| p.map(|x| (x / tol).floor() as i64);
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    // only the ends of free edges: vertices of closed shells keep their identity
    let u0 = uses(b);
    let mut on_free = vec![false; nv];
    for (e, ed) in b.edges.iter() { if u0[e.ix()] == 1 { for v in ed.v { on_free[v.ix()] = true; } } }
    for (i, v) in b.vertices.iter() { if on_free[i.ix()] { grid.entry(cell(v.p)).or_default().push(i.0); } }
    for (i, v) in b.vertices.iter() {
        if !on_free[i.ix()] { continue; }
        let c = cell(v.p);
        for dz in -1..=1 { for dy in -1..=1 { for dx in -1..=1 {
            if let Some(js) = grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) {
                for &j in js { if j > i.0 && norm(sub(b.vertices.get(Id::new(j as usize)).p, v.p)) <= tol { union(&mut par, i.0, j); } }
            }
        } } }
    }
    for i in 0..nv as u32 { if find(&mut par, i) != i { st.vertices_merged += 1; } }
    let ne = b.edges.len();
    for e in 0..ne {
        let ed = b.edges.get_mut(Id::new(e));
        ed.v = ed.v.map(|v| Id::new(find(&mut par, v.0) as usize));
    }
    st.edges_split = split_at_vertices(b, tol);
    let ne = b.edges.len();
    // free edges grouped by their end vertices
    let u = uses(b);
    st.free_before = u.iter().filter(|&&n| n == 1).count();
    let mut groups: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for e in 0..ne {
        if u[e] != 1 { continue; }
        let v = b.edges.get(Id::new(e)).v;
        groups.entry((v[0].0.min(v[1].0), v[0].0.max(v[1].0))).or_default().push(e);
    }
    // replacement: edge -> (kept edge, reversed)
    let mut repl: Vec<Option<(usize, bool)>> = vec![None; ne];
    // faces of each free edge; faces are joined (union-find) by their original shells and by each pairing.
    // Unambiguous pairs (two free edges between the same vertices) go first; where more edges coincide (e.g.
    // duplicated faces), an edge pairs with one whose face is already joined to its own, so that copies are sewn
    // to themselves and not to each other; left-over ambiguous edges pair in order at the end.
    let mut face_of = vec![u32::MAX; ne];
    let mut fpar: Vec<u32> = (0..b.faces.len() as u32).collect();
    for (_, sh) in b.shells.iter() {
        let mut first = None;
        for &f in &sh.faces {
            for &l in &b.faces.get(f).loops { for &c in &b.loops.get(l).coedges { face_of[b.coedges.get(c).edge.ix()] = f.0; } }
            match first { None => first = Some(f.0), Some(r) => union(&mut fpar, r, f.0) }
        }
    }
    let mut keys: Vec<_> = groups.keys().copied().collect();
    keys.sort_unstable();
    let mut taken = vec![false; ne];
    let pair = |b: &Body, fpar: &mut Vec<u32>, taken: &mut Vec<bool>, repl: &mut Vec<Option<(usize, bool)>>, i: usize, j: usize, same: bool| {
        taken[i] = true;
        taken[j] = true;
        repl[j] = Some((i, !same));
        if face_of[i] != u32::MAX && face_of[j] != u32::MAX { union(fpar, face_of[i], face_of[j]); }
        let _ = b;
    };
    for &k in &keys {
        let g = &groups[&k];
        if g.len() != 2 { continue; }
        if let Some(same) = coincide(b.edges.get(Id::new(g[0])), b.edges.get(Id::new(g[1])), tol) { pair(b, &mut fpar, &mut taken, &mut repl, g[0], g[1], same); st.edges_merged += 1; }
    }
    for strict in [true, true, true, false] {
        for &k in &keys {
            let g = &groups[&k];
            if g.len() <= 2 { continue; }
            for i in 0..g.len() {
                if taken[g[i]] { continue; }
                for j in i + 1..g.len() {
                    if taken[g[j]] || taken[g[i]] { continue; }
                    let joined = face_of[g[i]] != u32::MAX && face_of[g[j]] != u32::MAX && find(&mut fpar, face_of[g[i]]) == find(&mut fpar, face_of[g[j]]);
                    if strict && !joined { continue; }
                    if let Some(same) = coincide(b.edges.get(Id::new(g[i])), b.edges.get(Id::new(g[j])), tol) { pair(b, &mut fpar, &mut taken, &mut repl, g[i], g[j], same); st.edges_merged += 1; }
                }
            }
        }
    }
    let nc = b.coedges.len();
    for c in 0..nc {
        let ce = b.coedges.get_mut(Id::new(c));
        if let Some((k, rev)) = repl[ce.edge.ix()] { ce.edge = Id::new(k); if rev { ce.sense = !ce.sense; } }
    }
    // shells: faces joined by shared edges; closed shells are kept as they are
    let u = uses(b);
    st.free_after = u.iter().filter(|&&n| n == 1).count();
    let mut open_faces: Vec<Id<Face>> = vec![];
    let mut kept: Vec<Shell> = vec![];
    for (_, sh) in b.shells.iter() { if sh.closed { kept.push(sh.clone()); } else { open_faces.extend(&sh.faces); } }
    let mut fpar: Vec<u32> = (0..open_faces.len() as u32).collect();
    let mut by_edge: HashMap<u32, u32> = HashMap::new();
    for (i, &f) in open_faces.iter().enumerate() {
        for &l in &b.faces.get(f).loops {
            for &c in &b.loops.get(l).coedges {
                let e = b.coedges.get(c).edge.0;
                match by_edge.get(&e) { Some(&j) => union(&mut fpar, i as u32, j), None => { by_edge.insert(e, i as u32); } }
            }
        }
    }
    let mut comp: HashMap<u32, usize> = HashMap::new();
    let mut shells: Vec<Shell> = vec![];
    for (i, &f) in open_faces.iter().enumerate() {
        let r = find(&mut fpar, i as u32);
        let k = *comp.entry(r).or_insert_with(|| { shells.push(Shell { faces: vec![], closed: false }); shells.len() - 1 });
        shells[k].faces.push(f);
    }
    st.loops_reversed = loops_to_normal(b, &open_faces);
    for sh in &mut shells { st.flipped += orient(b, &sh.faces, &mut st.conflicts); }
    for sh in &mut shells {
        let mut n: HashMap<u32, u32> = HashMap::new();
        for &f in &sh.faces { for &l in &b.faces.get(f).loops { for &c in &b.loops.get(l).coedges { *n.entry(b.coedges.get(c).edge.0).or_default() += 1; } } }
        sh.closed = n.values().all(|&x| x == 2);
    }
    let mut out = Arena::default();
    for s in kept.into_iter().chain(shells) { out.push(s); }
    b.shells = out;
    st.shells_after = b.shells.len();
    st.closed_after = b.shells.iter().filter(|s| s.1.closed).count();
    st
}
