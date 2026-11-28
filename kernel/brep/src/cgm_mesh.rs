//! Triangle mesh of a CGM model (`cgm::Model`).
//!
//! Each edge is sampled once, along its first coedge's p-curve, until the chord deviation is below
//! `sag`; every face using the edge reuses those 3D points (watertight by construction) and maps the
//! samples to its own p-curve by the coedge's parameter range. A face's loops become UV polygons
//! (outer loop counter-clockwise, holes clockwise: CGM stores each loop's senses up to one direction
//! for the whole loop), are triangulated (`tri`) and refined to `sag` by `refine` (the STEP mesher's: Delaunay
//! flips in the surface metric, longest-edge bisection, so facets stay tangent to the surface).
use crate::cgm::*;
use crate::geom::*;
use crate::exact::orient2d;
use crate::trim::{area, P2};
use std::collections::HashMap;

pub struct Mesh {
    pub p: Vec<V>,
    pub tris: Vec<[u32; 3]>,
    /// per triangle: index into `faces`
    pub tri_face: Vec<u32>,
    /// per meshed face: (CGM id, face sense, shell sense, outer loop direction, outer loop UV area sign as chained)
    pub faces: Vec<(u32, i32, i32, i32, i32)>,
    /// faces that could not be triangulated, with the reason
    pub failed: Vec<(u32, String)>,
    /// per sampled B-rep edge: the edge index and its mesh points in order
    pub edge_pts: Vec<(usize, Vec<u32>)>,
    /// per sampled B-rep edge: how many coedges of the meshed faces use it
    pub edge_uses: HashMap<usize, u32>,
    /// faces whose refinement ran out of passes or triangles before meeting the chord (meshed, coarser than asked),
    /// with the largest deviation left
    pub coarse: Vec<(u32, f64)>,
}

/// open mesh edges by cause: on a B-rep edge used by one meshed coedge (its other face is not in the mesh), on a
/// B-rep edge used by two or more (the tessellation does not join), or not on a B-rep edge (inside a face)
#[derive(Debug, Default)]
pub struct OpenCauses { pub lone: usize, pub joined: usize, pub inside: usize, pub lone_edges: Vec<usize>, pub joined_edges: Vec<usize>, pub inside_faces: Vec<u32> }

pub fn open_causes(m: &Mesh) -> OpenCauses {
    let mut uses: HashMap<(u32, u32), i32> = HashMap::new();
    let mut face_of: HashMap<(u32, u32), u32> = HashMap::new();
    for (ti, t) in m.tris.iter().enumerate() { for k in 0..3 { let (i, j) = (t[k], t[(k + 1) % 3]); *uses.entry((i.min(j), i.max(j))).or_default() += 1; face_of.insert((i.min(j), i.max(j)), m.faces[m.tri_face[ti] as usize].0); } }
    let mut on: HashMap<(u32, u32), usize> = HashMap::new();
    for (e, pts) in &m.edge_pts { for w in pts.windows(2) { on.insert((w[0].min(w[1]), w[0].max(w[1])), *e); } }
    let mut c = OpenCauses::default();
    let mut keys: Vec<_> = uses.iter().filter(|x| *x.1 == 1).map(|x| *x.0).collect();
    keys.sort_unstable();
    for k in keys {
        match on.get(&k) {
            None => { c.inside += 1; c.inside_faces.push(face_of[&k]); }
            Some(&e) if m.edge_uses.get(&e).copied().unwrap_or(0) < 2 => { c.lone += 1; c.lone_edges.push(e); }
            Some(&e) => { c.joined += 1; c.joined_edges.push(e); }
        }
    }
    c.lone_edges.sort_unstable();
    c.lone_edges.dedup();
    c.joined_edges.sort_unstable();
    c.joined_edges.dedup();
    c.inside_faces.sort_unstable();
    c.inside_faces.dedup();
    c
}

/// `q` within the domain `d` up to 1e-3 of its width
fn inside(q: P2, d: [f64; 4]) -> bool {
    let (mu, mv) = (1e-3 * (d[1] - d[0]), 1e-3 * (d[3] - d[2]));
    q[0] >= d[0] - mu && q[0] <= d[1] + mu && q[1] >= d[2] - mv && q[1] <= d[3] + mv
}

/// largest distance between the end points of the two coedges of an edge (each evaluated through its p-curve and
/// surface), per pair of (surface kind, p-curve kind) of the two coedges: (edges, largest gap)
pub fn edge_gaps(m: &Model) -> std::collections::BTreeMap<[(u32, u32); 2], (usize, f64)> {
    let mut by_edge: HashMap<usize, Vec<usize>> = HashMap::new();
    for f in &m.faces { for l in &f.loops { for c in l.clone() { by_edge.entry(m.coedges[c].edge).or_default().push(c); } } }
    let face_of: HashMap<usize, usize> = m.faces.iter().flat_map(|f| f.loops.iter().flat_map(|l| l.clone()).map(move |c| (c, f.surface))).collect();
    let mut out = std::collections::BTreeMap::new();
    for cs in by_edge.values().filter(|cs| cs.len() == 2) {
        let ends = |c: usize| { let ce = &m.coedges[c]; [m.point(ce, ce.t[0]), m.point(ce, ce.t[1])] };
        let (a, b) = (ends(cs[0]), ends(cs[1]));
        let d = |x: V, y: V| norm(sub(x, y));
        let gap = d(a[0], b[0]).max(d(a[1], b[1])).min(d(a[0], b[1]).max(d(a[1], b[0])));
        let kind = |c: usize| (m.surface_kinds[face_of[&c]], m.pcurve_kinds[m.coedges[c].pcurve]);
        let (x, y) = (kind(cs[0]), kind(cs[1]));
        let e: &mut (usize, f64) = out.entry([x.min(y), x.max(y)]).or_default();
        e.0 += 1;
        e.1 = e.1.max(gap);
    }
    out
}

/// closure and size measures of a mesh
#[derive(Debug, Default)]
pub struct Stats { pub tris: usize, pub area: f64, pub volume: f64, pub open_edges: usize, pub nonmanifold_edges: usize, pub bbox: [f64; 6] }

pub fn stats(m: &Mesh) -> Stats {
    let mut s = Stats { tris: m.tris.len(), bbox: [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN], ..Default::default() };
    let mut uses: HashMap<(u32, u32), i32> = HashMap::new();
    for t in &m.tris {
        let (a, b, c) = (m.p[t[0] as usize], m.p[t[1] as usize], m.p[t[2] as usize]);
        s.area += norm(cross(sub(b, a), sub(c, a))) / 2.0;
        s.volume += dot(a, cross(b, c)) / 6.0;
        for k in 0..3 {
            let (i, j) = (t[k], t[(k + 1) % 3]);
            *uses.entry((i.min(j), i.max(j))).or_default() += 1;
        }
    }
    for &q in &m.p { for k in 0..3 { s.bbox[k] = s.bbox[k].min(q[k]); s.bbox[k + 3] = s.bbox[k + 3].max(q[k]); } }
    s.open_edges = uses.values().filter(|&&n| n == 1).count();
    s.nonmanifold_edges = uses.values().filter(|&&n| n > 2).count();
    s
}

/// per meshed face, whether it must be flipped so that the mesh is consistently oriented with
/// positive volume (None when inconsistent within the face or the mesh has open/shared edges)
pub fn orientation_fix(m: &Mesh) -> Option<Vec<bool>> {
    let nf = m.faces.len();
    // face adjacency through directed edges: same direction in two faces = opposite orientation
    let mut dir: HashMap<(u32, u32), Vec<(u32, bool)>> = HashMap::new();
    for (t, f) in m.tris.iter().zip(&m.tri_face) {
        for k in 0..3 { let (a, b) = (t[k], t[(k + 1) % 3]); dir.entry((a.min(b), a.max(b))).or_default().push((*f, a < b)); }
    }
    let mut adj: Vec<Vec<(usize, bool)>> = vec![vec![]; nf]; // (other face, same relative orientation)
    for uses in dir.values() {
        if uses.len() != 2 { return None; }
        let ((f, d), (g, e)) = (uses[0], uses[1]);
        if f == g { if d == e { return None; } continue; }
        adj[f as usize].push((g as usize, d != e));
    }
    let mut flip: Vec<Option<bool>> = vec![None; nf];
    for s in 0..nf {
        if flip[s].is_some() { continue; }
        flip[s] = Some(false);
        let mut stack = vec![s];
        let mut comp = vec![s];
        while let Some(f) = stack.pop() {
            for &(g, same) in &adj[f] {
                let want = if same { flip[f].unwrap() } else { !flip[f].unwrap() };
                match flip[g] { None => { flip[g] = Some(want); stack.push(g); comp.push(g); } Some(x) => if x != want { return None; } }
            }
        }
        // the component's volume sign decides its global flip
        let vol: f64 = m.tris.iter().zip(&m.tri_face).filter(|(_, f)| comp.contains(&(**f as usize))).map(|(t, f)| {
            let (a, b, c) = (m.p[t[0] as usize], m.p[t[1] as usize], m.p[t[2] as usize]);
            let v = dot(a, cross(b, c)) / 6.0;
            if flip[*f as usize].unwrap() { -v } else { v }
        }).sum();
        if vol < 0.0 { for &f in &comp { flip[f] = Some(!flip[f].unwrap()); } }
    }
    Some(flip.into_iter().map(|x| x.unwrap()).collect())
}

/// edge samples: parameters as fractions of the first coedge's range, and mesh point indices
struct EdgeSamples { f: Vec<f64>, idx: Vec<u32> }

/// the shells a part shows: those of the bodies in the first two lists of its CATSeeBodyCont
/// (solids, and sheets when `sheets`), or every solid when the part lists none. Corpus check against
/// CATIA's stored tessellation: the shown body is in list 0 (5-list form, with the origin planes) or
/// list 1 (3-list form; list 2 then holds the origin planes).
pub fn shown_shells(m: &Model, sheets: bool) -> Vec<usize> {
    let listed: Vec<&Body> = m.bodies.iter().filter(|b| b.flags & 0b110 != 0).collect();
    let mut out: Vec<usize> = if listed.is_empty() { (0..m.shells.len()).filter(|&i| m.shells[i].solid).collect() } else { listed.iter().flat_map(|b| b.shells.iter().copied()).filter(|&i| m.shells[i].solid || sheets).collect() };
    out.sort();
    out.dedup();
    out
}

/// bounding box of the shells' edges (each coedge's ends and middle): the part's size, for a tolerance relative to
/// it, without meshing
pub fn edge_bbox(m: &Model, shells: &[usize]) -> [f64; 6] {
    let mut b = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
    for &sh in shells {
        for &(f, _) in &m.shells[sh].faces {
            for c in m.faces[f].loops.iter().flat_map(|l| l.clone()) {
                let ce = &m.coedges[c];
                if m.bad_pcurve[ce.pcurve].is_some() || m.bad_surface[m.faces[f].surface].is_some() { continue; }
                for t in [ce.t[0], 0.5 * (ce.t[0] + ce.t[1]), ce.t[1]] {
                    let p = m.point(ce, t);
                    if p.iter().all(|x| x.is_finite()) { for k in 0..3 { b[k] = b[k].min(p[k]); b[k + 3] = b[k + 3].max(p[k]); } }
                }
            }
        }
    }
    b
}

/// meshes the faces of the given shells on `workers` threads (the result does not depend on their number)
pub fn tessellate(m: &Model, sag: f64, shells: &[usize], workers: usize) -> Mesh {
    let timing = crate::debug("time");
    let t_all = std::time::Instant::now();
    let mut out = Mesh { p: vec![], tris: vec![], tri_face: vec![], faces: vec![], failed: vec![], edge_pts: vec![], edge_uses: HashMap::new(), coarse: vec![] };
    let mut vtx: Vec<Option<u32>> = vec![None; m.vertices];
    let mut edges: HashMap<usize, EdgeSamples> = HashMap::new();
    // faces of the shells, each once
    let mut seen = vec![false; m.faces.len()];
    let faces: Vec<(usize, i32)> = shells.iter().flat_map(|&i| m.shells[i].faces.iter().copied()).filter(|&(f, _)| !std::mem::replace(&mut seen[f], true)).collect();
    let coedges: Vec<usize> = faces.iter().flat_map(|&(f, _)| m.faces[f].loops.iter().flat_map(|l| l.clone())).collect();
    // each edge is sampled along its first coedge whose p-curve stays in its surface's domain (ends and middle)
    let surf_of: HashMap<usize, usize> = faces.iter().flat_map(|&(f, _)| m.faces[f].loops.iter().flat_map(|l| l.clone()).map(move |c| (c, m.faces[f].surface))).collect();
    let in_dom = |c: usize| {
        let ce = &m.coedges[c];
        if m.bad_pcurve[ce.pcurve].is_some() || m.bad_surface[surf_of[&c]].is_some() { return false; }
        let d = m.surfaces[surf_of[&c]].dom();
        let pc = &m.pcurves[ce.pcurve].1;
        [ce.t[0], (ce.t[0] + ce.t[1]) / 2.0, ce.t[1]].iter().all(|&t| inside(pc.uv(t), d))
    };
    let mut along: HashMap<usize, usize> = HashMap::new();
    for &c in &coedges {
        let e = m.coedges[c].edge;
        match along.get(&e) { Some(&k) if in_dom(k) || !in_dom(c) => {} _ => { along.insert(e, c); } }
    }
    // edge samples as fractions of the sampling coedge's range; refined where a face's UV rings cross
    let mut fr: HashMap<usize, Vec<f64>> = HashMap::new();
    for &c in &coedges {
        *out.edge_uses.entry(m.coedges[c].edge).or_default() += 1;
        let ce = &m.coedges[along[&m.coedges[c].edge]];
        fr.entry(ce.edge).or_insert_with(|| sample(|t| m.point(ce, t), ce.t[0], ce.t[1], sag));
    }
    let cap: HashMap<usize, usize> = fr.iter().map(|(&e, f)| (e, 4 * f.len() + 64)).collect();
    let valid = |fi: usize| { let f = &m.faces[fi]; m.bad_surface[f.surface].is_none() && f.loops.iter().flat_map(|l| l.clone()).all(|c| m.bad_pcurve[m.coedges[c].pcurve].is_none()) };
    // a face keeps refining only while its crossings decrease (crossings of the curves themselves, i.e. bad loops,
    // do not go away)
    let mut last: HashMap<usize, usize> = HashMap::new();
    for _round in 0..10 {
        let mut split: HashMap<usize, Vec<usize>> = HashMap::new();
        let todo: Vec<usize> = faces.iter().map(|f| f.0).filter(|&fi| valid(fi) && last.get(&fi) != Some(&0)).collect();
        let found = crate::par::map(todo.len(), workers, |i| face_rings(m, todo[i], &fr, &along).ok().map(|r| crossing_segments(&r)));
        for (&fi, x) in todo.iter().zip(found) {
            let Some(x) = x else { continue };
            let before = last.insert(fi, x.len());
            if x.is_empty() || before.is_some_and(|b| x.len() >= b) { last.insert(fi, 0); continue; }
            for (e, k) in x { split.entry(e).or_default().push(k); }
        }
        if split.is_empty() { break; }
        let mut grew = false;
        for (e, mut ks) in split {
            ks.sort_unstable();
            ks.dedup();
            let f = fr.get_mut(&e).unwrap();
            // loops that really cross (bad input) cannot be cured by sampling: at most 4x the initial samples + 64
            if f.len() + ks.len() > cap[&e] { continue; }
            grew = true;
            for &k in ks.iter().rev() { let x = 0.5 * (f[k] + f[k + 1]); f.insert(k + 1, x); }
        }
        if !grew { break; }
    }
    // mesh points: one per sample, shared by every face using the edge; vertices shared across edges
    let mut keys: Vec<usize> = fr.keys().copied().collect();
    keys.sort_unstable();
    for e in keys {
        let f = fr.remove(&e).unwrap();
        let ce = &m.coedges[along[&e]];
        let ends = m.edges[e];
        let mut idx = Vec::with_capacity(f.len());
        for (k, &x) in f.iter().enumerate() {
            let at_end = if k == 0 { Some(ends[0]) } else if k + 1 == f.len() { Some(ends[1]) } else { None };
            let new = |out: &mut Mesh| { out.p.push(m.point(ce, ce.t[0] + x * (ce.t[1] - ce.t[0]))); (out.p.len() - 1) as u32 };
            idx.push(match at_end {
                Some(v) => *vtx[v].get_or_insert_with(|| new(&mut out)),
                None => new(&mut out),
            });
        }
        out.edge_pts.push((e, idx.clone()));
        edges.insert(e, EdgeSamples { f, idx });
    }
    let fr: HashMap<usize, Vec<f64>> = edges.iter().map(|(&e, es)| (e, es.f.clone())).collect();
    // B-rep edges each edge-sample mesh point lies on (a vertex: every edge meeting there)
    let mut on_edge: HashMap<u32, Vec<usize>> = HashMap::new();
    for (e, idx) in &out.edge_pts { for &i in idx { on_edge.entry(i).or_default().push(*e); } }
    let ropts = crate::refine::Opts { chord: sag, max_passes: 16, max_tris: 200_000 };
    let t_edges = t_all.elapsed();
    // faces are independent once the edge points are fixed: meshed in parallel, merged in face order
    let cx = FaceCtx { m, sag, fr: &fr, along: &along, edges: &edges, p: &out.p, on_edge: &on_edge, ropts };
    let faces_out = crate::par::map(faces.len(), workers, |i| mesh_face(&cx, faces[i].0, faces[i].1));
    let (mut t_rings, mut t_tri, mut t_ref) = (0.0, 0.0, 0.0);
    let mut slow: Vec<(f64, u32, usize)> = vec![];
    for fo in faces_out {
        t_rings += fo.t[0];
        t_tri += fo.t[1];
        t_ref += fo.t[2];
        out.failed.extend(fo.failed.into_iter().map(|w| (fo.id, w)));
        if let Some(dev) = fo.coarse { out.coarse.push((fo.id, dev)); }
        let Some((rec, index, new, tri, flip)) = fo.mesh else { continue };
        if timing { slow.push((fo.t.iter().sum(), fo.id, tri.len())); }
        let base = out.p.len() as u32;
        out.p.extend(new);
        let n0 = index.len();
        let at = |i: usize| if i < n0 { index[i] } else { base + (i - n0) as u32 };
        out.faces.push(rec);
        for t in tri {
            let (a, b, c) = (at(t[0]), at(t[1]), at(t[2]));
            // a triangle on a pole segment has the pole twice: no area, no part of the surface
            if a == b || b == c || c == a { continue; }
            out.tris.push(if flip { [a, c, b] } else { [a, b, c] });
            out.tri_face.push((out.faces.len() - 1) as u32);
        }
    }
    if timing {
        slow.sort_by(|x, y| y.0.total_cmp(&x.0));
        eprintln!("time: edges {:.1} ms, rings {:.1} ms, triangulation {:.1} ms, refinement {:.1} ms, total {:.1} ms; slowest faces (ms, id, triangles): {:?}",
            t_edges.as_secs_f64() * 1e3, t_rings * 1e3, t_tri * 1e3, t_ref * 1e3, t_all.elapsed().as_secs_f64() * 1e3,
            slow.iter().take(5).map(|x| (format!("{:.1}", x.0 * 1e3), x.1, x.2)).collect::<Vec<_>>());
    }
    out
}



/// what `mesh_face` reads: the model, the edge samples and their mesh points
struct FaceCtx<'a> {
    m: &'a Model,
    sag: f64,
    fr: &'a HashMap<usize, Vec<f64>>,
    along: &'a HashMap<usize, usize>,
    edges: &'a HashMap<usize, EdgeSamples>,
    /// mesh points of the edge samples
    p: &'a [V],
    on_edge: &'a HashMap<u32, Vec<usize>>,
    ropts: crate::refine::Opts,
}

/// one face's mesh: failure reasons, the deviation left when refinement ran out, seconds per phase (rings,
/// triangulation, refinement), and the mesh itself: the face record, global indices of the local points that are
/// edge points, the new interior points (local indices after those), triangles and whether to reverse them
struct FaceOut {
    id: u32,
    failed: Vec<String>,
    coarse: Option<f64>,
    t: [f64; 3],
    mesh: Option<((u32, i32, i32, i32, i32), Vec<u32>, Vec<V>, Vec<[usize; 3]>, bool)>,
}

fn mesh_face(cx: &FaceCtx, fi: usize, shell_sense: i32) -> FaceOut {
    let (m, sag) = (cx.m, cx.sag);
    let face = &m.faces[fi];
    let mut fo = FaceOut { id: face.id, failed: vec![], coarse: None, t: [0.0; 3], mesh: None };
    let fail = |mut fo: FaceOut, why: String| { fo.failed.push(why); fo };
    if let Some(why) = &m.bad_surface[face.surface] { return fail(fo, format!("surface kind {}: {why}", m.surface_kinds[face.surface])); }
    if let Some(c) = face.loops.iter().flat_map(|l| l.clone()).find(|&c| m.bad_pcurve[m.coedges[c].pcurve].is_some()) {
        let pc = m.coedges[c].pcurve;
        return fail(fo, format!("p-curve kind {}: {}", m.pcurve_kinds[pc], m.bad_pcurve[pc].as_ref().unwrap()));
    }
    let s = m.surfaces[face.surface].as_ref();
    // UV rings with their mesh point indices
    let t0 = std::time::Instant::now();
    let (mut rings, dirs) = match face_rings(m, fi, cx.fr, cx.along) {
        Ok(r) => {
            let dirs: Vec<i32> = r.iter().map(|x| x.2).collect();
            (r.into_iter().filter(|x| x.0.len() >= 3).map(|(uv, refs, _)| (uv, refs.iter().map(|&(e, k, _)| cx.edges[&e].idx[k]).collect::<Vec<u32>>())).collect::<Vec<_>>(), dirs)
        }
        Err(why) => return fail(fo, why),
    };
    // boundary within the surface's parameter domain (surfaces with finite data: knots, profile ranges), with
    // a margin of 1e-3 of its width
    let d = s.dom();
    if let Some(q) = rings.iter().flat_map(|r| &r.0).find(|&&q| !inside(q, d)) {
        return fail(fo, format!("boundary point ({:.6}, {:.6}) outside the surface domain {d:?}", q[0], q[1]));
    }
    if rings.is_empty() { return fail(fo, "no loop with 3 points".into()); }
    // the face's own p-curves must put its boundary on the shared edge points (within the chord tolerance)
    let off = rings.iter().flat_map(|r| r.0.iter().zip(&r.1)).map(|(q, &i)| norm(sub(s.ev(q[0], q[1]).p, cx.p[i as usize]))).fold(0.0, f64::max);
    if off > sag { return fail(fo, format!("boundary {off:.3e} off the shared edge points (p-curves disagree with the neighbouring faces)")); }
    // outer loop: largest area, counter-clockwise; the others clockwise
    let outer = (0..rings.len()).max_by(|&a, &b| area(&rings[a].0).abs().total_cmp(&area(&rings[b].0).abs())).unwrap();
    let chained = if area(&rings[outer].0) > 0.0 { 1 } else { -1 };
    for (i, r) in rings.iter_mut().enumerate() {
        let a = area(&r.0);
        if (i == outer) != (a > 0.0) { r.0.reverse(); r.1.reverse(); }
    }
    // outer ring first, then the holes
    rings.swap(0, outer);
    let uv_rings: Vec<Vec<P2>> = rings.iter().map(|r| r.0.clone()).collect();
    fo.t[0] = t0.elapsed().as_secs_f64();
    let t0 = std::time::Instant::now();
    let (pts, src, tris) = crate::tri::triangulate_src(&uv_rings);
    fo.t[1] = t0.elapsed().as_secs_f64();
    let index: Vec<u32> = src.iter().map(|&(r, k)| rings[r].1[k]).collect();
    let target = area(&rings[0].0).abs() - rings[1..].iter().map(|h| area(&h.0).abs()).sum::<f64>();
    let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]]).abs()).sum();
    let covered = (got - target).abs() <= 1e-6 * target.max(1e-300);
    if !covered { fo.failed.push(format!("triangulation covered {got:.6e} of {target:.6e} (UV area)")); }
    // every ring segment between distinct points must be a triangle edge
    let mut tedges: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for t in &tris { for k in 0..3 { let (x, y) = (index[t[k]], index[t[(k + 1) % 3]]); tedges.insert((x.min(y), x.max(y))); } }
    let lost = rings.iter().flat_map(|r| (0..r.1.len()).map(move |k| (r.1[k], r.1[(k + 1) % r.1.len()]))).filter(|&(x, y)| x != y && cx.p[x as usize] != cx.p[y as usize] && !tedges.contains(&(x.min(y), x.max(y)))).count();
    let dump = std::env::var("BREP_FACE").ok().and_then(|x| x.parse::<u32>().ok()) == Some(face.id);
    // BREP_FACE=ID: the face's UV rings and triangles
    if dump {
        for (ri, r) in rings.iter().enumerate() { for (k, q) in r.0.iter().enumerate() { eprintln!("ring {ri} {k} idx {} uv {:?} {:?}", r.1[k], q[0], q[1]); } }
        for t in &tris { eprintln!("tri {:?} idx {:?}", t, t.map(|i| index[i])); }
    }
    if lost > 0 { fo.failed.push(format!("{lost} boundary segments are not triangle edges")); }
    let mut uvs = pts;
    let mut tri: Vec<[usize; 3]> = tris;
    // refined to the chord with the STEP mesher's refinement (metric Delaunay flips, longest-edge bisection). Its
    // input must be a conforming triangulation of the face; one that overlaps or lost boundary segments (already
    // reported above) is emitted as triangulated
    let mut pos: Vec<V> = index.iter().map(|&i| cx.p[i as usize]).collect();
    let same_edge = |a: usize, c: usize| a < index.len() && c < index.len() && cx.on_edge.get(&index[a]).zip(cx.on_edge.get(&index[c])).is_some_and(|(x, y)| x.iter().any(|e| y.contains(e)));
    let n0 = uvs.len();
    if covered && lost == 0 {
        let t0 = std::time::Instant::now();
        let r = crate::refine::refine(&|q: P2| s.ev(q[0], q[1]), &mut uvs, &mut pos, &mut tri, &same_edge, &cx.ropts);
        fo.t[2] = t0.elapsed().as_secs_f64();
        if r.coarse { fo.coarse = Some(r.dev); }
    }
    // BREP_FACE=ID: also the refined triangles in UV
    if dump { for t in &tri { eprintln!("refined {:?} {:?} {:?}", uvs[t[0]], uvs[t[1]], uvs[t[2]]); } }
    // UV counter-clockwise gives the normal Su x Sv. A shell's sense +1 puts the material on the
    // surface-normal side, so the outward normal is -sense * (Su x Sv) (2,930 of 2,934 faces of
    // the closed meshes of the corpus; the face record's own i16 does not enter)
    pos.drain(..n0);
    fo.mesh = Some(((face.id, face.sense, shell_sense, dirs[outer], chained), index, pos, tri, shell_sense > 0));
    fo
}

/// UV rings of a face from edge sample fractions `fr`, each point with its (edge, sample index, interval of the
/// segment starting there), and the loop direction. On a periodic surface each coedge is shifted by whole periods to start where the previous one ended.
fn face_rings(m: &Model, fi: usize, fr: &HashMap<usize, Vec<f64>>, along: &HashMap<usize, usize>) -> Result<Vec<(Vec<P2>, Vec<(usize, usize, usize)>, i32)>, String> {
    let face = &m.faces[fi];
    // two parameter points of one vertex are a pole segment where the surface is singular there (Su x Sv vanishes:
    // a sphere's pole, a cone's or a spindle torus's apex) and they differ; elsewhere they differ by rounding only
    let surf = m.surfaces[face.surface].as_ref();
    let pole = |a: P2, b: P2| {
        if a == b { return false; }
        let d = surf.ev(a[0], a[1]);
        let (n, s2) = (norm(cross(d.u, d.v)), dot(d.u, d.u) + dot(d.v, d.v));
        n <= 1e-6 * s2
    };
    let per = m.surfaces[face.surface].periods();
    let mut rings = vec![];
    for l in &face.loops {
        let dir = loop_direction(m, l.clone()).ok_or("loop does not chain in either direction")?;
        let (mut uv, mut refs) = (vec![], vec![]);
        let mut last: Option<P2> = None;
        // the previous coedge's last point: dropped where the next coedge starts at it, kept where the next starts at
        // the vertex's other parameter point on a pole, so the ring keeps the pole segment
        let mut tail: Option<(P2, (usize, usize, usize))> = None;
        for c in l.clone() {
            let ce = &m.coedges[c];
            let f = &fr[&ce.edge];
            let pc = &m.pcurves[ce.pcurve].1;
            let n = f.len();
            let order: Vec<usize> = if ce.sense * dir > 0 { (0..n).collect() } else { (0..n).rev().collect() };
            // parameters on this coedge of the edge's sample points: the sampling coedge's own, else the projection
            // of the sample onto this coedge (coedges of one edge need not share a parametrization)
            let a = &m.coedges[along[&ce.edge]];
            let t_of = |x: f64| -> f64 {
                let t0 = ce.t[0] + x * (ce.t[1] - ce.t[0]);
                if along[&ce.edge] == c { return t0; }
                let target = m.point(a, a.t[0] + x * (a.t[1] - a.t[0]));
                // safeguarded Newton on the distance: t stays in the coedge's range and a step is taken only when it
                // brings the point closer (else halved), so an ill-conditioned projection (near a pole) cannot walk
                // the p-curve out of its range
                let (lo, hi) = (ce.t[0].min(ce.t[1]), ce.t[0].max(ce.t[1]));
                let h = 1e-7 * (hi - lo).max(1e-12);
                let dist = |t: f64| { let q = sub(m.point(ce, t), target); dot(q, q) };
                let (mut t, mut f) = (t0, dist(t0));
                for _ in 0..16 {
                    let q = sub(m.point(ce, t), target);
                    let d = mulv(sub(m.point(ce, t + h), m.point(ce, t - h)), 0.5 / h);
                    let dd = dot(d, d);
                    if dd <= 0.0 { break; }
                    let mut step = dot(q, d) / dd;
                    let mut moved = false;
                    for _ in 0..30 {
                        let tn = (t - step).clamp(lo, hi);
                        let fnew = dist(tn);
                        if fnew < f { t = tn; f = fnew; moved = true; break; }
                        step *= 0.5;
                    }
                    if !moved || step.abs() <= 1e-12 * (hi - lo) { break; }
                }
                t
            };
            let mut pts: Vec<P2> = order.iter().map(|&k| pc.uv(if k == 0 { ce.t[0] } else if k + 1 == n { ce.t[1] } else { t_of(f[k]) })).collect();
            if let Some(e) = last {
                for a in 0..2 {
                    if per[a] <= 0.0 { continue; }
                    let k = ((e[a] - pts[0][a]) / per[a]).round();
                    if k != 0.0 { for q in pts.iter_mut() { q[a] += k * per[a]; } }
                }
            }
            if let Some((q, r)) = tail.take() { if pole(q, pts[0]) { uv.push(q); refs.push(r); } }
            last = pts.last().copied();
            for (j, &k) in order.iter().enumerate().take(n - 1) { uv.push(pts[j]); refs.push((ce.edge, k, k.min(order[j + 1]))); }
            tail = Some((pts[n - 1], (ce.edge, order[n - 1], POLE)));
        }
        if let Some((q, r)) = tail { if !uv.is_empty() && pole(q, uv[0]) { uv.push(q); refs.push(r); } }
        rings.push((uv, refs, dir));
    }
    // a vertex met again in the face (a hole touching another loop) takes its first UV when the two agree within
    // rounding (1e-9 of the face's UV extent); seam copies a period apart stay
    let mut ext = 0.0f64;
    for (uv, _, _) in &rings { for q in uv { ext = ext.max(q[0].abs()).max(q[1].abs()); } }
    let mut first: HashMap<usize, P2> = HashMap::new();
    for (uv, refs, _) in rings.iter_mut() {
        for (q, &(e, k, _)) in uv.iter_mut().zip(refs.iter()) {
            let n = fr[&e].len();
            let v = if k == 0 { m.edges[e][0] } else if k + 1 == n { m.edges[e][1] } else { continue };
            match first.get(&v) {
                Some(f) if (f[0] - q[0]).hypot(f[1] - q[1]) <= 1e-9 * ext.max(1.0) => *q = *f,
                Some(_) => {}
                None => { first.insert(v, *q); }
            }
        }
    }
    Ok(rings)
}

/// sample interval of a ring segment that lies on no edge (a pole segment between two parameter points of one vertex)
const POLE: usize = usize::MAX;

/// ring segments that properly cross another ring segment (sampling too coarse where loops come close, e.g. a
/// hole tangent to the outer loop), as (edge, sample interval) to bisect
fn crossing_segments(rings: &[(Vec<P2>, Vec<(usize, usize, usize)>, i32)]) -> Vec<(usize, usize)> {
    // segments (a, b, (edge, interval))
    let mut seg: Vec<(P2, P2, (usize, usize))> = vec![];
    for (uv, refs, _) in rings { for k in 0..uv.len() { seg.push((uv[k], uv[(k + 1) % uv.len()], (refs[k].0, refs[k].2))); } }
    let on_edge = |i: usize| seg[i].2 .1 != POLE;
    let mut order: Vec<usize> = (0..seg.len()).collect();
    let lo = |i: usize| seg[i].0[0].min(seg[i].1[0]);
    let hi = |i: usize| seg[i].0[0].max(seg[i].1[0]);
    order.sort_by(|&a, &b| lo(a).total_cmp(&lo(b)));
    let o2 = |a: P2, b: P2, c: P2| orient2d(a, b, c).0;
    let mut out = vec![];
    for (x, &i) in order.iter().enumerate() {
        for &j in &order[x + 1..] {
            if lo(j) > hi(i) { break; }
            let ((a, b, ri), (c, d, rj)) = (seg[i], seg[j]);
            if a == c || a == d || b == c || b == d { continue; }
            // a pole segment is no edge interval: only its partner can be bisected
            if o2(a, b, c) * o2(a, b, d) < 0 && o2(c, d, a) * o2(c, d, b) < 0 {
                if on_edge(i) { out.push(ri); }
                if on_edge(j) { out.push(rj); }
            }
        }
    }
    out
}

/// +1 when the loop's coedges chain end to start with their stored senses, -1 when they chain with all
/// senses negated (CGM stores the senses up to the loop's direction), None otherwise
fn loop_direction(m: &Model, l: std::ops::Range<usize>) -> Option<i32> {
    let ends = |c: usize, dir: i32| { let ce = &m.coedges[c]; let e = m.edges[ce.edge]; if ce.sense * dir > 0 { (e[0], e[1]) } else { (e[1], e[0]) } };
    let cs: Vec<usize> = l.collect();
    [1, -1].into_iter().find(|&d| (0..cs.len()).all(|j| ends(cs[j], d).1 == ends(cs[(j + 1) % cs.len()], d).0))
}

/// fractions in [0, 1] of the range [t0, t1] where the chord deviation of `f` is below `sag`
fn sample(f: impl Fn(f64) -> V, t0: f64, t1: f64, sag: f64) -> Vec<f64> {
    fn rec(f: &dyn Fn(f64) -> V, a: (f64, V), b: (f64, V), sag: f64, depth: u32, out: &mut Vec<f64>) {
        let m = 0.5 * (a.0 + b.0);
        let pm = f(m);
        if depth < 14 && norm(sub(pm, mulv(addv(a.1, b.1), 0.5))) > sag {
            rec(f, a, (m, pm), sag, depth + 1, out);
            rec(f, (m, pm), b, sag, depth + 1, out);
        } else { out.push(b.0); }
    }
    let g = |x: f64| f(t0 + x * (t1 - t0));
    let mut out = vec![0.0];
    for k in 0..4 {
        let (a, b) = (k as f64 / 4.0, (k + 1) as f64 / 4.0);
        rec(&g, (a, g(a)), (b, g(b)), sag, 0, &mut out);
    }
    out
}
