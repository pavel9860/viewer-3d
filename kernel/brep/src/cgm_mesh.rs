//! Triangle mesh of a CGM model (`cgm::Model`).
//!
//! Each edge is sampled once, along its first coedge's p-curve, until the chord deviation is below
//! `sag`; every face using the edge reuses those 3D points (watertight by construction) and maps the
//! samples to its own p-curve by the coedge's parameter range. A face's loops become UV polygons
//! (outer loop counter-clockwise, holes clockwise: CGM stores each loop's senses up to one direction
//! for the whole loop), are ear-clipped, and interior edges are bisected until their 3D chord
//! deviation is below `sag`.
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
}

/// open mesh edges by cause: on a B-rep edge used by one meshed coedge (its other face is not in the mesh), on a
/// B-rep edge used by two or more (the tessellation does not join), or not on a B-rep edge (inside a face)
#[derive(Debug, Default)]
pub struct OpenCauses { pub lone: usize, pub joined: usize, pub inside: usize, pub lone_edges: Vec<usize>, pub joined_edges: Vec<usize> }

pub fn open_causes(m: &Mesh) -> OpenCauses {
    let mut uses: HashMap<(u32, u32), i32> = HashMap::new();
    for t in &m.tris { for k in 0..3 { let (i, j) = (t[k], t[(k + 1) % 3]); *uses.entry((i.min(j), i.max(j))).or_default() += 1; } }
    let mut on: HashMap<(u32, u32), usize> = HashMap::new();
    for (e, pts) in &m.edge_pts { for w in pts.windows(2) { on.insert((w[0].min(w[1]), w[0].max(w[1])), *e); } }
    let mut c = OpenCauses::default();
    let mut keys: Vec<_> = uses.iter().filter(|x| *x.1 == 1).map(|x| *x.0).collect();
    keys.sort_unstable();
    for k in keys {
        match on.get(&k) {
            None => c.inside += 1,
            Some(&e) if m.edge_uses.get(&e).copied().unwrap_or(0) < 2 => { c.lone += 1; c.lone_edges.push(e); }
            Some(&e) => { c.joined += 1; c.joined_edges.push(e); }
        }
    }
    c.lone_edges.sort_unstable();
    c.lone_edges.dedup();
    c.joined_edges.sort_unstable();
    c.joined_edges.dedup();
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

/// meshes the faces of the given shells
pub fn tessellate(m: &Model, sag: f64, shells: &[usize]) -> Mesh {
    let mut out = Mesh { p: vec![], tris: vec![], tri_face: vec![], faces: vec![], failed: vec![], edge_pts: vec![], edge_uses: HashMap::new() };
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
        for &(fi, _) in &faces {
            if !valid(fi) || last.get(&fi) == Some(&0) { continue; }
            let Ok(rings) = face_rings(m, fi, &fr, &along) else { continue };
            let x = crossing_segments(&rings);
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
    for &(fi, shell_sense) in &faces {
        let face = &m.faces[fi];
        if let Some(why) = &m.bad_surface[face.surface] { out.failed.push((face.id, format!("surface kind {}: {why}", m.surface_kinds[face.surface]))); continue; }
        if let Some(c) = face.loops.iter().flat_map(|l| l.clone()).find(|&c| m.bad_pcurve[m.coedges[c].pcurve].is_some()) {
            let pc = m.coedges[c].pcurve;
            out.failed.push((face.id, format!("p-curve kind {}: {}", m.pcurve_kinds[pc], m.bad_pcurve[pc].as_ref().unwrap())));
            continue;
        }
        let s = m.surfaces[face.surface].as_ref();
        // UV rings with their mesh point indices
        let (mut rings, dirs) = match face_rings(m, fi, &fr, &along) {
            Ok(r) => {
                let dirs: Vec<i32> = r.iter().map(|x| x.2).collect();
                (r.into_iter().filter(|x| x.0.len() >= 3).map(|(uv, refs, _)| (uv, refs.iter().map(|&(e, k, _)| edges[&e].idx[k]).collect::<Vec<u32>>())).collect::<Vec<_>>(), dirs)
            }
            Err(why) => { out.failed.push((face.id, why)); continue; }
        };

        // boundary within the surface's parameter domain (surfaces with finite data: knots, profile ranges), with
        // a margin of 1e-3 of its width
        let d = s.dom();
        if let Some(q) = rings.iter().flat_map(|r| &r.0).find(|&&q| !inside(q, d)) {
            out.failed.push((face.id, format!("boundary point ({:.6}, {:.6}) outside the surface domain {d:?}", q[0], q[1])));
            continue;
        }
        if rings.is_empty() { out.failed.push((face.id, "no loop with 3 points".into())); continue; }
        // the face's own p-curves must put its boundary on the shared edge points (within the chord tolerance)
        let off = rings.iter().flat_map(|r| r.0.iter().zip(&r.1)).map(|(q, &i)| norm(sub(s.ev(q[0], q[1]).p, out.p[i as usize]))).fold(0.0, f64::max);
        if off > sag { out.failed.push((face.id, format!("boundary {off:.3e} off the shared edge points (p-curves disagree with the neighbouring faces)"))); continue; }
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
        let (pts, src, tris) = crate::tri::triangulate_src(&uv_rings);
        let mut index: Vec<u32> = src.iter().map(|&(r, k)| rings[r].1[k]).collect();
        let target = area(&rings[0].0).abs() - rings[1..].iter().map(|h| area(&h.0).abs()).sum::<f64>();
        let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]]).abs()).sum();
        if (got - target).abs() > 1e-6 * target.max(1e-300) { out.failed.push((face.id, format!("triangulation covered {got:.6e} of {target:.6e} (UV area)"))); }
        // every ring segment between distinct points must be a triangle edge
        let mut tedges: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
        for t in &tris { for k in 0..3 { let (x, y) = (index[t[k]], index[t[(k + 1) % 3]]); tedges.insert((x.min(y), x.max(y))); } }
        let lost = rings.iter().flat_map(|r| (0..r.1.len()).map(move |k| (r.1[k], r.1[(k + 1) % r.1.len()]))).filter(|&(x, y)| x != y && out.p[x as usize] != out.p[y as usize] && !tedges.contains(&(x.min(y), x.max(y)))).count();
        // BREP_FACE=ID: the face's UV rings and triangles
        if std::env::var("BREP_FACE").ok().and_then(|x| x.parse::<u32>().ok()) == Some(face.id) {
            for (ri, r) in rings.iter().enumerate() { for (k, q) in r.0.iter().enumerate() { eprintln!("ring {ri} {k} idx {} uv {:?} {:?}", r.1[k], q[0], q[1]); } }
            for t in &tris { eprintln!("tri {:?} idx {:?}", t, t.map(|i| index[i])); }
        }
        if lost > 0 { out.failed.push((face.id, format!("{lost} boundary segments are not triangle edges"))); }
        let mut uvs = pts;
        let mut tri: Vec<[usize; 3]> = tris;
        refine(s, &mut uvs, &mut tri, &mut index, &mut out.p, sag);
        // UV counter-clockwise gives the normal Su x Sv. A shell's sense +1 puts the material on the
        // surface-normal side, so the outward normal is -sense * (Su x Sv) (2,930 of 2,934 faces of
        // the closed meshes of the corpus; the face record's own i16 does not enter)
        let flip = shell_sense > 0;
        out.faces.push((face.id, face.sense, shell_sense, dirs[outer], chained));
        for t in tri {
            let (a, b, c) = (index[t[0]], index[t[1]], index[t[2]]);
            out.tris.push(if flip { [a, c, b] } else { [a, b, c] });
            out.tri_face.push((out.faces.len() - 1) as u32);
        }
    }
    out
}


/// UV rings of a face from edge sample fractions `fr`, each point with its (edge, sample index, interval of the
/// segment starting there), and the loop direction. On a periodic surface each coedge is shifted by whole periods to start where the previous one ended.
fn face_rings(m: &Model, fi: usize, fr: &HashMap<usize, Vec<f64>>, along: &HashMap<usize, usize>) -> Result<Vec<(Vec<P2>, Vec<(usize, usize, usize)>, i32)>, String> {
    let face = &m.faces[fi];
    let per = m.surfaces[face.surface].periods();
    let mut rings = vec![];
    for l in &face.loops {
        let dir = loop_direction(m, l.clone()).ok_or("loop does not chain in either direction")?;
        let (mut uv, mut refs) = (vec![], vec![]);
        let mut last: Option<P2> = None;
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
                let (mut t, h) = (t0, 1e-7 * (ce.t[1] - ce.t[0]).abs().max(1e-12));
                for _ in 0..8 {
                    let q = m.point(ce, t);
                    let d = mulv(sub(m.point(ce, t + h), m.point(ce, t - h)), 0.5 / h);
                    let dd = dot(d, d);
                    if dd <= 0.0 { break; }
                    let step = dot(sub(q, target), d) / dd;
                    t -= step;
                    if step.abs() <= 1e-12 * (ce.t[1] - ce.t[0]).abs() { break; }
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
            last = pts.last().copied();
            for (j, &k) in order.iter().enumerate().take(n - 1) { uv.push(pts[j]); refs.push((ce.edge, k, k.min(order[j + 1]))); }
        }
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

/// ring segments that properly cross another ring segment (sampling too coarse where loops come close, e.g. a
/// hole tangent to the outer loop), as (edge, sample interval) to bisect
fn crossing_segments(rings: &[(Vec<P2>, Vec<(usize, usize, usize)>, i32)]) -> Vec<(usize, usize)> {
    // segments (a, b, (edge, interval))
    let mut seg: Vec<(P2, P2, (usize, usize))> = vec![];
    for (uv, refs, _) in rings { for k in 0..uv.len() { seg.push((uv[k], uv[(k + 1) % uv.len()], (refs[k].0, refs[k].2))); } }
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
            if o2(a, b, c) * o2(a, b, d) < 0 && o2(c, d, a) * o2(c, d, b) < 0 { out.push(ri); out.push(rj); }
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

/// bisects interior triangle edges (never boundary ones) whose midpoint lies more than `sag` off the
/// surface chord; `index` maps local UV points to mesh points, new points are appended to `p`
fn refine(s: &dyn Surf, uv: &mut Vec<P2>, tris: &mut Vec<[usize; 3]>, index: &mut Vec<u32>, p: &mut Vec<V>, sag: f64) {
    let nb = uv.len(); // boundary points come first; an edge between two of them that is a ring edge stays
    let key = |a: usize, b: usize| (a.min(b), a.max(b));
    let mut boundary: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    {
        let mut count: HashMap<(usize, usize), u32> = HashMap::new();
        for t in tris.iter() { for k in 0..3 { *count.entry(key(t[k], t[(k + 1) % 3])).or_default() += 1; } }
        for (e, n) in count { if n == 1 && e.0 < nb && e.1 < nb { boundary.insert(e); } }
    }
    let pt = |i: usize, p: &Vec<V>, index: &Vec<u32>| p[index[i] as usize];
    let budget = 400_000usize;
    for _round in 0..12 {
        let mut split: HashMap<(usize, usize), usize> = HashMap::new();
        for t in tris.iter() {
            for k in 0..3 {
                let e = key(t[k], t[(k + 1) % 3]);
                if boundary.contains(&e) || split.contains_key(&e) { continue; }
                let mid = [(uv[e.0][0] + uv[e.1][0]) / 2.0, (uv[e.0][1] + uv[e.1][1]) / 2.0];
                let q = s.ev(mid[0], mid[1]).p;
                if norm(sub(q, mulv(addv(pt(e.0, p, index), pt(e.1, p, index)), 0.5))) > sag && uv.len() < budget {
                    uv.push(mid);
                    p.push(q);
                    index.push((p.len() - 1) as u32);
                    split.insert(e, uv.len() - 1);
                }
            }
        }
        if split.is_empty() { break; }
        let mut next = Vec::with_capacity(tris.len() * 2);
        for t in tris.iter() {
            let m: Vec<Option<usize>> = (0..3).map(|k| split.get(&key(t[k], t[(k + 1) % 3])).copied()).collect();
            subdivide(*t, [m[0], m[1], m[2]], &mut next);
        }
        *tris = next;
    }
}

/// a triangle (a, b, c) with optional midpoints on ab, bc, ca, split into 1..4 triangles
fn subdivide(t: [usize; 3], m: [Option<usize>; 3], out: &mut Vec<[usize; 3]>) {
    let [a, b, c] = t;
    match m {
        [None, None, None] => out.push(t),
        [Some(x), None, None] => { out.push([a, x, c]); out.push([x, b, c]); }
        [None, Some(y), None] => { out.push([a, b, y]); out.push([a, y, c]); }
        [None, None, Some(z)] => { out.push([a, b, z]); out.push([z, b, c]); }
        [Some(x), Some(y), None] => { out.push([a, x, c]); out.push([x, b, y]); out.push([x, y, c]); }
        [None, Some(y), Some(z)] => { out.push([a, b, z]); out.push([z, b, y]); out.push([z, y, c]); }
        [Some(x), None, Some(z)] => { out.push([a, x, z]); out.push([x, b, c]); out.push([x, c, z]); }
        [Some(x), Some(y), Some(z)] => { out.push([a, x, z]); out.push([x, b, y]); out.push([z, y, c]); out.push([x, y, z]); }
    }
}
