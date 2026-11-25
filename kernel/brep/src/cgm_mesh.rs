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

/// meshes the faces of the solids, and of the sheets too when `sheets`. When a body named by the
/// stream's roots has a solid, only root bodies are meshed (the modeler also keeps intermediate
/// results); otherwise every solid is. Which body the part shows is decided by its feature tree,
/// not read yet, so parts keeping several solids without a root among them mesh all of them.
pub fn tessellate(m: &Model, sag: f64, sheets: bool) -> Mesh {
    let any_root = m.shells.iter().any(|s| s.root && s.solid);
    let mut out = Mesh { p: vec![], tris: vec![], tri_face: vec![], faces: vec![], failed: vec![] };
    let mut vtx: Vec<Option<u32>> = vec![None; m.vertices];
    let mut edges: HashMap<usize, EdgeSamples> = HashMap::new();
    // faces of the shells, each once
    let mut seen = vec![false; m.faces.len()];
    let faces: Vec<(usize, i32)> = m.shells.iter().filter(|s| (s.solid || sheets) && (s.root || !any_root)).flat_map(|s| s.faces.iter().copied()).filter(|&(f, _)| !std::mem::replace(&mut seen[f], true)).collect();
    let coedges: Vec<usize> = faces.iter().flat_map(|&(f, _)| m.faces[f].loops.iter().flat_map(|l| l.clone())).collect();
    for c in coedges {
        let ce = &m.coedges[c];
        if edges.contains_key(&ce.edge) { continue; }
        let f = sample(|t| m.point(ce, t), ce.t[0], ce.t[1], sag);
        let ends = m.edges[ce.edge];
        let mut idx = Vec::with_capacity(f.len());
        for (k, &x) in f.iter().enumerate() {
            let at_end = if k == 0 { Some(ends[0]) } else if k + 1 == f.len() { Some(ends[1]) } else { None };
            let new = |out: &mut Mesh| { out.p.push(m.point(ce, ce.t[0] + x * (ce.t[1] - ce.t[0]))); (out.p.len() - 1) as u32 };
            idx.push(match at_end {
                Some(v) => *vtx[v].get_or_insert_with(|| new(&mut out)),
                None => new(&mut out),
            });
        }
        edges.insert(ce.edge, EdgeSamples { f, idx });
    }
    for &(fi, shell_sense) in &faces {
        let face = &m.faces[fi];
        let s = m.surfaces[face.surface].as_ref();
        // UV rings with their mesh point indices
        let mut rings: Vec<(Vec<P2>, Vec<u32>)> = vec![];
        let mut bad_loop = false;
        let mut dirs = vec![];
        for l in &face.loops {
            let Some(dir) = loop_direction(m, l.clone()) else { bad_loop = true; break };
            dirs.push(dir);
            let (mut uv, mut ix) = (vec![], vec![]);
            for c in l.clone() {
                let ce = &m.coedges[c];
                let es = &edges[&ce.edge];
                let pc = &m.pcurves[ce.pcurve].1;
                let n = es.f.len();
                let order: Vec<usize> = if ce.sense * dir > 0 { (0..n - 1).collect() } else { (1..n).rev().collect() };
                for k in order {
                    uv.push(pc.uv(ce.t[0] + es.f[k] * (ce.t[1] - ce.t[0])));
                    ix.push(es.idx[k]);
                }
            }
            if uv.len() >= 3 { rings.push((uv, ix)); }
        }
        if bad_loop { out.failed.push((face.id, "loop does not chain in either direction".into())); continue; }
        if rings.is_empty() { out.failed.push((face.id, "no loop with 3 points".into())); continue; }
        // outer loop: largest area, counter-clockwise; the others clockwise
        let outer = (0..rings.len()).max_by(|&a, &b| area(&rings[a].0).abs().total_cmp(&area(&rings[b].0).abs())).unwrap();
        let chained = if area(&rings[outer].0) > 0.0 { 1 } else { -1 };
        for (i, r) in rings.iter_mut().enumerate() {
            let a = area(&r.0);
            if (i == outer) != (a > 0.0) { r.0.reverse(); r.1.reverse(); }
        }
        let holes: Vec<Vec<P2>> = rings.iter().enumerate().filter(|&(i, _)| i != outer).map(|(_, r)| r.0.clone()).collect();
        let (pts, tris) = triangulate(&rings[outer].0, &holes);
        let mut index: Vec<u32> = rings[outer].1.clone();
        for (i, r) in rings.iter().enumerate() { if i != outer { index.extend_from_slice(&r.1); } }
        let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]]).abs()).sum();
        let target = area(&rings[outer].0).abs() - holes.iter().map(|h| area(h).abs()).sum::<f64>();
        if (got - target).abs() > 1e-6 * target.max(1e-300) { out.failed.push((face.id, format!("ear clipping covered {got:.6e} of {target:.6e} (UV area)"))); }
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

fn o2(a: P2, b: P2, c: P2) -> i8 { orient2d(a, b, c).0 }
fn crosses(a: P2, b: P2, c: P2, d: P2) -> bool { o2(a, b, c) * o2(a, b, d) < 0 && o2(c, d, a) * o2(c, d, b) < 0 }

/// ear clipping of a counter-clockwise polygon with clockwise holes (each bridged to a visible
/// vertex). Every input point ends up in a triangle: collinear points, which a shared edge's
/// neighbour face also uses, are cut off as zero-area ears once no proper ear is left.
fn triangulate(outer: &[P2], holes: &[Vec<P2>]) -> (Vec<P2>, Vec<[usize; 3]>) {
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
    let mut flat_ok = false;
    while v.len() > 3 {
        let n = v.len();
        i %= n;
        let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
        let turn = o2(pts[a], pts[b], pts[c]);
        let ear = (turn > 0 || (flat_ok && turn == 0)) && (0..n).all(|k| {
            let q = v[k];
            if q == a || q == b || q == c || pts[q] == pts[a] || pts[q] == pts[b] || pts[q] == pts[c] { return true; }
            if o2(pts[v[(k + n - 1) % n]], pts[q], pts[v[(k + 1) % n]]) > 0 { return true; }
            !(o2(pts[a], pts[b], pts[q]) >= 0 && o2(pts[b], pts[c], pts[q]) >= 0 && o2(pts[c], pts[a], pts[q]) >= 0)
        });
        if ear { tris.push([a, b, c]); v.remove(i); misses = 0; flat_ok = false; i = i.saturating_sub(1); }
        else {
            i += 1;
            misses += 1;
            if misses > n { if flat_ok { break; } flat_ok = true; misses = 0; }
        }
    }
    if v.len() == 3 { tris.push([v[0], v[1], v[2]]); }
    (pts, tris)
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
