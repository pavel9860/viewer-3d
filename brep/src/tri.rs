//! Triangulation of a polygon with holes (de Berg et al., Computational Geometry, ch. 3): a sweep adds diagonals
//! that split the region into y-monotone pieces (holes need no bridges), each piece is triangulated with a stack.
//! Orientation decisions use exact predicates; the vertex order is total (y down, then x, then index), so
//! horizontal edges and repeated points are handled by the same rules.
use crate::exact::orient2d;
use crate::trim::P2;
use std::cmp::Ordering;

fn o2(a: P2, b: P2, c: P2) -> i8 { orient2d(a, b, c).0 }

/// rings: counter-clockwise outer boundaries and clockwise holes (region on the left of every edge).
/// Returns the points (ring points, consecutive repeats removed) and counter-clockwise triangles.
pub fn triangulate(rings: &[Vec<P2>]) -> (Vec<P2>, Vec<[usize; 3]>) {
    // vertices with ring neighbours
    let (mut p, mut prev, mut next): (Vec<P2>, Vec<usize>, Vec<usize>) = (vec![], vec![], vec![]);
    for r in rings {
        let mut q: Vec<P2> = vec![];
        for &x in r { if q.last() != Some(&x) { q.push(x); } }
        while q.len() > 1 && q[0] == q[q.len() - 1] { q.pop(); }
        if q.len() < 3 { continue; }
        let s = p.len();
        let n = q.len();
        for (k, x) in q.into_iter().enumerate() { p.push(x); prev.push(s + (k + n - 1) % n); next.push(s + (k + 1) % n); }
    }
    let n = p.len();
    // a above b: greater y, then smaller x, then smaller index
    let above = |a: usize, b: usize| -> bool {
        match p[b][1].total_cmp(&p[a][1]) { Ordering::Less => true, Ordering::Greater => false, Ordering::Equal => (p[a][0], a) < (p[b][0], b) }
    };
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| if a == b { Ordering::Equal } else if above(a, b) { Ordering::Less } else { Ordering::Greater });

    // sweep: status = edges i -> next[i] crossing the sweep line, with their helper vertex
    let mut status: Vec<(usize, usize)> = vec![];
    let mut merge = vec![false; n];
    let mut diags: Vec<(usize, usize)> = vec![];
    // edge e as (upper, lower) vertex
    let ends = |e: usize| if above(e, next[e]) { (e, next[e]) } else { (next[e], e) };
    // q east of edge e (on the left of the edge walked downward)
    let east = |e: usize, q: P2| { let (u, l) = ends(e); o2(p[u], p[l], q) > 0 };
    // edge f east of edge e on the sweep line (non-crossing edges): judged at an endpoint of one inside the other's span
    let within = |x: usize, e: usize| { let (u, l) = ends(e); !above(x, u) && !above(l, x) };
    let east_of = |f: usize, e: usize| -> bool {
        let (fu, fl) = ends(f);
        if within(fu, e) && p[fu] != p[ends(e).0] && p[fu] != p[ends(e).1] { return east(e, p[fu]); }
        if within(fl, e) && p[fl] != p[ends(e).0] && p[fl] != p[ends(e).1] { return east(e, p[fl]); }
        let (eu, el) = ends(e);
        if within(eu, f) && p[eu] != p[fu] && p[eu] != p[fl] { return !east(f, p[eu]); }
        if within(el, f) { return !east(f, p[el]); }
        east(e, p[fu])
    };
    // the status edge directly left (west) of v: the most eastern of the edges that v is east of
    let left_of = |status: &Vec<(usize, usize)>, v: usize| -> Option<usize> {
        let mut best: Option<usize> = None;
        for (k, &(e, _)) in status.iter().enumerate() {
            if !east(e, p[v]) { continue; }
            if best.map_or(true, |b| east_of(e, status[b].0)) { best = Some(k); }
        }
        best
    };
    let remove = |status: &mut Vec<(usize, usize)>, e: usize| if let Some(k) = status.iter().position(|s| s.0 == e) { status.swap_remove(k); };
    let helper = |status: &Vec<(usize, usize)>, e: usize| status.iter().find(|s| s.0 == e).map(|s| s.1);
    for &v in &order {
        let (pv, nv) = (prev[v], next[v]);
        let (pa, na) = (above(pv, v), above(nv, v));
        let convex = o2(p[pv], p[v], p[nv]) > 0;
        if !pa && !na {
            // start or split
            if !convex {
                if let Some(k) = left_of(&status, v) { diags.push((v, status[k].1)); status[k].1 = v; }
            }
            status.push((v, v));
        } else if pa && na {
            // end or merge
            if let Some(h) = helper(&status, pv) { if merge[h] { diags.push((v, h)); } }
            remove(&mut status, pv);
            if !convex {
                merge[v] = true;
                if let Some(k) = left_of(&status, v) { if merge[status[k].1] { diags.push((v, status[k].1)); } status[k].1 = v; }
            }
        } else if pa {
            // regular vertex with the region to its right (boundary runs downward)
            if let Some(h) = helper(&status, pv) { if merge[h] { diags.push((v, h)); } }
            remove(&mut status, pv);
            status.push((v, v));
        } else if let Some(k) = left_of(&status, v) {
            // regular vertex with the region to its left
            if merge[status[k].1] { diags.push((v, status[k].1)); }
            status[k].1 = v;
        }
    }

    // monotone pieces: faces of the ring edges plus both directions of every diagonal (region on the left)
    let mut he: Vec<(usize, usize)> = (0..n).map(|v| (v, next[v])).collect();
    for &(a, b) in &diags { if a != b { he.push((a, b)); he.push((b, a)); } }
    let mut out: Vec<Vec<usize>> = vec![vec![]; n];
    for (h, &(a, _)) in he.iter().enumerate() { out[a].push(h); }
    // next half-edge: at the end vertex b, the outgoing edge with the smallest clockwise turn from the reversed
    // edge b -> a (exact: groups by side of b -> a, then orientation within a group)
    let succ = |h: usize| -> usize {
        let (a, b) = he[h];
        let (pa, pb) = (p[a], p[b]);
        let group = |g: usize| { let q = p[he[g].1]; match o2(pb, pa, q) { -1 => 0, 1 => 2, _ => if (q[0] - pb[0]) * (pa[0] - pb[0]) + (q[1] - pb[1]) * (pa[1] - pb[1]) < 0.0 { 1 } else { 3 } } };
        *out[b].iter().min_by(|&&g, &&k| group(g).cmp(&group(k)).then_with(|| o2(pb, p[he[g].1], p[he[k].1]).cmp(&0))).unwrap()
    };
    let mut seen = vec![false; he.len()];
    let mut tris = vec![];
    for h0 in 0..he.len() {
        if seen[h0] { continue; }
        let mut face = vec![];
        let mut h = h0;
        while !seen[h] { seen[h] = true; face.push(he[h].0); h = succ(h); }
            monotone(&p, &face, &above, &mut tris);
        }
    (p, tris)
}

/// stack triangulation of a y-monotone counter-clockwise polygon
fn monotone(p: &[P2], face: &[usize], above: &dyn Fn(usize, usize) -> bool, tris: &mut Vec<[usize; 3]>) {
    let m = face.len();
    if m < 3 { return; }
    let mut emit = |a: usize, b: usize, c: usize| match o2(p[a], p[b], p[c]) { 1 => tris.push([a, b, c]), -1 => tris.push([a, c, b]), _ => {} };
    if m == 3 { emit(face[0], face[1], face[2]); return; }
    // left chain: from the top vertex along the boundary (counter-clockwise runs down the left side)
    let top = (0..m).min_by(|&a, &b| if above(face[a], face[b]) { Ordering::Less } else { Ordering::Greater }).unwrap();
    let bot = (0..m).min_by(|&a, &b| if above(face[a], face[b]) { Ordering::Greater } else { Ordering::Less }).unwrap();
    let mut left = vec![false; m];
    let mut k = top;
    while k != bot { left[k] = true; k = (k + 1) % m; }
    let mut u: Vec<usize> = (0..m).collect();
    u.sort_by(|&a, &b| if a == b { Ordering::Equal } else if above(face[a], face[b]) { Ordering::Less } else { Ordering::Greater });
    let mut st = vec![u[0], u[1]];
    for j in 2..m - 1 {
        let v = u[j];
        if left[v] != left[*st.last().unwrap()] {
            while st.len() > 1 { let a = st.pop().unwrap(); emit(face[v], face[a], face[*st.last().unwrap()]); }
            st.clear();
            st.push(u[j - 1]);
            st.push(v);
        } else {
            let mut last = st.pop().unwrap();
            while let Some(&t) = st.last() {
                // the diagonal v - t is inside when `last` is outside the triangle's side toward the chain
                let s = o2(p[face[t]], p[face[v]], p[face[last]]);
                if (left[v] && s < 0) || (!left[v] && s > 0) { emit(face[v], face[last], face[t]); last = st.pop().unwrap(); } else { break; }
            }
            st.push(last);
            st.push(v);
        }
    }
    let v = u[m - 1];
    while st.len() > 1 { let a = st.pop().unwrap(); emit(face[v], face[a], face[*st.last().unwrap()]); }
}
