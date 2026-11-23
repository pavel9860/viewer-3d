//! Primitive solids as exact B-rep bodies (one closed shell each).
use crate::geom::*;
use crate::nurbs::{ConicC, LineC};
use crate::topo::*;
use std::f64::consts::PI;
use std::sync::Arc;

const TOL: f64 = 1e-9;
const MADE: Origin = Origin { op: 0, src: 0 };

fn line(b: &mut Body, v0: Id<Vertex>, v1: Id<Vertex>) -> Id<Edge> {
    let (p0, p1) = (b.vertices.get(v0).p, b.vertices.get(v1).p);
    b.edges.push(Edge { curve: Some(Arc::new(LineC { o: p0, d: sub(p1, p0) })), v: [v0, v1], t: [0.0, 1.0], tol: TOL, origin: MADE })
}
fn face(b: &mut Body, s: Arc<dyn Surf>, loops: Vec<Vec<(Id<Edge>, bool)>>) -> Id<Face> {
    let loops = loops.into_iter().map(|l| { let coedges = l.into_iter().map(|(edge, sense)| b.coedges.push(Coedge { edge, sense })).collect(); b.loops.push(Loop { coedges }) }).collect();
    b.faces.push(Face { surface: s, sense: true, loops, tol: TOL, origin: MADE })
}

/// axis-aligned box: centre c, half sizes h
pub fn block(c: V, h: V) -> Body {
    let mut b = Body { tol: TOL, ..Default::default() };
    // vertex i: bit k set = +h[k]
    let vs: Vec<Id<Vertex>> = (0..8).map(|i| b.vertices.push(Vertex { p: [0, 1, 2].map(|k| c[k] + if i >> k & 1 == 1 { h[k] } else { -h[k] }), tol: TOL })).collect();
    let mut edges = std::collections::HashMap::new();
    for i in 0..8usize { for k in 0..3 { if i >> k & 1 == 0 { edges.insert((i, i | 1 << k), line(&mut b, vs[i], vs[i | 1 << k])); } } }
    let use_ = |a: usize, z: usize| if a < z { (edges[&(a, z)], true) } else { (edges[&(z, a)], false) };
    let mut faces = vec![];
    for k in 0..3 {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        for side in [0usize, 1] {
            // corners of the face in the (i, j) plane, counter-clockwise seen from outside
            let base = side << k;
            let mut q = [base, base | 1 << i, base | 1 << i | 1 << j, base | 1 << j];
            let (mut a, mut bb) = ([0.0; 3], [0.0; 3]);
            a[i] = 1.0; bb[j] = 1.0;
            if side == 0 { q.reverse(); std::mem::swap(&mut a, &mut bb); }
            let mut o = c;
            o[k] += if side == 1 { h[k] } else { -h[k] };
            let ext = 1.6 * h[i].max(h[j]);
            let lp = (0..4).map(|m| use_(q[m], q[(m + 1) % 4])).collect();
            faces.push(face(&mut b, Arc::new(Plane { o, a, b: bb, ext }), vec![lp]));
        }
    }
    b.shells.push(Shell { faces, closed: true });
    b
}

/// capped cylinder: centre c, radius r, half height h, axes (x, y, z); the seam runs at angle 0 (along x)
pub fn cylinder(c: V, r: f64, h: f64, x: V, y: V, z: V) -> Body {
    let mut b = Body { tol: TOL, ..Default::default() };
    let at = |zz: f64| addv(c, addv(mulv(x, r), mulv(z, zz)));
    let (v0, v1) = (b.vertices.push(Vertex { p: at(-h), tol: TOL }), b.vertices.push(Vertex { p: at(h), tol: TOL }));
    let circle = |b: &mut Body, v: Id<Vertex>, zz: f64| b.edges.push(Edge { curve: Some(Arc::new(ConicC { c: addv(c, mulv(z, zz)), x, y, a: r, b: r })), v: [v, v], t: [0.0, 2.0 * PI], tol: TOL, origin: MADE });
    let (bot, top) = (circle(&mut b, v0, -h), circle(&mut b, v1, h));
    let seam = line(&mut b, v0, v1);
    let lat: Arc<dyn Surf> = Arc::new(Moved { s: Box::new(Cyl { r, h, x, y, z }), r: frame([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]), t: c });
    let faces = vec![
        face(&mut b, lat, vec![vec![(bot, true), (seam, true), (top, false), (seam, false)]]),
        face(&mut b, Arc::new(Plane { o: addv(c, mulv(z, h)), a: x, b: y, ext: 1.6 * r }), vec![vec![(top, true)]]),
        face(&mut b, Arc::new(Plane { o: addv(c, mulv(z, -h)), a: y, b: x, ext: 1.6 * r }), vec![vec![(bot, false)]]),
    ];
    b.shells.push(Shell { faces, closed: true });
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn primitives_pass_the_checker() {
        for b in [block([0.0; 3], [1.0, 2.0, 3.0]), cylinder([0.0; 3], 0.5, 2.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])] {
            let r = crate::check::check(&b, 5, 1);
            assert_eq!(r.violations(), 0, "{r:?}");
        }
    }
}
