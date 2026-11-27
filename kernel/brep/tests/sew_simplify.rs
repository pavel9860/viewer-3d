//! Sewing restores a closed shell from loose faces; detail levels of a closed mesh stay closed and within their
//! tolerance.
use brep::geom::*;
use brep::tess::{open_edges, tessellate, Mesh, Opts};
use brep::topo::*;

/// every face of `b` as its own open shell with its own vertices (moved by up to `jitter`) and edges; faces with
/// odd index reversed (loops and sense)
fn explode(b: &Body, jitter: f64) -> Body {
    let mut o = Body { tol: b.tol, ..Default::default() };
    let mut k = 0u32;
    for (_, sh) in b.shells.iter() {
        for &f in &sh.faces {
            let fc = b.faces.get(f);
            let mut vmap = std::collections::HashMap::new();
            let mut loops = vec![];
            for &l in &fc.loops {
                let mut cs = vec![];
                for &c in &b.loops.get(l).coedges {
                    let ce = b.coedges.get(c);
                    let ed = b.edges.get(ce.edge);
                    let v = ed.v.map(|v| *vmap.entry(v.0).or_insert_with(|| {
                        k += 1;
                        let d = [(k as f64 * 0.37).sin(), (k as f64 * 0.73).cos(), (k as f64 * 1.11).sin()];
                        o.vertices.push(Vertex { p: addv(b.vertices.get(v).p, mulv(d, jitter)), tol: b.tol })
                    }));
                    let e = o.edges.push(Edge { curve: ed.curve.clone(), v, t: ed.t, tol: ed.tol, origin: ed.origin });
                    cs.push(o.coedges.push(Coedge { edge: e, sense: ce.sense }));
                }
                loops.push(cs);
            }
            let flip = o.faces.len() % 2 == 1;
            if flip { for cs in loops.iter_mut() { cs.reverse(); for &c in cs.iter() { let s = o.coedges.get(c).sense; o.coedges.get_mut(c).sense = !s; } } }
            let loops = loops.into_iter().map(|coedges| o.loops.push(Loop { coedges })).collect();
            let nf = o.faces.push(Face { surface: fc.surface.clone(), sense: fc.sense ^ flip, loops, tol: fc.tol, origin: fc.origin });
            o.shells.push(Shell { faces: vec![nf], closed: false });
        }
    }
    o
}

fn volume(m: &Mesh) -> f64 { m.tri.iter().map(|t| dot(m.pos[t[0] as usize], cross(m.pos[t[1] as usize], m.pos[t[2] as usize])) / 6.0).sum() }

fn mesh_all(b: &Body) -> Mesh {
    let faces: Vec<Id<Face>> = b.shell_faces().collect();
    tessellate(b, &faces, &Opts { chord: 1e-3, ..Default::default() }, 2).0
}

#[test]
fn loose_faces_sew_into_a_closed_shell() {
    for (name, solid) in [("block", brep::make::block([0.1, 0.0, -0.2], [0.3, 0.2, 0.25])), ("cylinder", brep::make::cylinder([0.0; 3], 0.2, 0.5, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))] {
        let want = volume(&mesh_all(&solid));
        let mut b = explode(&solid, 0.3 * 1e-6);
        let st = brep::sew::sew(&mut b, 1e-6);
        assert_eq!((st.free_after, st.shells_after, st.closed_after, st.conflicts), (0, 1, 1, 0), "{name}: {st:?}");
        let m = mesh_all(&b);
        assert_eq!(open_edges(&m), (0, 0), "{name}: sewn mesh not closed");
        assert!((volume(&m) - want).abs() < 1e-6 * want.abs(), "{name}: volume {} vs {want}", volume(&m));
    }
}

#[test]
fn levels_stay_closed_within_tolerance() {
    // UV sphere of radius 0.5 (closed, manifold)
    let (nu, nv) = (48, 24);
    let mut pos = vec![[0.0, 0.0, 0.5]];
    for i in 1..nv {
        let th = std::f64::consts::PI * i as f64 / nv as f64;
        for j in 0..nu { let ph = 2.0 * std::f64::consts::PI * j as f64 / nu as f64; pos.push([0.5 * th.sin() * ph.cos(), 0.5 * th.sin() * ph.sin(), 0.5 * th.cos()]); }
    }
    pos.push([0.0, 0.0, -0.5]);
    let ring = |i: usize, j: usize| (1 + (i - 1) * nu + j % nu) as u32;
    let bot = (pos.len() - 1) as u32;
    let mut tri = vec![];
    for j in 0..nu { tri.push([0, ring(1, j), ring(1, j + 1)]); }
    for i in 1..nv - 1 { for j in 0..nu { tri.push([ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)]); tri.push([ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)]); } }
    for j in 0..nu { tri.push([ring(nv - 1, j), bot, ring(nv - 1, j + 1)]); }
    let tag = vec![0u32; tri.len()];
    let tols = [1e-3, 4e-3, 1.6e-2];
    let lv = brep::simplify::levels(&pos, &tri, &tag, &tols, 2);
    assert!(!lv.is_empty());
    let mut last = tri.len();
    for l in &lv {
        let m = Mesh { pos: l.pos.clone(), tri: l.tri.clone(), face: l.tag.clone() };
        assert_eq!(open_edges(&m), (0, 0), "level of {} triangles not closed", l.tri.len());
        assert!(l.tri.len() < last);
        last = l.tri.len();
        assert!(volume(&m) > 0.0);
    }
    // each level's measured deviation is within the largest tolerance it was made with
    for (l, &t) in lv.iter().zip(&tols) { assert!(l.err <= t, "deviation {} over tolerance {t}", l.err); }
}
