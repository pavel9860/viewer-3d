//! Per listed face: the largest gap between consecutive coedges of its loops in 3D and in UV (loop direction as
//! the mesher chains it), and the largest distance from a coedge's midpoint to its partner coedge (the other use
//! of the same edge), one line per face: ID SURFACE_KIND GAP_3D GAP_UV PARTNER_OFF.
//! Usage: cgm_loops IR FACE_ID...
use brep::cgm::*;
use brep::geom::*;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let m = Model::read(&std::fs::read(&a[1]).unwrap()).unwrap();
    let mut uses: std::collections::HashMap<usize, Vec<usize>> = Default::default();
    for f in &m.faces { for l in &f.loops { for c in l.clone() { uses.entry(m.coedges[c].edge).or_default().push(c); } } }
    let poly = |c: usize| -> Vec<V> { let ce = &m.coedges[c]; (0..=256).map(|k| m.point(ce, ce.t[0] + (ce.t[1] - ce.t[0]) * k as f64 / 256.0)).collect() };
    let seg = |p: V, a: V, b: V| { let d = sub(b, a); let l = dot(d, d); let t = if l > 0.0 { (dot(sub(p, a), d) / l).clamp(0.0, 1.0) } else { 0.0 }; norm(sub(p, addv(a, mulv(d, t)))) };
    for id in a[2..].iter().map(|x| x.parse::<u32>().unwrap()) {
        let f = m.faces.iter().find(|f| f.id == id).unwrap();
        let (mut g3, mut guv) = (0.0f64, 0.0f64);
        for l in &f.loops {
            let ends = |dir: i32| -> Vec<((V, [f64; 2]), (V, [f64; 2]))> {
                l.clone().map(|c| {
                    let ce = &m.coedges[c];
                    let pc = &m.pcurves[ce.pcurve].1;
                    let at = |t: f64| (m.point(ce, t), pc.uv(t));
                    if ce.sense * dir > 0 { (at(ce.t[0]), at(ce.t[1])) } else { (at(ce.t[1]), at(ce.t[0])) }
                }).collect()
            };
            // the direction whose 3D chaining is tighter
            let gap = |e: &Vec<((V, [f64; 2]), (V, [f64; 2]))>| (0..e.len()).map(|k| norm(sub(e[k].1 .0, e[(k + 1) % e.len()].0 .0))).fold(0.0, f64::max);
            let (e1, e2) = (ends(1), ends(-1));
            let e = if gap(&e1) <= gap(&e2) { e1 } else { e2 };
            g3 = g3.max(gap(&e));
            // UV gaps after shifting each coedge by whole periods to start where the previous one ended
            let per = m.surfaces[f.surface].periods();
            let mut shift = [0.0f64; 2];
            let mut prev_end: Option<[f64; 2]> = None;
            let mut first: Option<[f64; 2]> = None;
            for k in 0..e.len() {
                let (mut st, mut en) = (e[k].0 .1, e[k].1 .1);
                for a in 0..2 {
                    if per[a] > 0.0 { if let Some(pe) = prev_end { shift[a] = ((pe[a] - st[a]) / per[a]).round() * per[a]; } }
                    st[a] += shift[a];
                    en[a] += shift[a];
                }
                if let Some(pe) = prev_end { guv = guv.max((pe[0] - st[0]).hypot(pe[1] - st[1])); }
                first.get_or_insert(st);
                prev_end = Some(en);
            }
            if let (Some(a), Some(b)) = (first, prev_end) { guv = guv.max((a[0] - b[0]).hypot(a[1] - b[1])); }
        }
        let mut off = 0.0f64;
        for c in f.loops.iter().flat_map(|l| l.clone()) {
            let ce = &m.coedges[c];
            let mid = m.point(ce, (ce.t[0] + ce.t[1]) / 2.0);
            for &o in uses[&ce.edge].iter().filter(|&&o| o != c) {
                let q = poly(o);
                off = off.max(q.windows(2).map(|w| seg(mid, w[0], w[1])).fold(f64::MAX, f64::min));
            }
        }
        println!("{id} {} {g3:.3e} {guv:.3e} {off:.3e}", m.surface_kinds[f.surface]);
    }
}
