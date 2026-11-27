//! One CGM face's loops in UV: per loop the direction, coedges (p-curve kind, edge, sense, UV start -> end) and
//! the UV gap where the loop does not close.
//! Usage: cgm_face IR FACE_ID
use brep::cgm::*;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let m = Model::read(&std::fs::read(&a[1]).unwrap()).unwrap();
    let id: u32 = a[2].parse().unwrap();
    let f = m.faces.iter().find(|f| f.id == id).unwrap();
    println!("face {id}: surface kind {} sense {} dom {:?} periods {:?}, {} loops", m.surface_kinds[f.surface], f.sense, m.surfaces[f.surface].dom(), m.surfaces[f.surface].periods(), f.loops.len());
    for (li, l) in f.loops.iter().enumerate() {
        let mut ends = vec![];
        for c in l.clone() {
            let ce = &m.coedges[c];
            let pc = &m.pcurves[ce.pcurve].1;
            let (s, e) = (pc.uv(ce.t[0]), pc.uv(ce.t[1]));
            println!("  loop {li} coedge {c}: p-curve kind {} edge {} {:?} sense {} uv ({:.5}, {:.5}) -> ({:.5}, {:.5})", m.pcurve_kinds[ce.pcurve], ce.edge, m.edges[ce.edge], ce.sense, s[0], s[1], e[0], e[1]);
            ends.push(if ce.sense > 0 { (s, e) } else { (e, s) });
        }
        let gaps: Vec<f64> = (0..ends.len()).map(|k| { let (a, b) = (ends[k].1, ends[(k + 1) % ends.len()].0); (a[0] - b[0]).hypot(a[1] - b[1]) }).collect();
        println!("  loop {li} UV gaps between consecutive coedges (stored senses): {:?}", gaps.iter().map(|g| format!("{g:.1e}")).collect::<Vec<_>>());
    }
}
