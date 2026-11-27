//! The uses of one CGM edge: per coedge its face (id, surface kind), p-curve kind, parameter range, UV and 3D
//! ends and middle.
//! Usage: cgm_edge IR EDGE
use brep::cgm::*;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let m = Model::read(&std::fs::read(&a[1]).unwrap()).unwrap();
    let e: usize = a[2].parse().unwrap();
    println!("edge {e}: vertices {:?}", m.edges[e]);
    for f in &m.faces {
        for (li, l) in f.loops.iter().enumerate() {
            for c in l.clone().filter(|&c| m.coedges[c].edge == e) {
                let ce = &m.coedges[c];
                let pc = &m.pcurves[ce.pcurve].1;
                let mid = (ce.t[0] + ce.t[1]) / 2.0;
                println!("  face {} (surface kind {}) loop {li} of {}, coedge {c} sense {} p-curve kind {} t {:?}", f.id, m.surface_kinds[f.surface], f.loops.len(), ce.sense, m.pcurve_kinds[ce.pcurve], ce.t);
                println!("    uv {:?} {:?} {:?}", pc.uv(ce.t[0]), pc.uv(mid), pc.uv(ce.t[1]));
                println!("    3D {:?} {:?} {:?}", m.point(ce, ce.t[0]), m.point(ce, mid), m.point(ce, ce.t[1]));
            }
        }
    }
}
