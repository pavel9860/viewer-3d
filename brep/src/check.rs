//! Body checker: topology invariants and geometry within tolerance.
//! Topology: loops close; in a closed shell every edge is used exactly twice with opposite senses.
//! Geometry: vertices lie on their edge curves, edge curves lie on the surfaces of the faces that use them.
use crate::geom::{dist, norm, sub};
use crate::topo::*;
use std::collections::HashMap;

#[derive(Default, Debug)]
pub struct Report {
    pub loops_open: usize,
    /// closed-shell edges not used exactly twice
    pub edge_use: usize,
    /// closed-shell edges used twice with the same sense (inconsistent orientation)
    pub orientation: usize,
    pub vertex_off_curve: usize,
    pub edge_off_surface: usize,
    pub edges_without_curve: usize,
    pub max_vertex_dev: f64,
    pub max_edge_dev: f64,
}
impl Report {
    pub fn violations(&self) -> usize { self.loops_open + self.edge_use + self.orientation + self.vertex_off_curve + self.edge_off_surface }
}

/// check with `samples` interior points per edge, using `workers` threads for the geometric part
pub fn check(b: &Body, samples: usize, workers: usize) -> Report {
    let mut r = Report::default();
    // loops close
    for (_, l) in b.loops.iter() {
        let n = l.coedges.len();
        r.loops_open += (0..n).filter(|&i| b.coedge_ends(l.coedges[i]).1 != b.coedge_ends(l.coedges[(i + 1) % n]).0).count();
    }
    // edge use in closed shells
    let mut uses: HashMap<Id<Edge>, Vec<bool>> = HashMap::new();
    let mut edge_faces: HashMap<Id<Edge>, Vec<Id<Face>>> = HashMap::new();
    for (_, sh) in b.shells.iter() {
        for &f in &sh.faces {
            for &l in &b.faces.get(f).loops {
                for &c in &b.loops.get(l).coedges {
                    let ce = b.coedges.get(c);
                    if sh.closed { uses.entry(ce.edge).or_default().push(ce.sense); }
                    edge_faces.entry(ce.edge).or_default().push(f);
                }
            }
        }
    }
    for s in uses.values() {
        if s.len() != 2 { r.edge_use += 1; } else if s[0] == s[1] { r.orientation += 1; }
    }
    // geometry, in parallel over edges
    let edges: Vec<(Id<Edge>, Vec<Id<Face>>)> = edge_faces.into_iter().collect();
    let debug = std::env::var("BREP_DEBUG").is_ok();
    let next = std::sync::atomic::AtomicUsize::new(0);
    let parts: Vec<(usize, usize, usize, f64, f64)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..workers.max(1)).map(|_| sc.spawn(|| {
            let mut acc = (0usize, 0usize, 0usize, 0.0f64, 0.0f64);
            loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= edges.len() { break; }
                let (eid, fs) = &edges[i];
                let e = b.edges.get(*eid);
                let Some(c) = &e.curve else { acc.2 += 1; continue };
                for k in 0..2 {
                    let d = norm(sub(c.ev(e.t[k]).0, b.vertices.get(e.v[k]).p));
                    acc.4 = acc.4.max(d);
                    if d > e.tol.max(b.vertices.get(e.v[k]).tol) { acc.0 += 1; if debug { eprintln!("vertex off curve: edge STEP #{} end {k} dev {d:.1e}", e.origin.src); } }
                }
                let mut fs = fs.clone();
                fs.sort_unstable_by_key(|f| f.0); fs.dedup();
                let dev = (1..=samples).map(|k| c.ev(e.t[0] + (e.t[1] - e.t[0]) * k as f64 / (samples + 1) as f64).0)
                    .flat_map(|p| fs.iter().map(move |f| (p, *f)))
                    .map(|(p, f)| dist(b.faces.get(f).surface.as_ref(), p)).fold(0.0, f64::max);
                acc.3 = acc.3.max(dev);
                if dev > e.tol { acc.1 += 1; if debug { let p = c.ev(0.5 * (e.t[0] + e.t[1])).0; for f in &fs { let s = b.faces.get(*f).surface.as_ref(); let uv = s.uv_of(p); eprintln!("edge off surface: edge STEP #{} dev {dev:.1e} face #{} dom {:?} uv {:?} d {:.1e} t {:?}", e.origin.src, b.faces.get(*f).origin.src, s.dom(), uv, dist(s, p), e.t); } } }
            }
            acc
        })).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for p in parts {
        r.vertex_off_curve += p.0;
        r.edge_off_surface += p.1;
        r.edges_without_curve += p.2;
        r.max_edge_dev = r.max_edge_dev.max(p.3);
        r.max_vertex_dev = r.max_vertex_dev.max(p.4);
    }
    r
}
