//! Time one face-pair SSI from a STEP file against the `brep bool` cylinder, counting surface evaluations.
//! Usage: ssi_pair FILE.stp STEP_FACE_ID
use brep::geom::*;
use brep::{make, mesh, ssi, step};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
struct Count<'a>(&'a dyn Surf, AtomicUsize);
impl Surf for Count<'_> {
    fn ev(&self, u: f64, v: f64) -> D2 { self.1.fetch_add(1, Relaxed); self.0.ev(u, v) }
    fn dom(&self) -> [f64; 4] { self.0.dom() }
    fn per_u(&self) -> bool { self.0.per_u() }
    fn per_v(&self) -> bool { self.0.per_v() }
    fn poles_v(&self) -> bool { self.0.poles_v() }
    fn implicit(&self, p: V) -> Option<f64> { self.0.implicit(p) }
    fn uv_of(&self, p: V) -> [f64; 2] { self.0.uv_of(p) }
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let im = step::import(&a[1]).unwrap();
    let id: u32 = a[2].parse().unwrap();
    let sh = im.body.shells.iter().map(|x| x.1).filter(|s| s.closed).max_by_key(|s| s.faces.len()).unwrap();
    let body = im.body.extract(&sh.faces);
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for (_, v) in body.vertices.iter() { for k in 0..3 { lo[k] = lo[k].min(v.p[k]); hi[k] = hi[k].max(v.p[k]); } }
    let r = 0.2 * (hi[0] - lo[0]).min(hi[1] - lo[1]).max(1e-3);
    let cyl = make::cylinder(mulv(addv(lo, hi), 0.5), r, (hi[2] - lo[2]).max(r), [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let (f, fc) = body.faces.iter().find(|(_, f)| f.origin.src == id).unwrap();
    let t = mesh::face_trim(&body, f);
    println!("face #{id} dom {:?} trim bb {:?}", fc.surface.dom(), t.bb);
    let (ca, cb) = (Count(fc.surface.as_ref(), AtomicUsize::new(0)), Count(cyl.faces.iter().next().unwrap().1.surface.as_ref(), AtomicUsize::new(0)));
    let t0 = std::time::Instant::now();
    let (cv, _) = ssi::intersect(&ca, &cb);
    println!("{} curves, {:.0} ms, evals A {} B {}", cv.len(), t0.elapsed().as_secs_f64() * 1e3, ca.1.load(Relaxed), cb.1.load(Relaxed));
    for c in &cv { let m = c.pts[c.pts.len() / 2]; println!("  mid uv on A {:?} (window form computed below)", fc.surface.uv_of(m)); }
    for c in &cv { println!("  full: {:?} {} pts {:?} .. {:?} ends {:?}", c.kind, c.pts.len(), c.pts[0], c.pts[c.pts.len() - 1], c.ends); }
    let (u0, u1) = { let x = t.raw([t.bb[0], 0.0])[0]; let y = t.raw([t.bb[1], 0.0])[0]; (x.min(y), x.max(y)) };
    let (pu, pv) = (0.02 * (u1 - u0) + 1e-9, 0.02 * (t.bb[3] - t.bb[2]) + 1e-9);
    let w = Window::new(fc.surface.as_ref(), [u0 - pu, u1 + pu, t.bb[2] - pv, t.bb[3] + pv]);
    println!("window {:?} per {} {}", w.dom(), w.per_u(), w.per_v());
    for c in ssi::intersect(&w, cb.0).0 { println!("  win: {:?} {} pts {:?} .. {:?} ends {:?}", c.kind, c.pts.len(), c.pts[0], c.pts[c.pts.len() - 1], c.ends); }
}
