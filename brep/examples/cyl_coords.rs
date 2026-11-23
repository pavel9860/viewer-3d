//! Faces' edge end points in the coordinates of the `brep bool` test cylinder (radius offset, angle, z).
//! Usage: cyl_coords FILE.stp STEP_FACE_ID...
use brep::geom::*;
use brep::step;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let im = step::import(&a[1]).unwrap();
    let sh = im.body.shells.iter().map(|x| x.1).filter(|s| s.closed).max_by_key(|s| s.faces.len()).unwrap();
    let b = im.body.extract(&sh.faces);
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for (_, v) in b.vertices.iter() { for k in 0..3 { lo[k] = lo[k].min(v.p[k]); hi[k] = hi[k].max(v.p[k]); } }
    let c = mulv(addv(lo, hi), 0.5);
    let r = 0.2 * (hi[0] - lo[0]).min(hi[1] - lo[1]).max(1e-3);
    println!("cylinder centre {c:?} r {r}");
    let cc = |p: V| { let q = sub(p, c); format!("(dr {:+.3e}, ang {:.5}, z {:+.5})", (q[0] * q[0] + q[1] * q[1]).sqrt() - r, q[1].atan2(q[0]).rem_euclid(std::f64::consts::TAU), q[2]) };
    for id in &a[2..] {
        let id: u32 = id.parse().unwrap();
        let (_, f) = b.faces.iter().find(|(_, f)| f.origin.src == id).unwrap();
        println!("face #{id}");
        for &l in &f.loops { for &ce in &b.loops.get(l).coedges { let e = b.edges.get(b.coedges.get(ce).edge); println!("  edge #{} {} -> {}", e.origin.src, cc(b.vertices.get(e.v[0]).p), cc(b.vertices.get(e.v[1]).p)); } }
    }
}
