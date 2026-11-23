//! Debug one imported face: loops in UV (areas, extents, items). Usage: face_dbg FILE.stp STEP_FACE_ID
use brep::{mesh, step, trim::area};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let im = step::import(&a[1]).unwrap();
    let id: u32 = a[2].parse().unwrap();
    let (f, fc) = im.body.faces.iter().find(|(_, f)| f.origin.src == id).unwrap();
    let t = mesh::face_trim(&im.body, f);
    println!("face #{id} sense {} dom {:?} per {:?} sg {}", fc.sense, fc.surface.dom(), t.per, t.sg);
    for (i, (l, p)) in t.loops.iter().zip(t.polygons()).enumerate() {
        println!(" loop {i}: area {:.4e}, {} pts, items {}", area(p), p.len(), l.len());
        for it in l { println!("   edge #{} sense {} {} pts uv {:?} -> {:?}", im.body.edges.get(it.e).origin.src, it.sense, it.uv.len(), it.uv[0], it.uv[it.uv.len() - 1]); }
    }
}
