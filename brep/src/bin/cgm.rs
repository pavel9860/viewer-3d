//! cgm IR... [--sag S] [--sheets] [--stl OUT.stl]: meshes CATIA CGM B-reps (IR from `cadcli cgm PART OUT`) and
//! prints one JSON line per part: faces, left-out faces, triangles, area, volume, closure.
use brep::cgm::Model;
use brep::cgm_mesh::{stats, tessellate};
use std::io::Write;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut files = vec![];
    let (mut sag_rel, mut stl, mut sheets, mut body) = (1e-3, None, false, None::<u32>);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--sag" => { sag_rel = args[i + 1].parse().expect("--sag S (fraction of the part's diagonal)"); i += 1; }
            "--sheets" => sheets = true,
            "--body" => { body = Some(args[i + 1].parse().expect("--body ID")); i += 1; }
            "--bodies" => body = Some(u32::MAX),
            "--stl" => { stl = Some(args[i + 1].clone()); i += 1; }
            f => files.push(f.to_string()),
        }
        i += 1;
    }
    if files.is_empty() { eprintln!("usage: cgm IR... [--sag S] [--sheets] [--stl OUT.stl]"); std::process::exit(2); }
    let mut bad = 0;
    for f in &files {
        let t0 = Instant::now();
        let m = match std::fs::read(f).map_err(|e| e.to_string()).and_then(|b| Model::read(&b)) {
            Ok(mut m) => {
                #[allow(unused_mut)]
                if body == Some(u32::MAX) {
                    // one line per body: its solid volume
                    let mut ids: Vec<u32> = m.shells.iter().filter(|s| s.solid).map(|s| s.body).collect();
                    ids.sort(); ids.dedup();
                    for b in ids {
                        let mut sub = Model::read(&std::fs::read(f).unwrap()).unwrap();
                        sub.shells.retain(|s| s.body == b);
                        let st = stats(&tessellate(&sub, f64::MAX, false));
                        println!("{{\"ir\":\"{f}\",\"body\":{b},\"volume\":{:.3},\"open_edges\":{}}}", st.volume, st.open_edges);
                    }
                    continue;
                }
                if let Some(b) = body { m.shells.retain(|s| s.body == b); }
                m
            }
            Err(e) => { println!("{{\"ir\":\"{f}\",\"error\":\"{e}\"}}"); bad += 1; continue; }
        };
        // sag relative to the size of the edges' end points
        let coarse = tessellate(&m, f64::MAX, sheets);
        let b = stats(&coarse).bbox;
        let diag = ((b[3] - b[0]).powi(2) + (b[4] - b[1]).powi(2) + (b[5] - b[2]).powi(2)).sqrt().max(1e-9);
        let mesh = tessellate(&m, sag_rel * diag, sheets);
        let s = stats(&mesh);
        println!("{{\"ir\":\"{f}\",\"solids\":{},\"faces\":{},\"skipped\":{},\"failed\":{},\"tris\":{},\"area\":{:.6},\"volume\":{:.6},\"open_edges\":{},\"nonmanifold_edges\":{},\"bbox\":[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}],\"ms\":{:.1}}}",
            m.shells.iter().filter(|s| s.solid).count(), m.faces.len(), m.skipped.len(), mesh.failed.len(), s.tris, s.area, s.volume, s.open_edges, s.nonmanifold_edges, b[0], b[1], b[2], b[3], b[4], b[5], t0.elapsed().as_secs_f64() * 1e3);
        if std::env::var("CGM_ORIENT").is_ok() {
            match brep::cgm_mesh::orientation_fix(&mesh) {
                Some(fix) => for (f, x) in mesh.faces.iter().zip(fix) { eprintln!("orient {} {} {} {} {} {}", f.0, f.1, f.2, f.3, f.4, x as i32); },
                None => eprintln!("orient none"),
            }
        }
        for (id, why) in &mesh.failed { eprintln!("{f}: face {id}: {why}"); }
        if let Some(out) = &stl {
            let mut w = std::io::BufWriter::new(std::fs::File::create(out).expect("cannot create the STL"));
            w.write_all(&[0u8; 80]).unwrap();
            w.write_all(&(mesh.tris.len() as u32).to_le_bytes()).unwrap();
            for t in &mesh.tris {
                w.write_all(&[0u8; 12]).unwrap();
                for &k in t { for c in mesh.p[k as usize] { w.write_all(&(c as f32).to_le_bytes()).unwrap(); } }
                w.write_all(&[0u8; 2]).unwrap();
            }
        }
    }
    std::process::exit((bad > 0) as i32);
}
