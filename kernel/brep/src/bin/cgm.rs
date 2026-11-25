//! cgm IR... [--sag S] [--sheets] [--all-solids] [--bodies] [--stl OUT.stl]
//!
//! Meshes CATIA CGM B-reps (IR from `cadcli cgm PART OUT`) and prints one JSON line per part:
//! faces, left-out faces, triangles, area, volume, closure. By default the shells the part shows
//! (`cgm_mesh::shown_shells`); `--all-solids` meshes every solid; `--bodies` prints one line per body
//! with the volume of its solids and its flags instead.
use brep::cgm::Model;
use brep::cgm_mesh::{shown_shells, stats, tessellate};
use std::io::Write;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut files = vec![];
    let (mut sag_rel, mut stl, mut sheets, mut all, mut bodies) = (1e-3, None, false, false, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--sag" => { sag_rel = args[i + 1].parse().expect("--sag S (fraction of the part's diagonal)"); i += 1; }
            "--stl" => { stl = Some(args[i + 1].clone()); i += 1; }
            "--sheets" => sheets = true,
            "--all-solids" => all = true,
            "--bodies" => bodies = true,
            f => files.push(f.to_string()),
        }
        i += 1;
    }
    if files.is_empty() { eprintln!("usage: cgm IR... [--sag S] [--sheets] [--all-solids] [--bodies] [--stl OUT.stl]"); std::process::exit(2); }
    let mut bad = 0;
    for f in &files {
        let t0 = Instant::now();
        let m = match std::fs::read(f).map_err(|e| e.to_string()).and_then(|b| Model::read(&b)) {
            Ok(m) => m,
            Err(e) => { println!("{{\"ir\":\"{f}\",\"error\":\"{e}\"}}"); bad += 1; continue; }
        };
        let solids: Vec<usize> = (0..m.shells.len()).filter(|&i| m.shells[i].solid).collect();
        // sag relative to the size of the solids' edge end points
        let b = stats(&tessellate(&m, f64::MAX, &solids)).bbox;
        let diag = ((b[3] - b[0]).powi(2) + (b[4] - b[1]).powi(2) + (b[5] - b[2]).powi(2)).sqrt();
        let sag = sag_rel * if diag.is_finite() && diag > 0.0 { diag } else { 1.0 };
        if bodies {
            for bd in &m.bodies {
                let sh: Vec<usize> = bd.shells.iter().copied().filter(|&i| m.shells[i].solid).collect();
                let st = stats(&tessellate(&m, sag, &sh));
                println!("{{\"ir\":\"{f}\",\"body\":{},\"flags\":{},\"solid_shells\":{},\"volume\":{:.6},\"open_edges\":{}}}", bd.id, bd.flags, sh.len(), st.volume, st.open_edges);
            }
            continue;
        }
        let shells = if all { solids.clone() } else { shown_shells(&m, sheets) };
        let mesh = tessellate(&m, sag, &shells);
        let s = stats(&mesh);
        let bb = s.bbox;
        println!("{{\"ir\":\"{f}\",\"solids\":{},\"shells\":{},\"faces\":{},\"skipped\":{},\"failed\":{},\"tris\":{},\"area\":{:.6},\"volume\":{:.6},\"open_edges\":{},\"nonmanifold_edges\":{},\"bbox\":[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}],\"ms\":{:.1}}}",
            solids.len(), shells.len(), m.faces.len(), m.skipped.len(), mesh.failed.len(), s.tris, s.area, s.volume, s.open_edges, s.nonmanifold_edges, bb[0], bb[1], bb[2], bb[3], bb[4], bb[5], t0.elapsed().as_secs_f64() * 1e3);
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
