//! brep CLI. Usage: brep check FILE.stp|LIST.tsv
use brep::{check, step};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 || args[1] != "check" { eprintln!("usage: brep check FILE.stp|LIST.tsv"); std::process::exit(2); }
    let files: Vec<String> = if args[2].ends_with(".tsv") { std::fs::read_to_string(&args[2]).unwrap().lines().filter(|l| !l.is_empty()).map(|l| l.rsplit('\t').next().unwrap().to_string()).collect() } else { vec![args[2].clone()] };
    let mut bad = 0;
    for f in &files {
        let t0 = Instant::now();
        let im = match step::import(f) { Ok(i) => i, Err(e) => { println!("{f}: {e}"); bad += 1; continue } };
        let ti = t0.elapsed().as_secs_f64() * 1e3;
        let t1 = Instant::now();
        let r = check::check(&im.body, 7, 4);
        let tc = t1.elapsed().as_secs_f64() * 1e3;
        let mm = |x: f64| x * im.diag;
        println!("{:<26} {} | import {ti:.0} ms, check {tc:.0} ms | violations {} (open loops {}, edge use {}, orientation {}, vertex off curve {}, edge off surface {}) | max dev vertex {:.1e} mm edge {:.1e} mm | no curve {} | unsupported {:?}",
            f.rsplit('/').next().unwrap(), im.body.stats(), r.violations(), r.loops_open, r.edge_use, r.orientation, r.vertex_off_curve, r.edge_off_surface, mm(r.max_vertex_dev), mm(r.max_edge_dev), r.edges_without_curve, im.unsupported);
        bad += r.violations();
    }
    std::process::exit((bad > 0) as i32);
}
