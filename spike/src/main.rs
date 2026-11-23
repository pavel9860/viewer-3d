//! Kernel spikes. Usage:
//!   spike ssi                 hand-made SSI cases
//!   spike stress [TRIALS]     randomized SSI stress (default 60 per family)
//!   spike step LIST [opts]    CATIA STEP edge check (see step.rs)
mod geom;
mod k;
mod ssi;
mod ssi_cases;
mod nurbs;
mod p21;
mod step;

use std::time::Instant;
pub fn ms(t: Instant) -> f64 { t.elapsed().as_secs_f64() * 1e3 }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let code = match args.get(1).map(|s| s.as_str()) {
        Some("ssi") => ssi_cases::cases(),
        Some("stress") => ssi_cases::stress(args.get(2).and_then(|v| v.parse().ok()).unwrap_or(60)),
        Some("step") => step::run(&args[2..]),
        _ => { eprintln!("usage: spike ssi | stress [N] | step LIST [--full] [--pair KIND] [--file NAME]"); 2 }
    };
    std::process::exit(code.min(1) as i32);
}
