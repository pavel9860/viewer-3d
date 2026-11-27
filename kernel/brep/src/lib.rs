//! B-rep kernel core.

/// diagnostics switch: BREP_DEBUG unset = off; empty, "1" or "all" = every topic; else a comma list of topics
/// (check, bool, time, mesh, tess, sew, orient, open, nonman, gaps, faces, tri)
pub fn debug(topic: &str) -> bool {
    match std::env::var("BREP_DEBUG") {
        Err(_) => false,
        Ok(v) => v.is_empty() || v == "1" || v == "all" || v.split(',').any(|t| t == topic),
    }
}
pub mod boolean;
pub mod cgm;
pub mod cgm_mesh;
pub mod check;
pub mod classify;
pub mod exact;
pub mod geom;
pub mod make;
pub mod mesh;
pub mod nurbs;
pub mod p21;
pub mod sew;
pub mod simplify;
pub mod par;
pub mod ssi;
pub mod step;
pub mod tess;
pub mod topo;
pub mod tri;
pub mod trim;
