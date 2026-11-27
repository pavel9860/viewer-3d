//! Point containment for a closed face set: the sign of (q - p) · n at the closest boundary point p, where n is
//! the face normal for a face-interior p and the sum of the adjacent face normals for an edge point (pseudonormal).
use crate::geom::*;
use crate::mesh::face_trim;
use crate::nurbs::tess;
use crate::topo::*;
use crate::trim::{Trim, SAG};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side { In, Out, On }

pub struct Solid<'a> {
    b: &'a Body,
    faces: Vec<(Id<Face>, Trim)>,
    /// edge polyline and the faces using the edge
    edges: Vec<(Vec<V>, Vec<Id<Face>>)>,
}

impl<'a> Solid<'a> {
    pub fn new(b: &'a Body, faces: &[Id<Face>]) -> Self {
        let mut users: HashMap<Id<Edge>, Vec<Id<Face>>> = HashMap::new();
        for &f in faces { for ce in b.face_coedges(f) { users.entry(ce.edge).or_default().push(f); } }
        let mut ids: Vec<_> = users.into_iter().collect();
        ids.sort_unstable_by_key(|x| x.0 .0);
        let edges = ids.into_iter().filter_map(|(e, fs)| {
            let ed = b.edges.get(e);
            let c = ed.curve.as_ref()?;
            Some((tess(c.as_ref(), ed.t[0], ed.t[1], SAG).into_iter().map(|t| c.ev(t).0).collect(), fs))
        }).collect();
        Solid { b, faces: faces.iter().map(|&f| (f, face_trim(b, f))).collect(), edges }
    }
    fn normal(&self, f: Id<Face>, p: V) -> V {
        let fc = self.b.faces.get(f);
        let uv = fc.surface.uv_of(p);
        let d = fc.surface.ev(uv[0], uv[1]);
        let n = unit(cross(d.u, d.v));
        if fc.sense { n } else { mulv(n, -1.0) }
    }
    /// side of q; `On` within `tol` of the boundary
    pub fn side(&self, q: V, tol: f64) -> Side {
        // (distance, outside?)
        let mut best = (f64::MAX, false);
        for (f, t) in &self.faces {
            let s = self.b.faces.get(*f).surface.as_ref();
            let p = { let uv = s.uv_of(q); s.ev(uv[0], uv[1]).p };
            let d = norm(sub(q, p));
            if d < best.0 && t.inside(t.uv(s, p)) { best = (d, dot(sub(q, p), self.normal(*f, p)) >= 0.0); }
        }
        for (pl, fs) in &self.edges {
            for w in pl.windows(2) {
                let ab = sub(w[1], w[0]);
                let u = (dot(sub(q, w[0]), ab) / dot(ab, ab).max(1e-300)).clamp(0.0, 1.0);
                let p = addv(w[0], mulv(ab, u));
                let d = norm(sub(q, p));
                if d < best.0 {
                    let n = fs.iter().fold([0.0; 3], |a, &f| addv(a, self.normal(f, p)));
                    best = (d, dot(sub(q, p), n) >= 0.0);
                }
            }
        }
        if best.0 <= tol { Side::On } else if best.1 { Side::Out } else { Side::In }
    }
}
