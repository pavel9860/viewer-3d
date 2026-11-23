//! B-rep topology.
//!
//! Entities live in chunked copy-on-write arenas: cloning a `Body` copies only chunk handles, and an edit copies
//! just the chunks it touches (structural sharing between versions; immutable bodies are safe to share across
//! threads). Ids are typed indices. Every face and edge carries its own tolerance and a provenance record
//! (which operation created it, from which input entity), which persistent naming builds on.
use crate::geom::{Surf, V};
use crate::nurbs::Crv;
use std::marker::PhantomData;
use std::sync::Arc;

const CHUNK: usize = 256;

/// typed index into an arena
pub struct Id<T>(pub u32, PhantomData<fn() -> T>);
impl<T> Id<T> { pub fn new(i: usize) -> Self { Id(i as u32, PhantomData) } pub fn ix(self) -> usize { self.0 as usize } }
impl<T> Clone for Id<T> { fn clone(&self) -> Self { *self } }
impl<T> Copy for Id<T> {}
impl<T> PartialEq for Id<T> { fn eq(&self, o: &Self) -> bool { self.0 == o.0 } }
impl<T> Eq for Id<T> {}
impl<T> std::hash::Hash for Id<T> { fn hash<H: std::hash::Hasher>(&self, h: &mut H) { self.0.hash(h) } }
impl<T> std::fmt::Debug for Id<T> { fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result { write!(f, "#{}", self.0) } }

/// chunked copy-on-write arena
pub struct Arena<T> { chunks: Vec<Arc<Vec<T>>>, len: usize }
impl<T> Default for Arena<T> { fn default() -> Self { Arena { chunks: vec![], len: 0 } } }
impl<T> Clone for Arena<T> { fn clone(&self) -> Self { Arena { chunks: self.chunks.clone(), len: self.len } } }
impl<T: Clone> Arena<T> {
    pub fn len(&self) -> usize { self.len }
    pub fn is_empty(&self) -> bool { self.len == 0 }
    pub fn push(&mut self, v: T) -> Id<T> {
        if self.len % CHUNK == 0 { self.chunks.push(Arc::new(Vec::with_capacity(CHUNK))); }
        Arc::make_mut(self.chunks.last_mut().unwrap()).push(v);
        self.len += 1;
        Id::new(self.len - 1)
    }
    pub fn get(&self, id: Id<T>) -> &T { &self.chunks[id.ix() / CHUNK][id.ix() % CHUNK] }
    /// mutable access; copies the chunk first if another body version shares it
    pub fn get_mut(&mut self, id: Id<T>) -> &mut T { &mut Arc::make_mut(&mut self.chunks[id.ix() / CHUNK])[id.ix() % CHUNK] }
    pub fn ids(&self) -> impl Iterator<Item = Id<T>> { (0..self.len).map(Id::new) }
    pub fn iter(&self) -> impl Iterator<Item = (Id<T>, &T)> { self.ids().map(move |i| (i, self.get(i))) }
    /// number of chunks physically shared with `other` (structural sharing between versions)
    pub fn shared_chunks(&self, other: &Self) -> usize { self.chunks.iter().zip(&other.chunks).filter(|(a, b)| Arc::ptr_eq(a, b)).count() }
}

/// which operation created an entity, and from which input entity (0 = none)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Origin { pub op: u32, pub src: u32 }

#[derive(Clone)] pub struct Vertex { pub p: V, pub tol: f64 }
#[derive(Clone)] pub struct Edge { pub curve: Option<Arc<dyn Crv>>, pub v: [Id<Vertex>; 2], pub t: [f64; 2], pub tol: f64, pub origin: Origin }
/// use of an edge by a loop; `sense` = traversed along the edge direction
#[derive(Clone)] pub struct Coedge { pub edge: Id<Edge>, pub sense: bool }
#[derive(Clone)] pub struct Loop { pub coedges: Vec<Id<Coedge>> }
/// `sense` = face normal agrees with the surface normal
#[derive(Clone)] pub struct Face { pub surface: Arc<dyn Surf>, pub sense: bool, pub loops: Vec<Id<Loop>>, pub tol: f64, pub origin: Origin }
#[derive(Clone)] pub struct Shell { pub faces: Vec<Id<Face>>, pub closed: bool }

/// a body: shells of faces; cheap to clone, edits copy only touched chunks
#[derive(Clone, Default)]
pub struct Body {
    pub vertices: Arena<Vertex>,
    pub edges: Arena<Edge>,
    pub coedges: Arena<Coedge>,
    pub loops: Arena<Loop>,
    pub faces: Arena<Face>,
    pub shells: Arena<Shell>,
    /// model tolerance (model units)
    pub tol: f64,
}
impl Body {
    /// start and end vertex of a coedge in traversal direction
    pub fn coedge_ends(&self, c: Id<Coedge>) -> (Id<Vertex>, Id<Vertex>) {
        let ce = self.coedges.get(c);
        let e = self.edges.get(ce.edge);
        if ce.sense { (e.v[0], e.v[1]) } else { (e.v[1], e.v[0]) }
    }
    pub fn stats(&self) -> String {
        format!("{} shells, {} faces, {} loops, {} edges, {} vertices", self.shells.len(), self.faces.len(), self.loops.len(), self.edges.len(), self.vertices.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn clone_shares_chunks_and_edit_copies_only_the_touched_chunk() {
        let mut a: Arena<Vertex> = Arena::default();
        for i in 0..10 * CHUNK { a.push(Vertex { p: [i as f64, 0.0, 0.0], tol: 0.0 }); }
        let mut b = a.clone();
        assert_eq!(b.shared_chunks(&a), 10);
        b.get_mut(Id::new(3 * CHUNK + 5)).p[1] = 1.0;
        assert_eq!(b.shared_chunks(&a), 9);
        assert_eq!(a.get(Id::new(3 * CHUNK + 5)).p[1], 0.0);
        assert_eq!(b.get(Id::new(3 * CHUNK + 5)).p[1], 1.0);
    }
}
