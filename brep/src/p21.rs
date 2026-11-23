//! Minimal ISO 10303-21 (STEP physical file) parser: DATA section, simple and complex entity instances.
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub enum Val { Num(f64), Str, Enum(String), Ref(u32), List(Vec<Val>), Typed, Null }
#[derive(Debug)]
pub enum Ent { Simple(String, Vec<Val>), Complex(Vec<(String, Vec<Val>)>) }

struct P<'a> { b: &'a [u8], i: usize }
impl P<'_> {
    fn ws(&mut self) { while self.i < self.b.len() && (self.b[self.i] as char).is_whitespace() { self.i += 1; } }
    fn val(&mut self) -> Val {
        self.ws();
        match self.b[self.i] {
            b'#' => { self.i += 1; let s = self.i; while self.b[self.i].is_ascii_digit() { self.i += 1; } Val::Ref(std::str::from_utf8(&self.b[s..self.i]).unwrap().parse().unwrap()) }
            b'\'' => { self.i += 1; loop { if self.b[self.i] == b'\'' { if self.b.get(self.i + 1) == Some(&b'\'') { self.i += 2; continue; } self.i += 1; break; } self.i += 1; } Val::Str }
            b'.' => { self.i += 1; let s = self.i; while self.b[self.i] != b'.' { self.i += 1; } let e = std::str::from_utf8(&self.b[s..self.i]).unwrap().to_string(); self.i += 1; Val::Enum(e) }
            b'(' => Val::List(self.list()),
            b'$' | b'*' => { self.i += 1; Val::Null }
            b'"' => { self.i += 1; while self.b[self.i] != b'"' { self.i += 1; } self.i += 1; Val::Str }
            c if c == b'-' || c == b'+' || c.is_ascii_digit() || c == b'.' => {
                let s = self.i; self.i += 1;
                while self.i < self.b.len() && (self.b[self.i].is_ascii_digit() || b".eE+-".contains(&self.b[self.i])) { self.i += 1; }
                Val::Num(std::str::from_utf8(&self.b[s..self.i]).unwrap().parse().unwrap_or(0.0))
            }
            _ => { let _ = self.name(); self.ws(); let _ = self.list(); Val::Typed }
        }
    }
    fn name(&mut self) -> String { self.ws(); let s = self.i; while self.i < self.b.len() && (self.b[self.i].is_ascii_alphanumeric() || self.b[self.i] == b'_') { self.i += 1; } std::str::from_utf8(&self.b[s..self.i]).unwrap().to_string() }
    fn list(&mut self) -> Vec<Val> {
        self.ws(); assert_eq!(self.b[self.i], b'('); self.i += 1;
        let mut v = vec![];
        loop { self.ws(); if self.b[self.i] == b')' { self.i += 1; break; } v.push(self.val()); self.ws(); if self.b[self.i] == b',' { self.i += 1; } }
        v
    }
}

/// entity instances of the DATA section by id
pub fn parse(raw: &[u8]) -> HashMap<u32, Ent> {
    // strip comments outside strings
    let mut b = Vec::with_capacity(raw.len());
    let (mut i, mut in_str) = (0, false);
    while i < raw.len() {
        if !in_str && raw[i] == b'/' && raw.get(i + 1) == Some(&b'*') { i += 2; while i + 1 < raw.len() && !(raw[i] == b'*' && raw[i + 1] == b'/') { i += 1; } i += 2; continue; }
        if raw[i] == b'\'' { in_str = !in_str; }
        b.push(raw[i]); i += 1;
    }
    let s = b.windows(5).position(|w| w == b"DATA;").unwrap() + 5;
    let mut p = P { b: &b, i: s };
    let mut m = HashMap::new();
    loop {
        p.ws();
        if p.b[p.i..].starts_with(b"ENDSEC") { break; }
        let Val::Ref(id) = p.val() else { panic!("entity id expected") };
        p.ws(); assert_eq!(p.b[p.i], b'='); p.i += 1; p.ws();
        let e = if p.b[p.i] == b'(' {
            p.i += 1;
            let mut parts = vec![];
            loop { p.ws(); if p.b[p.i] == b')' { p.i += 1; break; } let n = p.name(); let a = p.list(); parts.push((n, a)); }
            Ent::Complex(parts)
        } else { let n = p.name(); let a = p.list(); Ent::Simple(n, a) };
        p.ws(); assert_eq!(p.b[p.i], b';'); p.i += 1;
        m.insert(id, e);
    }
    m
}
pub fn num(v: &Val) -> f64 { if let Val::Num(x) = v { *x } else { panic!("number expected: {v:?}") } }
pub fn rf(v: &Val) -> u32 { if let Val::Ref(x) = v { *x } else { panic!("reference expected: {v:?}") } }
pub fn lst(v: &Val) -> &Vec<Val> { if let Val::List(x) = v { x } else { panic!("list expected: {v:?}") } }
