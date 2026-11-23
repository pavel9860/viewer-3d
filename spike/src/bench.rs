//! Architecture benchmarks: thread scaling + determinism reference (compared with the browser), worker transfer,
//! incremental recompute, content hashing, persistent naming, assembly instancing (BVH, pick, cull).
use crate::ms;
use kcore::*;
use std::collections::HashMap;
use std::hint::black_box;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::time::Instant;

pub const WORKERS: usize = 2;
pub fn as_bytes<T>(v: &[T]) -> &[u8] { unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) } }

/// run `f(i)` for i in 0..n on `t` threads (work stealing by atomic counter)
pub fn par_for(n: usize, t: usize, f: impl Fn(usize) + Sync) {
    let next = AtomicUsize::new(0);
    std::thread::scope(|s| for _ in 0..t { s.spawn(|| loop { let i = next.fetch_add(1, Relaxed); if i >= n { break; } f(i); }); });
}

/// thread scaling and the determinism reference hashes (identical code runs in the browser)
pub fn determinism() {
    println!("== thread scaling + determinism reference (kcore) ==");
    let mut base = 0.0;
    for t in 1..=WORKERS {
        let res: Vec<AtomicU64> = (0..NT).map(|_| AtomicU64::new(0)).collect();
        let t0 = Instant::now();
        par_for(NT, t, |i| res[i].store(patch_work(i as u64, PN), Relaxed));
        let e = ms(t0);
        if t == 1 { base = e; }
        let v: Vec<u64> = res.iter().map(|x| x.load(Relaxed)).collect();
        println!("  threads={t} {e:7.1} ms speedup {:.2}x  patch hash {:016x}", base / e, hash_bytes(as_bytes(&v)));
    }
    let mut x = vec![0f64; SUM_N];
    fill(&mut x, 42);
    let s = x.chunks(SUM_CH).map(sum_seq).fold(0.0, |a, b| a + b);
    println!("  chunked sum {:016x} | fill {:016x} | orient2d {:016x} | own sin {:016x}", s.to_bits(), hash_bytes(as_bytes(&x)), orient_hash(7, 1_000_000), sin_hash(9, 1_000_000));
    let mut r = Rng(9);
    let (mut diff, mut maxulp) = (0, 0i64);
    for _ in 0..1_000_000 { let v = r.f() * 1000.0; let (a, b) = (v.sin(), own_sin(v)); if a.to_bits() != b.to_bits() { diff += 1; maxulp = maxulp.max((a.to_bits() as i64 - b.to_bits() as i64).abs()); } }
    println!("  platform sin vs own sin: {diff} / 1M differ (max {maxulp} ulp) -> own math functions are required for determinism");
}

/// worker transfer of a typical part body (5K faces: 32 B topology + 512 B geometry)
pub fn transfer() {
    use std::sync::{mpsc, Arc};
    println!("== worker transfer (5K-face body, 2.72 MB) ==");
    let n = 5000;
    let topo: Vec<[u64; 4]> = (0..n).map(|i| [i as u64; 4]).collect();
    let geom: Vec<[f64; 64]> = (0..n).map(|i| [i as f64; 64]).collect();
    let body = Arc::new((topo, geom));
    let it = 200;
    let round_trip = |send: &dyn Fn(&mpsc::Sender<Box<dyn FnOnce() -> f64 + Send>>)| {
        let (tx, rx) = mpsc::channel::<Box<dyn FnOnce() -> f64 + Send>>();
        let (btx, brx) = mpsc::channel::<f64>();
        let h = std::thread::spawn(move || for job in rx { btx.send(job()).unwrap(); });
        let t0 = Instant::now();
        for _ in 0..it { send(&tx); black_box(brx.recv().unwrap()); }
        let e = ms(t0) / it as f64;
        drop(tx); h.join().unwrap();
        e
    };
    let shared = round_trip(&|tx| { let b = body.clone(); tx.send(Box::new(move || b.1[n - 1][0])).unwrap(); });
    let copy = round_trip(&|tx| { let b = (body.0.clone(), body.1.clone()); tx.send(Box::new(move || b.1[n - 1][0])).unwrap(); });
    let ser = round_trip(&|tx| {
        let mut v = Vec::with_capacity(n * 544);
        v.extend_from_slice(as_bytes(&body.0)); v.extend_from_slice(as_bytes(&body.1));
        tx.send(Box::new(move || { let mut g = vec![[0f64; 64]; n]; unsafe { std::ptr::copy_nonoverlapping(v[n * 32..].as_ptr(), g.as_mut_ptr() as *mut u8, n * 512) }; g[n - 1][0] })).unwrap();
    });
    println!("  shared handle {shared:.4} ms | copy {copy:.3} ms | serialize+deserialize {ser:.3} ms (round trip)");
}

/// incremental recompute over a 50-feature dependency graph with a content-hash cache
pub fn incremental() {
    println!("== incremental recompute (50 features, ~1 ms each) ==");
    let nf = 50;
    let mut r = Rng(3);
    let deps: Vec<Vec<usize>> = (0..nf).map(|i| { let mut d = vec![]; if i > 0 { d.push(i - 1); } if i > 5 && i % 5 == 0 { d.push((r.next() % (i as u64 - 2)) as usize); } d }).collect();
    let mut params: Vec<u64> = (0..nf as u64).collect();
    let mut cache: HashMap<u64, u64> = HashMap::new();
    let build = |params: &[u64], cache: &mut HashMap<u64, u64>| {
        let mut out = vec![0u64; nf];
        let mut computed = 0;
        for i in 0..nf {
            let key = deps[i].iter().fold(mix(0xF00D, params[i]), |k, &d| mix(k, out[d]));
            out[i] = *cache.entry(key).or_insert_with(|| { computed += 1; patch_work(key, 180) });
        }
        computed
    };
    let t0 = Instant::now();
    println!("  full build {:.1} ms ({} computed)", { build(&params, &mut cache); ms(t0) }, nf);
    for f in [0usize, 25, 49] {
        params[f] += 1000;
        let t0 = Instant::now();
        let c = build(&params, &mut cache);
        println!("  edit feature {f:2}: {:6.1} ms, {c:2} recomputed", ms(t0));
    }
    params[25] -= 1000;
    let t0 = Instant::now();
    let c = build(&params, &mut cache);
    println!("  undo (old value): {:.2} ms, {c} recomputed", ms(t0));
    let worst = (0..30).map(|k| { params[40] = 5000 + k; let t0 = Instant::now(); build(&params, &mut cache); ms(t0) }).fold(0.0, f64::max);
    println!("  drag feature 40: worst frame {worst:.1} ms (budget 16.7 ms) -> {}", if worst < 16.7 { "OK" } else { "TOO SLOW" });
}

pub fn hashing() {
    println!("== content hashing ==");
    let mut x = vec![0f64; 16 << 20];
    fill(&mut x, 1);
    let b = as_bytes(&x);
    let t0 = Instant::now();
    let h = black_box(hash_bytes(b));
    let e = ms(t0);
    println!("  {} MB in {e:.1} ms = {:.2} GB/s ({h:016x}); 2.72 MB part ~{:.3} ms", b.len() >> 20, b.len() as f64 / e / 1e6, 2.72e6 / (b.len() as f64 / e));
}

/// persistent naming mechanism: provenance names vs face indices under random program edits (toy topology)
pub fn naming() {
    type Name = Vec<(u32, u8)>;
    #[derive(Clone)]
    enum Op { Box(u32), Split(u32, Name, bool), Hole(u32, Name) }
    let id = |o: &Op| match o { Op::Box(i) | Op::Split(i, ..) | Op::Hole(i, _) => *i };
    let target = |o: &Op| match o { Op::Box(_) => None, Op::Split(_, t, _) | Op::Hole(_, t) => Some(t.clone()) };
    enum Res { Found(usize), Ambiguous, Missing }
    let resolve = |n: &Name, f: &[Name]| -> Res {
        if let Some(i) = f.iter().position(|x| x == n) { return Res::Found(i); }
        let d: Vec<usize> = (0..f.len()).filter(|&i| f[i].starts_with(n)).collect();
        match d.len() { 0 => Res::Missing, 1 => Res::Found(d[0]), _ => Res::Ambiguous }
    };
    let build = |p: &[Op]| -> Vec<Name> {
        let mut f: Vec<Name> = vec![];
        for o in p {
            match o {
                Op::Box(i) => f.extend((0..6).map(|r| vec![(*i, r)])),
                Op::Split(i, t, on) => if let (Res::Found(k), true) = (resolve(t, &f), *on) { let base = f[k].clone(); f[k] = [base.clone(), vec![(*i, 0)]].concat(); f.push([base, vec![(*i, 1)]].concat()); },
                Op::Hole(i, t) => if let Res::Found(_) = resolve(t, &f) { f.push(vec![(*i, 2)]); },
            }
        }
        f
    };
    println!("== persistent naming: provenance vs index, 100 cumulative random edits ==");
    let mut r = Rng(5);
    let mut next = 1;
    let mut prog = vec![Op::Box(0)];
    for k in 0..15 { let f = build(&prog); let t = f[(r.next() % f.len() as u64) as usize].clone(); prog.push(if k % 3 == 2 { Op::Hole(next, t) } else { Op::Split(next, t, true) }); next += 1; }
    let refs: Vec<(Name, usize)> = build(&prog).into_iter().enumerate().map(|(i, n)| (n, i)).collect();
    let (mut ok, mut wrong, mut gone, mut found, mut amb, mut miss) = (0, 0, 0, 0, 0, 0);
    for _ in 0..100 {
        match r.next() % 4 {
            0 => { let s: Vec<usize> = (0..prog.len()).filter(|&i| matches!(prog[i], Op::Split(..))).collect(); if !s.is_empty() { let i = s[(r.next() % s.len() as u64) as usize]; if let Op::Split(_, _, on) = &mut prog[i] { *on = r.next() % 3 != 0; } } }
            1 => { let at = 1 + (r.next() % prog.len() as u64) as usize; let f = build(&prog[..at]); let t = f[(r.next() % f.len() as u64) as usize].clone(); prog.insert(at, Op::Split(next, t, true)); next += 1; }
            2 => if prog.len() > 2 { let i = 1 + (r.next() % (prog.len() - 1) as u64) as usize; prog.remove(i); },
            _ => if prog.len() > 3 {
                let i = 1 + (r.next() % (prog.len() - 2) as u64) as usize;
                let uses = |o: &Op, by: u32| target(o).is_some_and(|t| t.iter().any(|x| x.0 == by));
                if !uses(&prog[i + 1], id(&prog[i])) && !uses(&prog[i], id(&prog[i + 1])) { prog.swap(i, i + 1); }
            },
        }
        let f = build(&prog);
        for (n, i) in &refs {
            let truth = resolve(n, &f);
            match truth { Res::Found(_) => found += 1, Res::Ambiguous => amb += 1, Res::Missing => miss += 1 }
            match f.get(*i) { None => gone += 1, Some(x) => if x.starts_with(n) { ok += 1 } else { wrong += 1 } }
        }
    }
    let tot = (ok + wrong + gone) as f64;
    println!("  index naming: correct {:.1}%  SILENTLY WRONG {:.1}%  missing {:.1}%", 100.0 * ok as f64 / tot, 100.0 * wrong as f64 / tot, 100.0 * gone as f64 / tot);
    println!("  provenance:   resolved {:.1}%  flagged ambiguous {:.1}%  flagged missing {:.1}%  silently wrong 0 (by construction)", 100.0 * found as f64 / tot, 100.0 * amb as f64 / tot, 100.0 * miss as f64 / tot);
}

/// axis-aligned bounding-volume hierarchy over boxes (median split, pre-order: left child = id + 1)
pub struct Bvh { pub lo: Vec<[f32; 3]>, pub hi: Vec<[f32; 3]>, pub right: Vec<u32>, pub first: Vec<u32>, pub count: Vec<u32>, pub items: Vec<u32> }
impl Bvh {
    pub fn build(lo: &[[f32; 3]], hi: &[[f32; 3]], leaf: usize) -> Bvh {
        let mut b = Bvh { lo: vec![], hi: vec![], right: vec![], first: vec![], count: vec![], items: (0..lo.len() as u32).collect() };
        b.node(lo, hi, 0, lo.len(), leaf);
        b
    }
    fn node(&mut self, lo: &[[f32; 3]], hi: &[[f32; 3]], a: usize, e: usize, leaf: usize) -> usize {
        let id = self.lo.len();
        let (mut l, mut h) = ([f32::MAX; 3], [f32::MIN; 3]);
        for &i in &self.items[a..e] { for k in 0..3 { l[k] = l[k].min(lo[i as usize][k]); h[k] = h[k].max(hi[i as usize][k]); } }
        self.lo.push(l); self.hi.push(h); self.right.push(0); self.first.push(a as u32); self.count.push((e - a) as u32);
        if e - a <= leaf { return id; }
        let ax = (0..3).max_by(|&x, &y| (h[x] - l[x]).total_cmp(&(h[y] - l[y]))).unwrap();
        let m = (a + e) / 2;
        self.items[a..e].select_nth_unstable_by(m - a, |&x, &y| (lo[x as usize][ax] + hi[x as usize][ax]).total_cmp(&(lo[y as usize][ax] + hi[y as usize][ax])));
        self.count[id] = 0;
        self.node(lo, hi, a, m, leaf);
        self.right[id] = self.node(lo, hi, m, e, leaf) as u32;
        id
    }
    /// visit leaf items whose box passes `test`
    pub fn query(&self, test: &impl Fn([f32; 3], [f32; 3]) -> bool, mut f: impl FnMut(u32)) {
        let mut st = vec![0usize];
        while let Some(n) = st.pop() {
            if !test(self.lo[n], self.hi[n]) { continue; }
            if self.count[n] > 0 { for &i in &self.items[self.first[n] as usize..(self.first[n] + self.count[n]) as usize] { f(i); } } else { st.push(self.right[n] as usize); st.push(n + 1); }
        }
    }
}
pub fn slab(lo: [f32; 3], hi: [f32; 3], o: [f32; 3], inv: [f32; 3]) -> f32 {
    let (mut t0, mut t1) = (0f32, f32::MAX);
    for k in 0..3 { let (a, b) = ((lo[k] - o[k]) * inv[k], (hi[k] - o[k]) * inv[k]); t0 = t0.max(a.min(b)); t1 = t1.min(a.max(b)); }
    if t0 <= t1 { t0 } else { f32::MAX }
}

/// assembly scale: 10K instances of 200 unique parts; instance BVH, culling checked against brute force
pub fn assembly() {
    println!("== assembly instancing (10K instances, 200 unique parts) ==");
    #[derive(Clone, Copy)]
    struct Inst { part: u32, rot: [f32; 9], t: [f32; 3], s: f32 }
    let mut r = Rng(9);
    let mut rf = || r.f() as f32;
    let part_box: Vec<([f32; 3], [f32; 3])> = (0..200).map(|_| ([-1.0 + 0.1 * rf(); 3], [1.0 + 0.1 * rf(); 3])).collect();
    let inst: Vec<Inst> = (0..10_000).map(|_| { let a = 3.14159 * rf(); let (s, c) = a.sin_cos(); Inst { part: ((rf() + 1.0) * 99.9) as u32, rot: [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0], t: [rf() * 100.0, rf() * 100.0, rf() * 5.0], s: 1.0 + rf().abs() } }).collect();
    let world = |i: &Inst| {
        let (bl, bh) = part_box[i.part as usize];
        let (mut lo, mut hi) = ([f32::MAX; 3], [f32::MIN; 3]);
        for m in 0..8 {
            let p = [if m & 1 == 0 { bl[0] } else { bh[0] }, if m & 2 == 0 { bl[1] } else { bh[1] }, if m & 4 == 0 { bl[2] } else { bh[2] }];
            for k in 0..3 { let q = (i.rot[k * 3] * p[0] + i.rot[k * 3 + 1] * p[1] + i.rot[k * 3 + 2] * p[2]) * i.s + i.t[k]; lo[k] = lo[k].min(q); hi[k] = hi[k].max(q); }
        }
        (lo, hi)
    };
    let boxes: Vec<_> = inst.iter().map(world).collect();
    let (lo, hi): (Vec<_>, Vec<_>) = boxes.iter().copied().unzip();
    let t0 = Instant::now();
    let bvh = Bvh::build(&lo, &hi, 4);
    let tb = ms(t0);
    let view = |k: usize| { let (cx, cy, w) = ((k * 7 % 80) as f32, (k * 13 % 80) as f32, 30.0f32); move |l: [f32; 3], h: [f32; 3]| !(h[0] < cx - w || l[0] > cx + w || h[1] < cy - w || l[1] > cy + w) };
    let cull = |k: usize| { let inside = view(k); let mut n = 0; bvh.query(&inside, |i| if inside(lo[i as usize], hi[i as usize]) { n += 1; }); n };
    let t0 = Instant::now();
    let visible: usize = (0..100).map(cull).sum();
    let tc = ms(t0) / 100.0;
    let brute: usize = (0..20).map(|k| { let inside = view(k); boxes.iter().filter(|b| inside(b.0, b.1)).count() }).sum();
    let ok = (0..20).map(cull).sum::<usize>() == brute;
    println!("  {} B/instance -> {:.2} MB | BVH {} nodes in {tb:.2} ms | cull {tc:.3} ms/view (~{} visible) | BVH vs brute force, 20 views: {}", std::mem::size_of::<Inst>(), std::mem::size_of::<Inst>() as f64 * 1e4 / 1e6, bvh.lo.len(), visible / 100, if ok { "MATCH" } else { "MISMATCH" });
}
