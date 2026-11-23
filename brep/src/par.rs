//! Deterministic parallel map: work-stealing over indices, results in index order.
use std::sync::atomic::{AtomicUsize, Ordering};

pub fn map<T: Send>(n: usize, workers: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let next = AtomicUsize::new(0);
    let mut parts: Vec<(usize, T)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..workers.clamp(1, n.max(1))).map(|_| sc.spawn(|| {
            let mut out = vec![];
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n { break out; }
                out.push((i, f(i)));
            }
        })).collect();
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    parts.sort_unstable_by_key(|p| p.0);
    parts.into_iter().map(|p| p.1).collect()
}
