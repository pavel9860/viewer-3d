//! Shared allocation-free helpers (deterministic PRNG, content hash). Also compiled for wasm32 later.
pub struct Rng(pub u64);
impl Rng {
    /// splitmix64
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// uniform in [-1, 1)
    pub fn f(&mut self) -> f64 { (self.next() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) * 2.0 - 1.0 }
}
