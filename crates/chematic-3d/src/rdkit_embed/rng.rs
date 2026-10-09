//! RDKit's random source: `boost::minstd_rand` driven through
//! `boost::uniform_real<>(0, 1)` (`RDKit::double_source_type`).

/// `boost::random::minstd_rand` (Park–Miller, multiplier 48271).
#[derive(Clone, Debug)]
pub(crate) struct MinstdRand {
    x: u32,
}

const MODULUS: u64 = 2_147_483_647;
const MULTIPLIER: u64 = 48_271;

impl MinstdRand {
    /// `rng_type gen(42u); gen.seed(seed)`.
    pub(crate) fn new(seed: u32) -> Self {
        let mut x = (seed as u64 % MODULUS) as u32;
        if x == 0 {
            x = 1;
        }
        MinstdRand { x }
    }

    fn next_u32(&mut self) -> u32 {
        self.x = ((self.x as u64 * MULTIPLIER) % MODULUS) as u32;
        self.x
    }

    /// One draw of `boost::uniform_real<double>(0.0, 1.0)`.
    pub(crate) fn next_f64(&mut self) -> f64 {
        // generate_uniform_real: numerator = eng() - min, divisor =
        // (max - min) + 1 with min = 1, max = 2^31 - 2; retry when the
        // result reaches max_value (cannot happen here).
        loop {
            let numerator = (self.next_u32() - 1) as f64;
            let divisor = (2_147_483_646u32 - 1) as f64 + 1.0;
            let result = numerator / divisor * (1.0 - 0.0) + 0.0;
            if result < 1.0 {
                return result;
            }
        }
    }
}
