//! The two pieces of CPython (3.8+, 64-bit) that RDKit's Python
//! `EnumerateStereoisomers` uses to sample isomers when there are more flip
//! combinations than `maxIsomers`: `hash()` of a tuple of small-int pairs
//! (its default seed) and `random.Random(seed).getrandbits(k)` (MT19937
//! seeded by `init_by_array`).

const XXPRIME_1: u64 = 11400714785074694791;
const XXPRIME_2: u64 = 14029467366897019727;
const XXPRIME_5: u64 = 2870177450012600261;

/// CPython `hash(int)` for an `i64` (exact for |x| < 2^61 - 1).
fn int_hash(x: i64) -> u64 {
    const MODULUS: u64 = (1 << 61) - 1;
    let mut h = (x.unsigned_abs() % MODULUS) as i64;
    if x < 0 {
        h = -h;
    }
    if h == -1 {
        h = -2;
    }
    h as u64
}

/// CPython `tuplehash` over item hashes (`Py_uhash_t` lanes).
fn tuple_hash(lanes: impl ExactSizeIterator<Item = u64>) -> i64 {
    let len = lanes.len() as u64;
    let mut acc = XXPRIME_5;
    for lane in lanes {
        acc = acc.wrapping_add(lane.wrapping_mul(XXPRIME_2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(XXPRIME_1);
    }
    acc = acc.wrapping_add(len ^ (XXPRIME_5 ^ 3527539));
    if acc == u64::MAX {
        return 1546275796;
    }
    acc as i64
}

/// CPython `hash(tuple(pairs))` where `pairs` are 2-tuples of ints.
pub(crate) fn hash_pair_tuple(pairs: &[(i64, i64)]) -> i64 {
    tuple_hash(
        pairs
            .iter()
            .map(|&(a, b)| tuple_hash([int_hash(a), int_hash(b)].into_iter()) as u64),
    )
}

/// CPython's `random.Random` Mersenne Twister.
pub(crate) struct PyRandom {
    mt: [u32; 624],
    mti: usize,
}

impl PyRandom {
    fn init_genrand(s: u32) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = s;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        PyRandom { mt, mti: 624 }
    }

    /// `random.Random(seed)` for an int seed: `init_by_array` over the
    /// 32-bit little-endian words of `abs(seed)`.
    pub(crate) fn from_int_seed(seed: i64) -> Self {
        let n = seed.unsigned_abs();
        let key: Vec<u32> = if n >> 32 != 0 {
            vec![n as u32, (n >> 32) as u32]
        } else {
            vec![n as u32]
        };
        let mut r = Self::init_genrand(19650218);
        let mt = &mut r.mt;
        let (mut i, mut j) = (1usize, 0usize);
        let mut k = 624.max(key.len());
        while k > 0 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
            k -= 1;
        }
        k = 623;
        while k > 0 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1566083941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
            k -= 1;
        }
        mt[0] = 0x8000_0000;
        r
    }

    fn genrand_uint32(&mut self) -> u32 {
        const MAG01: [u32; 2] = [0, 0x9908_b0df];
        if self.mti >= 624 {
            let mt = &mut self.mt;
            for kk in 0..624 {
                let y = (mt[kk] & 0x8000_0000) | (mt[(kk + 1) % 624] & 0x7fff_ffff);
                mt[kk] = mt[(kk + 397) % 624] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
            }
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// `getrandbits(k)` as little-endian 32-bit words (`k >= 1`).
    pub(crate) fn getrandbits(&mut self, k: usize) -> Vec<u32> {
        if k <= 32 {
            return vec![self.genrand_uint32() >> (32 - k)];
        }
        let words = (k - 1) / 32 + 1;
        let mut out = Vec::with_capacity(words);
        let mut left = k;
        for _ in 0..words {
            let mut r = self.genrand_uint32();
            if left < 32 {
                r >>= 32 - left;
            }
            out.push(r);
            left = left.saturating_sub(32);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_u128(w: &[u32]) -> u128 {
        w.iter()
            .enumerate()
            .map(|(i, &x)| (x as u128) << (32 * i))
            .sum()
    }

    /// Values from CPython 3.13.
    #[test]
    fn matches_cpython() {
        let t = [(1, 6), (1, 8), (1, 17), (2, 6), (3, 6), (4, 6)];
        let h = hash_pair_tuple(&t);
        assert_eq!(h, -2720979148460318484);
        assert_eq!(tuple_hash(std::iter::empty()), 5740354900026072187);
        assert_eq!(hash_pair_tuple(&[(0, 0)]), 2983095501974941876);
        let mut r = PyRandom::from_int_seed(h);
        let first: Vec<u128> = (0..5).map(|_| to_u128(&r.getrandbits(11))).collect();
        assert_eq!(first, [808, 10, 305, 1658, 118]);
        assert_eq!(to_u128(&r.getrandbits(40)), 548527612548);
        assert_eq!(to_u128(&r.getrandbits(64)), 8176386211759365718);
        assert_eq!(to_u128(&r.getrandbits(70)), 1175514432062034866642);
        let mut z = PyRandom::from_int_seed(0);
        let v: Vec<u128> = (0..3).map(|_| to_u128(&z.getrandbits(32))).collect();
        assert_eq!(v, [3626764237, 1654615998, 3255389356]);
    }
}
