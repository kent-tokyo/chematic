//! Shared folding for RDKit Atom Pair and Topological Torsion counts.

use crate::bitvec::{BitVec2048, BitVecN};

pub(crate) const COUNT_BOUNDS: [u32; 4] = [1, 2, 4, 8];
pub(crate) const MAX_CONFIGURED_BITS: usize = 1_048_576;

pub(crate) fn valid_count_bounds(bounds: Option<&[u32]>, nbits: usize) -> bool {
    bounds.is_none_or(|bounds| {
        !bounds.is_empty()
            && bounds.len() < nbits
            && bounds.iter().all(|&bound| bound <= i32::MAX as u32)
    })
}

pub(crate) fn folded_elements(counts: Vec<u32>) -> Vec<(u32, u32)> {
    counts
        .into_iter()
        .enumerate()
        .filter(|&(_, count)| count > 0)
        .map(|(bucket, count)| (bucket as u32, count))
        .collect()
}

/// Pack the fixed-width default thresholds [1, 2, 4, 8].
pub(crate) fn default_count_bits(counts: &[u32]) -> BitVec2048 {
    debug_assert_eq!(counts.len(), 512);
    const BITS_BY_COUNT: [u64; 9] = [0, 1, 3, 3, 7, 7, 7, 7, 15];
    let mut words = [0u64; 32];
    for (word, buckets) in words.iter_mut().zip(counts.as_chunks::<16>().0) {
        for (position, &count) in buckets.iter().enumerate() {
            *word |= BITS_BY_COUNT[count.min(8) as usize] << (4 * position);
        }
    }
    BitVec2048::from_words(words)
}

pub(crate) struct ConfiguredCounts<'a> {
    folded: Vec<u32>,
    simulated: Option<Vec<u32>>,
    bounds: &'a [u32],
}

impl<'a> ConfiguredCounts<'a> {
    pub(crate) fn new(nbits: usize, bounds: Option<&'a [u32]>) -> Self {
        let bounds = bounds.unwrap_or(&[1]);
        Self {
            folded: vec![0; nbits],
            simulated: (bounds.len() != 1).then(|| vec![0; nbits / bounds.len()]),
            bounds,
        }
    }

    pub(crate) fn add(&mut self, hash: u32) {
        let hash = hash as usize;
        let width = self.folded.len();
        let count = &mut self.folded[hash % width];
        *count = count.saturating_add(1);
        if let Some(simulated) = &mut self.simulated {
            let width = simulated.len();
            let count = &mut simulated[hash % width];
            *count = count.saturating_add(1);
        }
    }

    pub(crate) fn finish(self) -> (BitVecN, Vec<(u32, u32)>) {
        let mut bits = BitVecN::new(self.folded.len());
        let counts = self.simulated.as_deref().unwrap_or(&self.folded);
        for (bucket, &count) in counts.iter().enumerate().filter(|&(_, &count)| count > 0) {
            for (i, &bound) in self.bounds.iter().enumerate() {
                if count >= bound {
                    bits.set(bucket * self.bounds.len() + i);
                }
            }
        }
        (bits, folded_elements(self.folded))
    }
}
