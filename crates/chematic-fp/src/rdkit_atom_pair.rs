//! RDKit-compatible Atom Pair fingerprint (`GetHashedAtomPairFingerprintAsBitVect`).
//!
//! Opt-in `RdkitExact`-mode addition, kept fully separate from
//! [`crate::atom_pair::atom_pair_fp`] (chematic's own native scheme). Second of
//! the fingerprint-parity series (Track A / 99-point directive Phase 6), after
//! [`crate::rdkit_torsion_fp`] -- reuses that module's atom-invariant primitive
//! ([`crate::rdkit_torsion::atom_code`]) since RDKit's own `AtomPairAtomInvGenerator`
//! is the shared invariant generator for both fingerprints, just constructed with
//! its `topologicalTorsionCorrection` flag `false` here instead of `true`.
//!
//! Derived from RDKit's C++ source (`Code/GraphMol/Fingerprints/AtomPairGenerator.{h,cpp}`,
//! `AtomPairGenerator.h`'s `AtomPairArguments` defaults, `FingerprintGenerator.cpp`'s
//! `getFingerprintHelper`), fetched at research time against `master`. Every
//! constant/formula below is a direct port -- see each function's doc comment for
//! the exact C++ source it mirrors.
//!
//! Key differences from [`crate::rdkit_torsion_fp`], confirmed against source
//! before implementing (not assumed from the superficial code-family similarity):
//! - No `-2` torsion correction on atom invariants -- `getAtomCode(atom, 0, false)`,
//!   unmodified.
//! - Atom codes are reduced by a plain `% ((1 << codeSize) - 1)`, with no `+1`
//!   offset and no interior-position decrement (that offset is
//!   `getTopologicalTorsionHash`-specific, not shared).
//! - The bit id is `hash_combine(hash_combine(hash_combine(0, min(codeI, codeJ)),
//!   dist), max(codeI, codeJ))` (`AtomPairAtomEnv::getBitId`'s `hashResults` branch,
//!   always taken here since `GetHashedAtomPairFingerprintAsBitVect` always sets a
//!   nonzero `fpSize`) -- a 3-element [`hash_vec`] fold, not the path-direction-
//!   canonicalizing scheme torsion uses.
//! - Pairs range over topological distance `[1, 30]` (`AtomPairArguments`'s
//!   `minDistance=1`, `maxDistance=maxPathLen-1=30`), each unordered pair counted
//!   once -- there is no 3-membered-ring-closure analogue (that was a torsion-path
//!   artifact of revisiting the pivot atom, not applicable to plain pairwise
//!   distance).
//!
//! The final 512-bucket count-simulation encoding (`nBitsPerEntry=4`, bounds
//! `[1,2,4,8]`, folded into a 2048-bit output) is identical in structure to
//! [`crate::rdkit_torsion_fp`]'s.
//!
//! Against RDKit 2026.09.1's generator, the pinned CLI comparison is
//! 4,072/4,072 exact on `stress-v1.smi`, 10,000/10,000 on the browser corpus,
//! and 5,000/5,000 on `descriptor_census_corpus.smi`. These fixed-corpus
//! results do not establish universal bit-exactness.

use std::fmt;
use std::ops::RangeInclusive;

use chematic_core::{AtomIdx, Molecule};
use rustc_hash::FxHashMap;

use crate::bitvec::{BitVec2048, BitVecN};
use crate::rdkit_morgan_hash::hash_vec;
use crate::rdkit_torsion::{CODE_SIZE, atom_code};

/// `AtomPairArguments`: `minDistance = 1`, `maxDistance = maxPathLen - 1 = 30`
/// (`numPathBits = 5` => `maxPathLen = (1 << 5) - 1 = 31`). Pairs outside this
/// range are skipped entirely, not clamped.
/// `numPathBits`.
const NUM_PATH_BITS: u32 = 5;
const MIN_DISTANCE: u32 = 1;
const MAX_DISTANCE: u32 = 30;

const COUNT_BOUNDS: [u32; 4] = [1, 2, 4, 8];
const N_BITS_PER_ENTRY: usize = 4;
const MAX_CONFIGURED_BITS: usize = 1_048_576;

/// Visit each undirected pair once. The queue also records which distances
/// to clear before the next source, so distant and disconnected atoms cost
/// no work per source.
fn visit_pairs(
    mol: &Molecule,
    distances: RangeInclusive<u32>,
    mut visit: impl FnMut(usize, usize, u32),
) {
    const UNVISITED: u8 = u8::MAX;
    let mut distance = vec![UNVISITED; mol.atom_count()];
    let mut queue = Vec::new();
    let min = *distances.start();
    let max = *distances.end();
    if max == 0 {
        return;
    }
    for start in 0..mol.atom_count() {
        distance[start] = 0;
        queue.push(AtomIdx(start as u32));
        let mut head = 0;
        while let Some(&current) = queue.get(head) {
            head += 1;
            let next = u32::from(distance[current.0 as usize]) + 1;
            if next > max {
                continue;
            }
            for (neighbor, _) in mol.neighbors(current) {
                let end = neighbor.0 as usize;
                if distance[end] != UNVISITED {
                    continue;
                }
                distance[end] = next as u8;
                queue.push(neighbor);
                if end > start && next >= min {
                    visit(start, end, next);
                }
            }
        }
        for atom in queue.drain(..) {
            distance[atom.0 as usize] = UNVISITED;
        }
    }
}

fn atom_codes(mol: &Molecule, chirality: &[Option<u8>], code_size: u32) -> Vec<u32> {
    let limit = (1 << code_size) - 1;
    (0..mol.atom_count())
        .map(|idx| {
            let stereo_bits = match chirality.get(idx).copied().flatten() {
                Some(b'R') => 1 << CODE_SIZE,
                Some(b'S') => 2 << CODE_SIZE,
                _ => 0,
            };
            (atom_code(mol, AtomIdx(idx as u32), 0) | stereo_bits) % limit
        })
        .collect()
}

fn pair_codes(codes: &[u32], i: usize, j: usize) -> (u32, u32) {
    (codes[i].min(codes[j]), codes[i].max(codes[j]))
}

fn sparse_code(lo: u32, hi: u32, distance: u32, code_size: u32) -> u32 {
    distance | (lo << NUM_PATH_BITS) | (hi << (NUM_PATH_BITS + code_size))
}

fn sparse_elements(counts: FxHashMap<u32, u32>) -> Vec<(u32, u32)> {
    let mut out: Vec<_> = counts.into_iter().collect();
    out.sort_unstable_by_key(|&(code, _)| code);
    out
}

fn folded_elements(counts: Vec<u32>) -> Vec<(u32, u32)> {
    counts
        .into_iter()
        .enumerate()
        .filter(|&(_, count)| count > 0)
        .map(|(bucket, count)| (bucket as u32, count))
        .collect()
}

/// RDKit-compatible Atom Pair fingerprint, matching
/// `rdkit.Chem.AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)`
/// (`includeChirality=False`, `nBitsPerEntry=4`, the Python API's own defaults)
/// on the validated corpora. See this module's corpus results and the shared
/// atom-code implementation at [`crate::rdkit_torsion::num_pi_electrons`].
pub fn rdkit_atom_pair_fp(mol: &Molecule) -> BitVec2048 {
    chematic_smiles::with_rdkit_model_molecule(mol, rdkit_atom_pair_fp_impl)
}

fn rdkit_atom_pair_fp_impl(mol: &Molecule) -> BitVec2048 {
    // RDKit fingerprints sanitized molecules, i.e. after aromaticity
    // perception; Kekule input must not see a different graph than the
    // equivalent aromatic spelling. The perceived view is memoized on `mol`.
    chematic_perception::with_rdkit_parity_view(mol, |view| {
        rdkit_atom_pair_fp_prepared(view.unwrap_or(mol))
    })
}

/// RDKit's hashed atom-pair *count* fingerprint
/// (`rdMolDescriptors.GetHashedAtomPairFingerprint(mol, nBits)`): the
/// nonzero `(bucket, count)` elements, sorted by bucket.
pub fn rdkit_atom_pair_counts(mol: &Molecule, n_bits: u32) -> Vec<(u32, u32)> {
    chematic_smiles::with_rdkit_model_molecule(mol, |m| rdkit_atom_pair_counts_impl(m, n_bits))
}

fn rdkit_atom_pair_counts_impl(mol: &Molecule, n_bits: u32) -> Vec<(u32, u32)> {
    let n_bits = n_bits.max(1);
    chematic_perception::with_rdkit_parity_view(mol, |view| {
        let counts = rdkit_atom_pair_bucket_counts(view.unwrap_or(mol), n_bits);
        folded_elements(counts)
    })
}

/// RDKit's unhashed atom-pair count fingerprint
/// (`rdFingerprintGenerator.GetAtomPairGenerator()
/// .GetSparseCountFingerprint(mol)`, distances 1..=30, no chirality):
/// `getAtomPairCode(min code, max code, distance)` per pair, the nonzero
/// `(code, count)` elements sorted by code.
pub fn rdkit_atom_pair_sparse_counts(mol: &Molecule) -> Vec<(u32, u32)> {
    chematic_perception::with_rdkit_parity_view(mol, |view| {
        let mol = view.unwrap_or(mol);
        let codes = atom_codes(mol, &[], CODE_SIZE);
        let mut counts = FxHashMap::<u32, u32>::default();
        visit_pairs(mol, MIN_DISTANCE..=MAX_DISTANCE, |i, j, distance| {
            let (lo, hi) = pair_codes(&codes, i, j);
            let count = counts
                .entry(sparse_code(lo, hi, distance, CODE_SIZE))
                .or_default();
            *count = count.saturating_add(1);
        });
        sparse_elements(counts)
    })
}

/// Per-bucket pair counts of `getHashedAtomPairFingerprint(mol, n_buckets)`
/// on a molecule whose aromaticity is already RDKit-perceived.
fn rdkit_atom_pair_bucket_counts(mol: &Molecule, n_buckets: u32) -> Vec<u32> {
    let codes = atom_codes(mol, &[], CODE_SIZE);
    let mut counts = vec![0u32; n_buckets as usize];
    visit_pairs(mol, MIN_DISTANCE..=MAX_DISTANCE, |i, j, distance| {
        let (lo, hi) = pair_codes(&codes, i, j);
        let bucket = (hash_vec(&[lo, distance, hi]) % n_buckets) as usize;
        counts[bucket] = counts[bucket].saturating_add(1);
    });
    counts
}

/// [`rdkit_atom_pair_fp`] on a molecule whose aromaticity is already RDKit-perceived.
fn rdkit_atom_pair_fp_prepared(mol: &Molecule) -> BitVec2048 {
    const BLOCK_LENGTH: u32 = 2048 / N_BITS_PER_ENTRY as u32;
    let counts = rdkit_atom_pair_bucket_counts(mol, BLOCK_LENGTH);
    let mut fp = BitVec2048::new();
    for (bucket, &count) in counts.iter().enumerate() {
        for (i, &bound) in COUNT_BOUNDS.iter().enumerate() {
            if count >= bound {
                fp.set(bucket * N_BITS_PER_ENTRY + i);
            }
        }
    }
    fp
}

/// Options for the two-dimensional RDKit Atom Pair generator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitAtomPairConfig {
    /// Inclusive topological distance limits, between 0 and 30.
    pub min_distance: u32,
    pub max_distance: u32,
    /// Encode legacy RDKit R/S atom labels in the atom invariant.
    pub include_chirality: bool,
    /// Width of the packed and folded count fingerprints.
    pub nbits: usize,
    /// Thresholds in bit order. `None` disables count simulation.
    pub count_bounds: Option<Vec<u32>>,
}

impl Default for RdkitAtomPairConfig {
    fn default() -> Self {
        Self {
            min_distance: MIN_DISTANCE,
            max_distance: MAX_DISTANCE,
            include_chirality: false,
            nbits: 2048,
            count_bounds: Some(COUNT_BOUNDS.to_vec()),
        }
    }
}

/// The packed bits and nonzero counts from one Atom Pair calculation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitAtomPairResult {
    pub fingerprint: BitVecN,
    /// Hashed counts folded to `nbits`, independent of count simulation.
    pub folded_counts: Vec<(u32, u32)>,
    /// Unhashed counts, independent of `nbits` and count simulation.
    pub sparse_counts: Vec<(u32, u32)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitAtomPairError {
    InvalidDistances,
    InvalidBitWidth,
    InvalidCountBounds,
    Chirality(chematic_smiles::RdkitSmilesError),
}

impl fmt::Display for RdkitAtomPairError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDistances => f.write_str("Atom Pair requires 0 <= min_distance <= max_distance <= 30"),
            Self::InvalidBitWidth => write!(f, "Atom Pair bit width must be between 1 and {MAX_CONFIGURED_BITS}"),
            Self::InvalidCountBounds => f.write_str("Atom Pair count simulation needs nonempty bounds shorter than nbits and at most 2147483647"),
            Self::Chirality(error) => write!(f, "Atom Pair chirality: {error}"),
        }
    }
}

impl std::error::Error for RdkitAtomPairError {}

/// Compute the configured 2D Atom Pair outputs on the RDKit input model.
/// A missing legacy stereo calculation is reported when chirality is enabled.
pub fn rdkit_atom_pair_config(
    mol: &Molecule,
    config: &RdkitAtomPairConfig,
) -> Result<RdkitAtomPairResult, RdkitAtomPairError> {
    if config.min_distance > config.max_distance || config.max_distance > MAX_DISTANCE {
        return Err(RdkitAtomPairError::InvalidDistances);
    }
    if !(1..=MAX_CONFIGURED_BITS).contains(&config.nbits) {
        return Err(RdkitAtomPairError::InvalidBitWidth);
    }
    if config.count_bounds.as_ref().is_some_and(|bounds| {
        bounds.is_empty()
            || bounds.len() >= config.nbits
            || bounds.iter().any(|&bound| bound > i32::MAX as u32)
    }) {
        return Err(RdkitAtomPairError::InvalidCountBounds);
    }
    chematic_smiles::with_rdkit_model_molecule(mol, |model| {
        let chirality = if config.include_chirality
            && model
                .atoms()
                .any(|(_, atom)| atom.chirality.is_tetrahedral())
        {
            chematic_smiles::rdkit_legacy_stereo(model)
                .map_err(RdkitAtomPairError::Chirality)?
                .atom_cip
        } else {
            Vec::new()
        };
        Ok(chematic_perception::with_rdkit_parity_view(model, |view| {
            configured_pairs(view.unwrap_or(model), config, &chirality)
        }))
    })
}

fn configured_pairs(
    mol: &Molecule,
    config: &RdkitAtomPairConfig,
    chirality: &[Option<u8>],
) -> RdkitAtomPairResult {
    let code_size = CODE_SIZE + if config.include_chirality { 2 } else { 0 };
    let codes = atom_codes(mol, chirality, code_size);
    let bounds = config.count_bounds.as_deref().unwrap_or(&[1]);
    let mut folded = vec![0u32; config.nbits];
    let sim_width = config.nbits / bounds.len();
    let mut simulated = (sim_width != config.nbits).then(|| vec![0u32; sim_width]);
    let mut sparse = FxHashMap::<u32, u32>::default();
    visit_pairs(
        mol,
        config.min_distance..=config.max_distance,
        |i, j, distance| {
            let (lo, hi) = pair_codes(&codes, i, j);
            let hash = hash_vec(&[lo, distance, hi]) as usize;
            let count = &mut folded[hash % config.nbits];
            *count = count.saturating_add(1);
            if let Some(simulated) = &mut simulated {
                let count = &mut simulated[hash % sim_width];
                *count = count.saturating_add(1);
            }
            let count = sparse
                .entry(sparse_code(lo, hi, distance, code_size))
                .or_default();
            *count = count.saturating_add(1);
        },
    );
    let mut fingerprint = BitVecN::new(config.nbits);
    let bit_counts = simulated.as_deref().unwrap_or(&folded);
    for (bucket, &count) in bit_counts
        .iter()
        .enumerate()
        .filter(|&(_, &count)| count > 0)
    {
        for (i, &bound) in bounds.iter().enumerate() {
            if count >= bound {
                fingerprint.set(bucket * bounds.len() + i);
            }
        }
    }
    RdkitAtomPairResult {
        fingerprint,
        folded_counts: folded_elements(folded),
        sparse_counts: sparse_elements(sparse),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn deterministic() {
        let m = mol("CCO");
        assert_eq!(rdkit_atom_pair_fp(&m), rdkit_atom_pair_fp(&m));
    }

    #[test]
    fn different_molecules_differ() {
        assert_ne!(
            rdkit_atom_pair_fp(&mol("CCO")),
            rdkit_atom_pair_fp(&mol("CCN"))
        );
    }

    #[test]
    fn single_atom_has_no_pairs() {
        assert_eq!(rdkit_atom_pair_fp(&mol("C")), BitVec2048::new());
    }

    #[test]
    fn configured_chiral_pairs_match_rdkit_in_a_non_power_of_two_width() {
        let config = RdkitAtomPairConfig {
            min_distance: 2,
            max_distance: 3,
            include_chirality: true,
            nbits: 17,
            count_bounds: Some(vec![3, 1, 2]),
        };
        // RDKit 2026.09.2 GetAtomPairGenerator, with a disconnected chloride.
        let left = rdkit_atom_pair_config(&mol("N[C@@H](C)C(=O)O.Cl"), &config).unwrap();
        let right = rdkit_atom_pair_config(&mol("N[C@H](C)C(=O)O.Cl"), &config).unwrap();
        assert_eq!(left.fingerprint.to_le_bytes(), [146, 125, 0]);
        assert_eq!(right.fingerprint.to_le_bytes(), [146, 126, 0]);
        assert_eq!(
            left.folded_counts,
            [(0, 1), (2, 1), (3, 2), (7, 2), (8, 1), (11, 3)]
        );
        assert_eq!(
            left.sparse_counts
                .iter()
                .map(|&(id, _)| id)
                .collect::<Vec<_>>(),
            [
                2819106, 4260898, 4261218, 6358051, 6359075, 6882339, 6883363, 6884386, 69405730,
                69405986
            ]
        );
        assert!(left.sparse_counts.iter().all(|&(_, count)| count == 1));
    }

    #[test]
    fn configured_empty_distance_range_and_invalid_options() {
        let molecule = mol("C1CC1.F");
        let mut config = RdkitAtomPairConfig {
            min_distance: 0,
            max_distance: 0,
            nbits: 1,
            count_bounds: None,
            ..Default::default()
        };
        let result = rdkit_atom_pair_config(&molecule, &config).unwrap();
        assert_eq!(result.fingerprint.to_le_bytes(), [0]);
        assert!(result.sparse_counts.is_empty());
        assert!(result.folded_counts.is_empty());

        config.min_distance = 1;
        assert_eq!(
            rdkit_atom_pair_config(&molecule, &config),
            Err(RdkitAtomPairError::InvalidDistances)
        );
        config.max_distance = 1;
        config.nbits = 0;
        assert_eq!(
            rdkit_atom_pair_config(&molecule, &config),
            Err(RdkitAtomPairError::InvalidBitWidth)
        );
        config.nbits = 1;
        config.count_bounds = Some(vec![1]);
        assert_eq!(
            rdkit_atom_pair_config(&molecule, &config),
            Err(RdkitAtomPairError::InvalidCountBounds)
        );
    }

    #[test]
    fn visits_pairs_beyond_the_native_atom_pair_cap() {
        // A 10-carbon chain has a pair at distance 9, well past the native
        // `atom_pair.rs` MAX_DIST=7 cap -- must not be silently dropped here.
        let m = mol("CCCCCCCCCC");
        let mut pairs = Vec::new();
        visit_pairs(&m, 8..=9, |i, j, distance| pairs.push((i, j, distance)));
        assert_eq!(pairs, [(0, 8, 8), (0, 9, 9), (1, 9, 8)]);
    }
}
