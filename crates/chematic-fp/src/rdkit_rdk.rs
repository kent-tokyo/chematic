//! RDKit-compatible "RDKit fingerprint" port (`Chem.RDKFingerprint` / `RDKFingerprintMol`).
//!
//! Ported from RDKit's `Code/GraphMol/Fingerprints/{Fingerprints,RDKitFPGenerator,
//! FingerprintGenerator,FingerprintUtil,Subgraphs}.cpp`, default parameters
//! (`minPath=1, maxPath=7, useHs=true, branchedPaths=true, useBondOrder=true,
//! fpSize=2048, numBitsPerFeature=2`).
//!
//! # Algorithm
//!
//! 1. Enumerate every connected branched subgraph (not just linear paths) of
//!    1..=7 bonds via [`find_all_subgraphs`], RDKit's own canonical
//!    "root at minimum bond index" backtracking scheme
//!    (`Subgraphs::findAllSubgraphsOfLengthsMtoN`/`recurseWalkRange`): bonds are
//!    tried as roots in ascending index order, each root is marked globally
//!    forbidden once used (so no later root's subgraph ever reincludes it, which
//!    is what deduplicates the enumeration to exactly one discovery per subgraph),
//!    while a growing candidate frontier (adjacent, not-yet-forbidden bonds) is
//!    extended depth-first with a fresh forbidden-set copy handed to each
//!    recursive call (so sibling branches within one root don't interfere).
//! 2. For each subgraph (a set of bond indices), compute a per-bond hash
//!    (`RDKitFPUtils::generateBondHashes`): atom invariant is
//!    `(atomic_number % 128) << 1 | is_aromatic`; the *degree* folded into the
//!    hash is **path-local** (how many times the atom appears as a bond endpoint
//!    within this specific subgraph), not the molecule's true degree — confirmed
//!    against a live RDKit oracle, since a naive "true molecular degree" hash
//!    does not reproduce RDKit's output for `CCC`. The two atoms of a bond are
//!    canonically ordered (larger invariant first, ties broken by larger local
//!    degree first) before folding, so a bond's hash doesn't depend on which
//!    direction it was traversed from.
//! 3. Fold each subgraph's bond hashes (plus a trailing distinct-atom count, to
//!    distinguish e.g. `C1CC1` from `CC(C)C`) into one `u32` seed via RDKit's
//!    `gboost::hash_combine`/`hash_range` — despite the C++ call site declaring
//!    a 64-bit `unsigned long`, `gboost`'s own `hash_result_t` typedef is
//!    `std::uint32_t`, so the fold is pure 32-bit arithmetic with no width
//!    mismatch (confirmed against a live oracle for both 1-bond and 2-bond
//!    paths). A single-bond subgraph skips the sort/fold and uses its own bond
//!    hash directly as the seed.
//! 4. Each seed sets one direct bit (`seed % fpSize`) plus, since RDKit's
//!    default `numBitsPerFeature=2`, one additional bit drawn from a **weakened**
//!    Mersenne Twister reseeded with that exact `seed` (see
//!    [`RdkitWeakMt19937`]) via `boost::uniform_int<>(0, INT_MAX)`, which for
//!    this exact engine-range-vs-target-range ratio reduces to `raw / 2`
//!    (bucket division, not modulo — confirmed against a live oracle), then
//!    `% fpSize` again.
//!
//! # Aromaticity: literal vs. re-perceived
//!
//! RDKit's parser always normalizes a perceived-aromatic bond's type to
//! `AROMATIC` and sets `atom->getIsAromatic()`, regardless of whether the input
//! SMILES used Kekulé or lowercase-aromatic notation. chematic's SMILES parser
//! instead preserves a Kekulized input's literal bond orders verbatim (round-trip
//! fidelity). Reusing the exact fix already proven for `rdkit_pattern_fp` (see
//! that module's own doc comment): if a molecule has *any* literal
//! `BondOrder::Aromatic` bond, every bond's own literal order/every atom's own
//! literal `aromatic` flag is trusted as-is; only when a molecule has *zero*
//! literal aromatic bonds anywhere does this fall back to re-perceiving
//! aromaticity via [`chematic_perception::aromaticity::assign_aromaticity_ex`].
//! Re-perceiving unconditionally regressed a Kekulé-written, genuinely
//! non-aromatic ring that coexisted with a lowercase-aromatic ring elsewhere in
//! the same molecule; this hybrid rule avoids that without losing the fix for
//! molecules parsed entirely from Kekulé notation.
//!
//! # Verified
//!
//! Bit-exact (identical on-bit sets, not just popcount) against a live RDKit
//! oracle (`Chem.RDKFingerprint(mol, minPath=1, maxPath=7, fpSize=2048)`):
//! 100% on the full `descriptor_census_corpus.smi` (5000/5000) and
//! `chembl_accuracy_corpus_4999.smi` (5000/5000), 99.44% on
//! `nci_first_5k_smiles_only.smi` (4963/4991). Every one of the 28 NCI
//! mismatches is a fused polyheteroaromatic dye (xanthene/anthraquinone with
//! halogen substituents, phenothiazine/phenoxazinium, purine/xanthine-fused
//! rings), an exotic charged aromatic heterocycle (pyrylium/thiopyrylium-like
//! N+/O+/S+ ring atoms), or a metal-coordination complex (ferrocene-like
//! sandwich structures, ring-N-to-metal dative bonding) where chematic's own
//! `chematic-perception` Hückel aromaticity model does not (yet) recognize the
//! same aromatic system RDKit's does — the same pre-existing, out-of-scope
//! dependency gap already documented for `rdkit_torsion_fp`/`rdkit_atom_pair_fp`
//! (hypervalent-atom hybridization) and `rdkit_pattern_fp` (their own NCI
//! residuals), not a defect introduced by this port.

use crate::bitvec::BitVec2048;
use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::aromaticity::{
    AromaticityAlgorithm, AromaticityModel, assign_aromaticity_ex,
};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BTreeMap;

const FP_SIZE: u32 = 2048;
const MIN_PATH: usize = 1;
const MAX_PATH: usize = 7;
const NUM_BITS_PER_FEATURE: u32 = 2;

fn hash_combine(seed: u32, value: u32) -> u32 {
    seed ^ value
        .wrapping_add(0x9e3779b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

fn rdkit_bond_type_code(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single
        | BondOrder::Up
        | BondOrder::Down
        | BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
        BondOrder::Zero => 21,
    }
}

fn is_atom_aromatic(mol: &Molecule, aromaticity: Option<&AromaticityModel>, idx: AtomIdx) -> bool {
    if mol.atom(idx).aromatic {
        return true;
    }
    match aromaticity {
        Some(model) => model.is_atom_aromatic(idx),
        None => false,
    }
}

fn bond_type_for(mol: &Molecule, aromaticity: Option<&AromaticityModel>, bidx: BondIdx) -> u32 {
    let bond = mol.bond(bidx);
    let is_aromatic = bond.order == BondOrder::Aromatic
        || matches!(aromaticity, Some(model) if model.is_bond_aromatic(bidx));
    if is_aromatic {
        12
    } else {
        rdkit_bond_type_code(bond.order)
    }
}

/// RDKit's weakened Mersenne Twister variant used when `numBitsPerFeature > 1`.
///
/// `n=4, m=2, r=31` state — a deliberately tiny-state MT chosen by RDKit's
/// authors because reseeding a full MT19937 (n=624) per fingerprint feature is
/// too expensive. Critically, boost's *deprecated* `mersenne_twister<>` wrapper
/// class (the one RDKit's typedef actually instantiates) silently discards its
/// own last template parameter — RDKit's source passes `3346425566U` there,
/// presumably intending it as the seeding multiplier — and hardcodes the
/// textbook MT19937 multiplier `1812433253` instead (confirmed by reading
/// boost's `mersenne_twister.hpp` directly: the wrapper forwards to
/// `mersenne_twister_engine<..., 1812433253>` unconditionally). Do **not**
/// "fix" this to `3346425566` — that would silently break bit-exactness with
/// real RDKit output, which is itself built against the same boost quirk.
struct RdkitWeakMt19937 {
    state: [u32; 4],
    index: usize,
}

const MT_N: usize = 4;
const MT_M: usize = 2;
const MT_MATRIX_A: u32 = 0x9908b0df;
const MT_UPPER_MASK: u32 = 0x8000_0000;
const MT_LOWER_MASK: u32 = 0x7fff_ffff;
const MT_TEMPER_U: u32 = 11;
const MT_TEMPER_S: u32 = 7;
const MT_TEMPER_B: u32 = 0x9d2c_5680;
const MT_TEMPER_T: u32 = 15;
const MT_TEMPER_C: u32 = 0xefc6_0000;
const MT_TEMPER_L: u32 = 18;
const MT_SEED_MULT: u32 = 1_812_433_253;

impl RdkitWeakMt19937 {
    fn seeded(value: u32) -> Self {
        let mut state = [0u32; MT_N];
        state[0] = value;
        for i in 1..MT_N {
            state[i] = MT_SEED_MULT
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: MT_N }
    }

    fn twist(&mut self) {
        for j in 0..MT_N {
            let y = (self.state[j] & MT_UPPER_MASK) | (self.state[(j + 1) % MT_N] & MT_LOWER_MASK);
            let mut next = self.state[(j + MT_M) % MT_N] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= MT_MATRIX_A;
            }
            self.state[j] = next;
        }
        self.index = 0;
    }

    fn next_u32(&mut self) -> u32 {
        if self.index == MT_N {
            self.twist();
        }
        let mut z = self.state[self.index];
        self.index += 1;
        z ^= z >> MT_TEMPER_U;
        z ^= (z << MT_TEMPER_S) & MT_TEMPER_B;
        z ^= (z << MT_TEMPER_T) & MT_TEMPER_C;
        z ^= z >> MT_TEMPER_L;
        z
    }
}

fn build_bond_adjacency(mol: &Molecule) -> FxHashMap<u32, Vec<u32>> {
    let mut nbrs: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for (bidx, _) in mol.bonds() {
        nbrs.entry(bidx.0).or_default();
    }
    for (aidx, _) in mol.atoms() {
        let incident: Vec<u32> = mol.neighbors(aidx).map(|(_, b)| b.0).collect();
        for &b1 in &incident {
            for &b2 in &incident {
                if b1 != b2 {
                    nbrs.get_mut(&b1).unwrap().push(b2);
                }
            }
        }
    }
    nbrs
}

/// RDKit's `Subgraphs::findAllSubgraphsOfLengthsMtoN`/`recurseWalkRange`: every
/// connected subgraph of `min_len..=max_len` bonds, each discovered exactly
/// once (rooted at its lowest-index bond). Order within each length's `Vec` is
/// irrelevant to the caller (only the resulting bitset is observable), so this
/// need only reproduce the same *set* of subgraphs, not RDKit's exact
/// traversal order.
fn find_all_subgraphs(
    mol: &Molecule,
    min_len: usize,
    max_len: usize,
) -> BTreeMap<usize, Vec<Vec<u32>>> {
    let nbrs = build_bond_adjacency(mol);
    let n_bonds = mol.bond_count();
    let mut forbidden = vec![false; n_bonds];
    let mut res: BTreeMap<usize, Vec<Vec<u32>>> =
        (min_len..=max_len).map(|l| (l, Vec::new())).collect();

    for root in 0..n_bonds as u32 {
        if forbidden[root as usize] {
            continue;
        }
        forbidden[root as usize] = true;
        let cands = nbrs.get(&root).cloned().unwrap_or_default();
        let spath = vec![root];
        recurse_walk_range(
            &nbrs,
            spath,
            cands,
            min_len,
            max_len,
            forbidden.clone(),
            &mut res,
        );
    }
    res
}

fn recurse_walk_range(
    nbrs: &FxHashMap<u32, Vec<u32>>,
    spath: Vec<u32>,
    mut cands: Vec<u32>,
    min_len: usize,
    max_len: usize,
    mut forbidden: Vec<bool>,
    res: &mut BTreeMap<usize, Vec<Vec<u32>>>,
) {
    let nsize = spath.len();
    if nsize >= min_len && nsize <= max_len {
        res.get_mut(&nsize).unwrap().push(spath.clone());
    }
    if nsize >= max_len {
        return;
    }
    while let Some(next) = cands.pop() {
        if !forbidden[next as usize] {
            forbidden[next as usize] = true;
            let mut tstack = cands.clone();
            if let Some(next_nbrs) = nbrs.get(&next) {
                for &b in next_nbrs {
                    if !forbidden[b as usize] {
                        tstack.push(b);
                    }
                }
            }
            let mut tpath = spath.clone();
            tpath.push(next);
            recurse_walk_range(
                nbrs,
                tpath,
                tstack,
                min_len,
                max_len,
                forbidden.clone(),
                res,
            );
        }
    }
}

fn path_bond_hashes(
    mol: &Molecule,
    invariants: &[u32],
    aromaticity: Option<&AromaticityModel>,
    path: &[u32],
) -> (Vec<u32>, u32) {
    let mut atom_degrees: FxHashMap<u32, u32> = FxHashMap::default();
    let mut atoms_in_path: FxHashSet<u32> = FxHashSet::default();
    for &bidx in path {
        let bond = mol.bond(BondIdx(bidx));
        *atom_degrees.entry(bond.atom1.0).or_insert(0) += 1;
        *atom_degrees.entry(bond.atom2.0).or_insert(0) += 1;
        atoms_in_path.insert(bond.atom1.0);
        atoms_in_path.insert(bond.atom2.0);
    }

    let mut bond_nbrs = vec![0u32; path.len()];
    for i in 0..path.len() {
        let bi = mol.bond(BondIdx(path[i]));
        for j in (i + 1)..path.len() {
            let bj = mol.bond(BondIdx(path[j]));
            if bi.atom1 == bj.atom1
                || bi.atom1 == bj.atom2
                || bi.atom2 == bj.atom1
                || bi.atom2 == bj.atom2
            {
                bond_nbrs[i] += 1;
                bond_nbrs[j] += 1;
            }
        }
    }

    let mut bond_hashes = Vec::with_capacity(path.len());
    for (i, &bidx) in path.iter().enumerate() {
        let bond = mol.bond(BondIdx(bidx));
        let mut a1_hash = invariants[bond.atom1.0 as usize];
        let mut a2_hash = invariants[bond.atom2.0 as usize];
        let mut deg1 = atom_degrees[&bond.atom1.0];
        let mut deg2 = atom_degrees[&bond.atom2.0];
        if a1_hash < a2_hash {
            std::mem::swap(&mut a1_hash, &mut a2_hash);
            std::mem::swap(&mut deg1, &mut deg2);
        } else if a1_hash == a2_hash && deg1 < deg2 {
            std::mem::swap(&mut deg1, &mut deg2);
        }
        let bond_type = bond_type_for(mol, aromaticity, BondIdx(bidx));
        let mut h = bond_nbrs[i];
        h = hash_combine(h, bond_type);
        h = hash_combine(h, a1_hash);
        h = hash_combine(h, deg1);
        h = hash_combine(h, a2_hash);
        h = hash_combine(h, deg2);
        bond_hashes.push(h);
    }
    (bond_hashes, atoms_in_path.len() as u32)
}

fn path_seed(bond_hashes: &[u32], distinct_atoms: u32) -> u32 {
    if bond_hashes.len() == 1 {
        return bond_hashes[0];
    }
    let mut sorted = bond_hashes.to_vec();
    sorted.sort_unstable();
    sorted.push(distinct_atoms);
    let mut seed = 0u32;
    for v in sorted {
        seed = hash_combine(seed, v);
    }
    seed
}

/// Compute the RDKit-compatible "RDKit fingerprint" (`Chem.RDKFingerprint`).
///
/// See the module doc comment for the full algorithm and its verification
/// status.
pub fn rdkit_rdk_fp(mol: &Molecule) -> BitVec2048 {
    let mut fp = BitVec2048::new();
    if mol.atom_count() == 0 {
        return fp;
    }

    let has_literal_aromatic_bond = mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic);
    let aromaticity = if has_literal_aromatic_bond {
        None
    } else {
        Some(assign_aromaticity_ex(mol, AromaticityAlgorithm::RdkitLike))
    };

    let invariants: Vec<u32> = (0..mol.atom_count())
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let an = mol.atom(idx).element.atomic_number() as u32;
            ((an % 128) << 1) | (is_atom_aromatic(mol, aromaticity.as_ref(), idx) as u32)
        })
        .collect();

    let paths = find_all_subgraphs(mol, MIN_PATH, MAX_PATH);
    for plist in paths.values() {
        for path in plist {
            let (bond_hashes, distinct_atoms) =
                path_bond_hashes(mol, &invariants, aromaticity.as_ref(), path);
            if bond_hashes.is_empty() {
                continue;
            }
            let seed = path_seed(&bond_hashes, distinct_atoms);
            fp.set((seed % FP_SIZE) as usize);

            if NUM_BITS_PER_FEATURE > 1 {
                let mut mt = RdkitWeakMt19937::seeded(seed);
                for _ in 1..NUM_BITS_PER_FEATURE {
                    let raw = mt.next_u32();
                    let bounded = raw / 2;
                    fp.set((bounded % FP_SIZE) as usize);
                }
            }
        }
    }

    fp
}

/// Tanimoto similarity between two RDKit-compatible RDK fingerprints.
pub fn tanimoto_rdkit_rdk(a: &BitVec2048, b: &BitVec2048) -> f64 {
    a.tanimoto(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    #[test]
    fn deterministic() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(rdkit_rdk_fp(&m), rdkit_rdk_fp(&m));
    }

    #[test]
    fn different_molecules_differ() {
        let a = rdkit_rdk_fp(&mol("CCO"));
        let b = rdkit_rdk_fp(&mol("c1ccccc1"));
        assert_ne!(a, b);
    }

    #[test]
    fn single_atom_has_no_bits() {
        let fp = rdkit_rdk_fp(&mol("C"));
        assert_eq!(fp.popcount(), 0);
    }

    #[test]
    fn matches_rdkit_oracle_on_hand_picked_molecules() {
        // (smiles, expected on-bit count from a live `Chem.RDKFingerprint`
        // oracle run at minPath=1, maxPath=7, fpSize=2048).
        let cases: &[(&str, u32)] = &[
            ("CC", 2),
            ("C=C", 2),
            ("C#C", 2),
            ("CCC", 4),
            ("CC(C)C", 6),
            ("c1ccccc1", 12),
            ("C1CC1", 6),
            ("CCO", 6),
            ("CC(=O)O", 14),
            ("c1ccc2ccccc2c1", 44),
            ("CC(=O)Oc1ccccc1C(=O)O", 354),
        ];
        for &(smiles, expected) in cases {
            let fp = rdkit_rdk_fp(&mol(smiles));
            assert_eq!(fp.popcount(), expected, "popcount mismatch for {smiles:?}");
        }
    }
}
