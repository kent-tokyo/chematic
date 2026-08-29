//! RDKit-bit-exact Pattern fingerprint (`Chem.PatternFingerprint`).
//!
//! Opt-in `RdkitExact`-mode addition, kept fully separate from
//! [`crate::pattern::pattern_fp`] (chematic's own native scheme). Third of the
//! fingerprint-parity series (Track A / 99-point directive Phase 6), after
//! [`crate::rdkit_torsion_fp`]/[`crate::rdkit_atom_pair_fp`]. Structurally
//! unrelated to those two -- Pattern uses SMARTS substructure matching against
//! 13 fixed patterns, not a path/pair enumeration + shared atom invariant.
//!
//! Derived from RDKit's C++ source
//! (`Code/GraphMol/Fingerprints/PatternFingerprints.cpp`, function
//! `updatePatternFingerprint` / `PatternFingerprintMol`), fetched at research
//! time against `master`. Every constant/formula below is a direct port.
//!
//! Algorithm, confirmed against source before implementing:
//! - [`PATTERNS`] is RDKit's own `pqs[]` array verbatim (13 SMARTS strings,
//!   1-based index `pIdx` used as a hash seed component).
//! - For each pattern, find every match via chematic's own
//!   `chematic_smarts::find_matches_with_config` with `uniquify: false` --
//!   this already reproduces RDKit's own `SubstructMatch(..., uniquify=false)`
//!   semantics exactly, *including* returning one entry per automorphism
//!   (e.g. a 6-membered-ring pattern on benzene yields 12 matches, not 1) --
//!   verified directly against a live RDKit oracle on 9 pattern/molecule pairs
//!   before writing this module, not assumed from the API's docstring.
//! - **Occurrence-count bits**: RDKit's C++ mutates a single seed variable
//!   `mIdx` (initialized to `pIdx + numPatternAtoms + numPatternBonds`) via
//!   one `hash_combine(mIdx, 0xBEEF)` step *per match*, carrying the updated
//!   value forward to the next match of the *same* pattern (not reset per
//!   match) -- a running accumulator, not a fresh hash per match. Since the
//!   combined value (`0xBEEF`) never depends on which atoms matched, this
//!   produces a sequence of `k` distinct bits purely as a function of the
//!   pattern's own seed and the match *count* `k` -- independent of match
//!   *order*, so this module doesn't need to replicate RDKit's own match
//!   enumeration order for this part (confirmed by construction, not by
//!   trial and error).
//! - **Content bits**: `bitId` starts fresh at `pIdx` per match, then folds in
//!   (via [`hash_combine`]) each matched atom's atomic number in *query-atom-
//!   index* order (0..n), followed by each pattern bond's *matched* bond type
//!   in the pattern's own bond order (query bonds and atoms are `QueryMolecule`
//!   `Vec`s, i.e. already in SMARTS-text parse order -- the same order RDKit's
//!   own query mol exposes). Aromatic bonds always hash RDKit's `AROMATIC`
//!   code (`12`) regardless of Kekule order -- determined via
//!   [`chematic_perception::aromaticity::assign_aromaticity_ex`]
//!   (`AromaticityAlgorithm::RdkitLike`), **not** a raw `bond.order ==
//!   BondOrder::Aromatic` check: chematic's parser preserves a Kekulized
//!   *input*'s literal single/double bond orders verbatim even when the ring
//!   is perceived as aromatic, while RDKit always normalizes such bonds to
//!   `AROMATIC` type regardless of input notation. An earlier version of this
//!   module used the raw-order check and was 100% bit-exact on two corpora
//!   but only 30.4% on a third (NCI) before this was caught -- every
//!   Kekule-notation heteroaromatic (e.g. `S1C2=CC3=CC=CC=C3C=C2N=C1...`)
//!   silently hashed the wrong bond code. See [`matched_bond_code`]'s own doc
//!   comment for the confirmed repro.
//! - `fpSize` buckets directly (`bitId % fpSize`) -- **no** count-simulation
//!   tail (unlike torsion/atom-pair): `PatternFingerprintMol`'s signature has
//!   no `nBitsPerEntry` parameter at all, confirmed from source before
//!   assuming the torsion/atom-pair tail would apply here too.
//!
//! **Not implemented** (mirrors this series' established chirality/query
//! carve-outs): `tautomericFingerprint=True` (a second, alternate bond-hash
//! path for tautomer-invariant fingerprinting), and the `isQueryAtom`/
//! `isQueryBond` skip logic (only reachable when fingerprinting a
//! `QueryMolecule`/SMARTS pattern itself -- this module, like
//! `rdkit_torsion_fp`/`rdkit_atom_pair_fp`, targets fingerprinting concrete
//! molecules parsed from SMILES, which never carry query atoms/bonds, so that
//! branch is unreachable dead code for this scope, not a gap).
//!
//! **Status, measured against a live RDKit oracle (1000-molecule samples,
//! three corpora with different chemical distributions)**: 100% bit-exact on
//! `scripts/descriptor_census_corpus.smi`, 99.6% on
//! `scripts/nci_first_5k_smiles_only.smi`, 94.1% on
//! `scripts/chembl_accuracy_corpus_4999.smi`. Every mismatch traces to
//! chematic's [`chematic_perception::aromaticity::assign_aromaticity_ex`]
//! disagreeing with RDKit's own aromaticity perception on a specific ring
//! system -- confirmed by direct atom/bond-level comparison against a live
//! RDKit oracle for representative failures on each corpus (a benzothiazole
//! fused to a naphtho ring; a benzo-fused tetrahydropyrimidine with two
//! pendant, non-fused phenyl substituents that chematic perceives as part of
//! one aromatic system while RDKit does not; brominated anthraquinone/xanthone
//! cores; S/N-containing bicyclic heterocycles; a Zn coordination complex).
//! This is a pre-existing, already partially-tracked gap in chematic's own
//! aromaticity model (see `AromaticityAlgorithm::RdkitLike`'s own doc comment
//! for its known P/keto-lactam gaps) surfacing through this fingerprint's
//! sensitivity to bond-level aromaticity -- not a defect in this module's own
//! match-enumeration or hashing logic, both independently confirmed correct
//! (per-pattern match *counts* were verified identical to RDKit's own
//! `GetSubstructMatches(..., uniquify=False)` across all 13 patterns on a
//! repro molecule before the aromaticity gap was even found). Full
//! aromaticity-perception parity with RDKit is a separate, much larger
//! undertaking, out of scope for this round.

use rustc_hash::FxHashMap;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::aromaticity::{
    AromaticityAlgorithm, AromaticityModel, assign_aromaticity_ex,
};
use chematic_smarts::{MatchConfig, QueryMolecule, find_matches_with_config, parse_smarts};

use crate::bitvec::BitVec2048;
use crate::rdkit_morgan_hash::hash_combine;

/// RDKit's `pqs[]` (`PatternFingerprints.cpp`), verbatim and in order --
/// 1-based position feeds directly into each pattern's hash seed.
const PATTERNS: [&str; 13] = [
    "[*]~[*]",
    "[*]~[*]~[*]",
    "[R]~1~[R]~[R]~1",
    "[*]~[*](~[*])~[*]",
    "[R]~1[R]~[R]~[R]~1",
    "[*]~[*]~[*](~[*])~[*]",
    "[R]~1~[R]~[R]~[R]~[R]~1",
    "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
    "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
    "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
    "[*]~[R](@[R])@[R](@[R])~[*]",
    "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
    "[*]",
];

/// RDKit's raw `Bond::BondType` enum integer, for the bond types that can
/// occur on a concrete (non-query) `Molecule`'s own bonds. `Aromatic` is
/// handled by the caller (always `12`, regardless of Kekule order -- see
/// module doc comment); `Up`/`Down` are chematic's own stereo-annotated
/// single bonds, hashed as plain `SINGLE` like every other single-bond
/// context in this crate (matching [`crate::rdkit_torsion::num_pi_electrons`]'s
/// own treatment).
fn rdkit_bond_type_code(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::QueryAny => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
        BondOrder::Zero => 21,
        // Query-only bond orders never occur on a concrete Molecule's own
        // bonds (only on QueryMolecule bonds, which are never hashed here).
        BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1,
    }
}

/// Uses [`AromaticityModel::is_bond_aromatic`] rather than checking
/// `bond.order == BondOrder::Aromatic` directly: chematic's SMILES parser
/// preserves a Kekulized *input*'s literal alternating single/double bond
/// orders verbatim (round-trip fidelity), even for a ring the atoms
/// themselves are perceived as aromatic -- unlike RDKit, which always
/// normalizes a perceived-aromatic bond's type to `AROMATIC` regardless of
/// how the input SMILES wrote it. A bond-order check alone therefore misses
/// every aromatic ring bond in Kekule-notation input (confirmed via a live
/// oracle repro: a benzothiazole-fused system parsed from explicit `C1=CC=..`
/// notation kept `Single`/`Double` bond orders in chematic while RDKit
/// reported `AROMATIC` for the same bonds) -- re-perceiving aromaticity here
/// is the fix, not a Kekulized-order-only check.
fn matched_bond_code(
    mol: &Molecule,
    model: &AromaticityModel,
    m: &FxHashMap<usize, AtomIdx>,
    a1: usize,
    a2: usize,
) -> u32 {
    let (bidx, bond) = mol
        .bond_between(m[&a1], m[&a2])
        .expect("a matched pattern bond must exist between its mapped target atoms");
    if model.is_bond_aromatic(bidx) {
        12
    } else {
        rdkit_bond_type_code(bond.order)
    }
}

/// RDKit-bit-exact Pattern fingerprint, matching
/// `rdkit.Chem.PatternFingerprint(mol, fpSize=2048)` (the Python API's own
/// default, `tautomericFingerprint=False`).
pub fn rdkit_pattern_fp(mol: &Molecule) -> BitVec2048 {
    const FP_SIZE: u32 = 2048;
    let mut fp = BitVec2048::new();
    let cfg = MatchConfig {
        uniquify: false,
        ..MatchConfig::default()
    };
    let aromaticity = assign_aromaticity_ex(mol, AromaticityAlgorithm::RdkitLike);

    for (i, smarts) in PATTERNS.iter().enumerate() {
        let p_idx = (i + 1) as u32;
        let query: QueryMolecule =
            parse_smarts(smarts).unwrap_or_else(|e| panic!("built-in pattern '{smarts}': {e}"));
        let matches = find_matches_with_config(&query, mol, &cfg);

        let mut m_idx = p_idx
            .wrapping_add(query.atoms.len() as u32)
            .wrapping_add(query.bonds.len() as u32);
        for _ in &matches {
            m_idx = hash_combine(m_idx, 0xBEEF);
            fp.set((m_idx % FP_SIZE) as usize);
        }

        for m in &matches {
            let mut bit_id = p_idx;
            for qi in 0..query.atoms.len() {
                let an = mol.atom(m[&qi]).element.atomic_number() as u32;
                bit_id = hash_combine(bit_id, an);
            }
            for qb in &query.bonds {
                let code = matched_bond_code(mol, &aromaticity, m, qb.atom1, qb.atom2);
                bit_id = hash_combine(bit_id, code);
            }
            fp.set((bit_id % FP_SIZE) as usize);
        }
    }
    fp
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
        assert_eq!(rdkit_pattern_fp(&m), rdkit_pattern_fp(&m));
    }

    #[test]
    fn different_molecules_differ() {
        assert_ne!(rdkit_pattern_fp(&mol("CCO")), rdkit_pattern_fp(&mol("CCN")));
    }

    #[test]
    fn single_atom_sets_at_least_one_bit() {
        // The "[*]" fragment pattern alone guarantees a nonzero fingerprint
        // for any non-empty molecule.
        assert_ne!(rdkit_pattern_fp(&mol("C")), BitVec2048::new());
    }
}
