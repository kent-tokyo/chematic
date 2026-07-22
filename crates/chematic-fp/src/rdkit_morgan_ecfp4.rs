//! Production, fallible, RDKit-bit-exact ECFP4 path.
//!
//! Promotes [`crate::rdkit_morgan_hash`]'s source-verified hash port (Milestone M4-A0, PR #124)
//! to a real public API, restricted to the single option envelope M4-A0 actually verified
//! numerically against a live RDKit oracle: radius = 2 (ECFP4), fpSize = 2048,
//! `includeRedundantEnvironments = false` (RDKit's default, suppressed lifecycle),
//! `useChirality = false`, `useBondTypes = true`. **Not** ECFP6/radius = 3 — M4-A0 never
//! compared radius 3 against the oracle, so claiming bit-exactness there would be an
//! unverified extrapolation.
//!
//! Uses [`apply_aromaticity_rdkit_parity_experimental`] internally as a fallible `Result`
//! step. There is no fallback to production Hückel aromaticity anywhere in this module's
//! public path — see the project's own
//! `feedback_fallback_pooling_measurement_error` lesson (M4-A0's original report pooled a
//! Hückel-fallback result into an "RDKit-parity" success count and had to be corrected): this
//! module's whole claim is bit-exactness against real RDKit, so silently substituting a
//! different aromaticity engine on `Err` would silently invalidate that claim on exactly the
//! inputs where a caller most needs to know it doesn't hold.
//!
//! No entry point in this module accepts a pre-aromatized [`Molecule`] — aromaticity
//! perception always happens internally, via the same RDKit-parity engine every call, so a
//! caller cannot bypass it with different (possibly Hückel-derived) flags and silently receive
//! a fingerprint that no longer carries the bit-exactness guarantee.

use chematic_core::{BondIdx, BondOrder, Molecule};
use chematic_perception::AromaticityError;
use rustc_hash::FxHashMap;

use crate::bitvec::BitVec2048;
use crate::rdkit_morgan_hash::{checked_bond_invariant, expand_one_pass};

const ECFP4_RADIUS: u32 = 2;
const ECFP4_FP_SIZE: usize = 2048;

/// Every RDKit-hash-exact ECFP4 view of a molecule, computed from one shared expansion pass
/// (RDKit's `includeRedundantEnvironments = false` lifecycle) — not independently recomputed
/// per field.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RdkitMorganEcfp4 {
    /// 2048-bit folded fingerprint (`raw_identifier % 2048`, OR-combined).
    pub fingerprint: BitVec2048,
    /// Raw (unfolded) identifier → emission count (RDKit's `GetSparseCountFingerprint`
    /// shape). An identifier's count is how many distinct `(atom, radius)` environments
    /// emitted it, which can exceed 1 on an accidental hash collision.
    pub sparse_counts: FxHashMap<u32, u32>,
    /// Raw identifier → the `(atom_idx, radius)` environments that produced it (RDKit's
    /// `AdditionalOutput.GetBitInfoMap()`, unfolded).
    pub raw_bit_info: FxHashMap<u32, Vec<(u32, u32)>>,
    /// Folded bit (`0..2048`) → the `(atom_idx, radius)` environments that set it (RDKit's
    /// `AdditionalOutput.GetBitInfoMap()` on the folded fingerprint).
    pub folded_bit_info: FxHashMap<usize, Vec<(u32, u32)>>,
}

/// Why [`rdkit_morgan_ecfp4_experimental`] could not produce a result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitMorganError {
    /// RDKit-parity aromaticity preprocessing failed — see the wrapped
    /// [`AromaticityError`] for the specific reason (kekulization failure or an internal
    /// invariant violation). No fallback to another aromaticity engine is ever attempted; see
    /// the module docs.
    Aromaticity(AromaticityError),
    /// `bond_idx`'s `order` has no real RDKit `Bond::BondType` counterpart (only chematic's
    /// SMARTS-query-only `BondOrder` variants — cannot occur for a SMILES-parsed molecule, but
    /// is checked explicitly rather than assumed unreachable).
    UnsupportedBondOrder { bond_idx: BondIdx, order: BondOrder },
    /// A post-computation sanity check failed that should never happen for chemically valid
    /// input — surfaced as an error rather than a panic or a silently wrong fingerprint.
    InternalInvariantViolation { reason: String },
}

impl std::fmt::Display for RdkitMorganError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RdkitMorganError::Aromaticity(e) => write!(f, "rdkit-exact ecfp4: aromaticity: {e}"),
            RdkitMorganError::UnsupportedBondOrder { bond_idx, order } => write!(
                f,
                "rdkit-exact ecfp4: bond {bond_idx:?} has no RDKit BondType counterpart: {order:?}"
            ),
            RdkitMorganError::InternalInvariantViolation { reason } => {
                write!(
                    f,
                    "rdkit-exact ecfp4: internal invariant violation: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for RdkitMorganError {}

impl From<AromaticityError> for RdkitMorganError {
    fn from(e: AromaticityError) -> Self {
        RdkitMorganError::Aromaticity(e)
    }
}

/// RDKit-bit-exact ECFP4 (radius = 2, 2048 bits, `includeRedundantEnvironments = false`,
/// `useChirality = false`, `useBondTypes = true`) — bit-exact against real RDKit for every
/// input where preprocessing succeeds, verified on the full M4-A0 corpus (5,046/5,046
/// `rdkit_parity_success` rows, 100% agreement across raw identifiers, sparse counts, folded
/// bits, and bitInfo). Kekulization/aromaticity failures are reported as `Err`, never silently
/// degraded to a different, non-exact result — see the module docs.
pub fn rdkit_morgan_ecfp4_experimental(
    mol: &Molecule,
) -> Result<RdkitMorganEcfp4, RdkitMorganError> {
    let aromatized = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol)?;

    let mut result = RdkitMorganEcfp4::default();
    if aromatized.atom_count() == 0 {
        return Ok(result);
    }

    let ring_set = chematic_perception::find_sssr(&aromatized);
    let bond_count = aromatized.bond_count();
    let mut bond_invariants = Vec::with_capacity(bond_count);
    for b in 0..bond_count {
        let bond_idx = BondIdx(b as u32);
        let order = aromatized.bond(bond_idx).order;
        let invariant = checked_bond_invariant(order)
            .ok_or(RdkitMorganError::UnsupportedBondOrder { bond_idx, order })?;
        bond_invariants.push(invariant);
    }

    let emitted = expand_one_pass(&aromatized, &ring_set, &bond_invariants, ECFP4_RADIUS, true);

    for ((atom_idx, radius), raw_id) in emitted {
        let folded = (raw_id as usize) % ECFP4_FP_SIZE;
        result.fingerprint.set(folded);
        *result.sparse_counts.entry(raw_id).or_insert(0) += 1;
        result
            .raw_bit_info
            .entry(raw_id)
            .or_default()
            .push((atom_idx, radius));
        result
            .folded_bit_info
            .entry(folded)
            .or_default()
            .push((atom_idx, radius));
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::MoleculeBuilder;
    use chematic_smiles::parse;

    #[test]
    fn benzene_matches_rdkit_ground_truth_radius0_and_folds_correctly() {
        let mol = parse("c1ccccc1").unwrap();
        let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        // Radius-0 ground truth pinned in `rdkit_morgan_hash`'s own test module.
        assert!(result.raw_bit_info.contains_key(&3218693969));
        for &raw_id in result.raw_bit_info.keys() {
            let folded = (raw_id as usize) % ECFP4_FP_SIZE;
            assert!(result.fingerprint.get(folded));
        }
    }

    #[test]
    fn hermetic_equivalence_to_diagnostic_default_lifecycle() {
        // The production path's raw_bit_info (inverted) must equal
        // `rdkit_morgan_raw_trace`'s `raw_identifier_default`-Some entries on the same
        // already-aromatized molecule -- proves the promotion didn't drift from the path
        // M4-A0 already validated at 5,046/5,046.
        for smi in ["c1ccncc1", "CC(=O)[O-]", "c1ccc2ccccc2c1", "CC"] {
            let mol = parse(smi).unwrap();
            let aromatized =
                chematic_perception::apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
            let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();

            let trace = crate::rdkit_morgan_hash::rdkit_morgan_raw_trace(&aromatized, 2);
            let mut expected: Vec<((u32, u32), u32)> = trace
                .iter()
                .filter_map(|e| {
                    e.raw_identifier_default
                        .map(|rid| ((e.atom_idx, e.radius), rid))
                })
                .collect();
            expected.sort_unstable();

            let mut got: Vec<((u32, u32), u32)> = result
                .raw_bit_info
                .iter()
                .flat_map(|(&raw_id, envs)| envs.iter().map(move |&(a, r)| ((a, r), raw_id)))
                .collect();
            got.sort_unstable();

            assert_eq!(got, expected, "mismatch for {smi}");
        }
    }

    /// Was `kekule_pyridinium_reports_kekulization_failed_not_a_fallback_result`, pinned to
    /// pyridinium's `c1cc[nH+]cc1` as its "kekulize fails" example. `fix/kekulize-charge-aware-k1`
    /// (chematic-core's `atom_must_be_matched`) fixed pyridinium's kekulization -- it now
    /// succeeds and is bit-exact against RDKit (see
    /// `validation/results/ecfp4_bitexact_matrix_summary.json`'s `charged_kekulize_fail` bucket,
    /// 6/6 `verified_bit_exact`) -- so it's no longer a valid "kekulize fails" example. This test
    /// checks a general contract (a real aromaticity/kekulization failure must surface as
    /// `Err`, never silently succeed), not a pyridinium-specific fact, so it's renamed and
    /// re-pointed at a molecule that still genuinely fails: the same bridgehead-N purine-like
    /// ring `chematic_perception::rdkit_parity`'s own
    /// `production_api_reports_kekulize_failure_not_panic` test uses, confirmed still
    /// `KekulizationFailed` after K1 (also confirmed by the 5,000-molecule corpus diff in
    /// `validation/results/kekulize_charge_aware_k1_corpus_diff.json`: this is the one
    /// pre-existing failure that's unchanged before/after K1).
    #[test]
    fn kekule_bridgehead_n_purine_reports_kekulization_failed_not_a_fallback_result() {
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = parse(smi).unwrap();
        match rdkit_morgan_ecfp4_experimental(&mol) {
            Err(RdkitMorganError::Aromaticity(AromaticityError::KekulizationFailed { .. })) => {}
            other => panic!("expected Aromaticity(KekulizationFailed) for {smi}, got {other:?}"),
        }
    }

    #[test]
    fn degree_zero_atom_emits_only_its_radius0_identifier() {
        // Radius 0 is unconditional for every atom (see `rdkit_morgan_hash`'s own
        // `degree_zero_atom_never_appears_past_radius_zero`) -- only radius >= 1 is suppressed
        // by degree-0 death.
        let mol = parse("[Cl-]").unwrap();
        let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        assert_eq!(result.sparse_counts.len(), 1);
        let envs: Vec<_> = result.raw_bit_info.values().flatten().copied().collect();
        assert_eq!(envs, vec![(0, 0)]);
    }

    #[test]
    fn empty_molecule_yields_empty_result_not_an_error() {
        let mol = MoleculeBuilder::new().build();
        let result = rdkit_morgan_ecfp4_experimental(&mol).unwrap();
        assert_eq!(result, RdkitMorganEcfp4::default());
    }

    /// `BondOrder::Query*` has no real RDKit `Bond::BondType` counterpart (see
    /// [`checked_bond_invariant`]'s doc comment) and cannot arise from `parse()` -- built
    /// programmatically here specifically to prove the explicit-`Err` path, not a guessed
    /// mapping, actually fires.
    #[test]
    fn query_bond_order_is_an_explicit_unsupported_bond_order_error_not_a_guess() {
        use chematic_core::{Atom, Element};

        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let b = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(a, b, BondOrder::QueryAny).unwrap();
        let mol = builder.build();

        match rdkit_morgan_ecfp4_experimental(&mol) {
            Err(RdkitMorganError::UnsupportedBondOrder { bond_idx, order }) => {
                assert_eq!(bond_idx, BondIdx(0));
                assert_eq!(order, BondOrder::QueryAny);
            }
            other => panic!("expected UnsupportedBondOrder, got {other:?}"),
        }
    }

    /// Positive control #9 (Phase B spec): a reintroduced Hückel fallback must be caught by a
    /// test. Simulated here without editing production code — this test independently proves
    /// that on the known kekulization-gap molecule, plain (Hückel-based) `apply_aromaticity`
    /// produces a *different* atom/bond aromaticity outcome than what an `Err` from this
    /// module's real path implies, so a hypothetical silent fallback substituting the former
    /// for the latter would be numerically detectable, not just contractually forbidden.
    ///
    /// Fixture swapped from pyridinium's `c1cc[nH+]cc1` (same reason as
    /// `kekule_bridgehead_n_purine_reports_kekulization_failed_not_a_fallback_result` above:
    /// `fix/kekulize-charge-aware-k1` made pyridinium kekulize successfully, so it's no longer a
    /// molecule this test's premise -- "the real path fails" -- holds for) to the same
    /// still-failing bridgehead-N purine-like ring. The expected Hückel pattern below is the
    /// actual observed output of `chematic_perception::apply_aromaticity` on this molecule (20
    /// atoms), captured via a throwaway probe run against this exact commit, not guessed.
    #[test]
    fn hueckel_fallback_would_be_detectable_if_silently_reintroduced() {
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = parse(smi).unwrap();

        // The real, fallible path must fail -- no result to compare against RDKit at all.
        let real = rdkit_morgan_ecfp4_experimental(&mol);
        assert!(matches!(real, Err(RdkitMorganError::Aromaticity(_))));

        // A hypothetical silent fallback would instead run production Hückel aromaticity and
        // report success. Prove that path is reachable and produces a concrete, observable
        // aromatic-atom partition, so such a substitution is not merely "different code path"
        // but "numerically distinguishable, and thus catchable" if it were ever reintroduced.
        let hueckel_fallback_mol = chematic_perception::apply_aromaticity(&mol);
        let hueckel_aromatic_atoms: Vec<bool> = (0..hueckel_fallback_mol.atom_count())
            .map(|i| {
                hueckel_fallback_mol
                    .atom(chematic_core::AtomIdx(i as u32))
                    .aromatic
            })
            .collect();
        // Hückel perceives 10 of this molecule's 20 atoms as aromatic (a mixed pattern, not
        // uniformly true/false) -- i.e. a silent fallback would have returned Ok(..) with a
        // fingerprint reflecting this partition here, directly contradicting the real path's
        // Err. This assertion is what would fail if a fallback were reintroduced and this test
        // were updated to call the (currently nonexistent) fallback path instead of the real one.
        assert_eq!(
            hueckel_aromatic_atoms,
            vec![
                false, true, true, true, true, false, true, true, true, true, false, false, false,
                false, false, true, true, true, true, false
            ]
        );
    }
}
