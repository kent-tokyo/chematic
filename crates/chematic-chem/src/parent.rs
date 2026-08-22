//! Parent identity (ROADMAP.md Phase 2, round 2B).
//!
//! A **Parent** is an idempotent, deterministic reduction of one axis of
//! molecular variability to one representative structure, meant to be used
//! as a grouping/dedup key -- see
//! `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` section 4.3.
//!
//! [`fragment_parent`], [`charge_parent`], [`isotope_parent`], and
//! [`stereo_parent`] live in [`crate::standardize`], next to the low-level
//! primitives they wrap. [`tautomer_parent`](crate::tautomer::tautomer_parent)
//! lives in [`crate::tautomer`]. This module holds the shared result types
//! and [`super_parent`], which composes all five in one fixed order.

#![forbid(unsafe_code)]

use chematic_core::Molecule;

use crate::standardize::{
    TransformationRecord, charge_parent, fragment_parent, isotope_parent, stereo_parent,
};
use crate::tautomer::{TautomerLimits, tautomer_parent};

/// Why a Parent computation could not reach a definite answer.
///
/// A new, orthogonal type -- not an extension of [`crate::standardize::PipelineStatus`],
/// which answers "did the molecule change," a different question from "did
/// the computation reach a definite answer, and if not, why not."
/// `#[non_exhaustive]`: new statuses may be added without a breaking change.
///
/// There is no `Canceled` variant: no cancellation mechanism (token or
/// callback) exists yet, and a state nothing can ever produce would look
/// like a supported feature to a caller who matches on it. Add it back only
/// alongside the mechanism that produces it.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ParentComputationStatus {
    Completed,
    MaxTransformsReached,
    MaxTautomersReached,
    /// Outside the determinism guarantee `max_transforms`/`max_tautomers`
    /// carry -- see [`TautomerLimits`]'s `timeout_ms` doc comment.
    TimedOut,
    Abstained(AbstainReason),
    InvalidInput(InvalidInputReason),
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AbstainReason {
    NoConfidentOrganicParent,
    AmbiguousFragmentSelection,
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum InvalidInputReason {
    EmptyMolecule,
    UnparsableStructure,
}

/// Result of a Parent computation that has a result-state question
/// ([`tautomer_parent`](crate::tautomer::tautomer_parent) and
/// [`super_parent`]). `fragment_parent`/`charge_parent`/`isotope_parent`/
/// `stereo_parent` keep returning `(Molecule, TransformationRecord)`
/// directly (Phase 1's shape, unchanged) -- they are simple mechanical
/// transforms with no budget/timeout/abstain question.
///
/// `#[non_exhaustive]` struct: fields may be added without a breaking
/// change. Does not implement `Serialize`/`Deserialize` even under the
/// `serde` feature: it embeds a raw `Molecule`, which does not implement
/// them.
#[non_exhaustive]
#[derive(Clone)]
pub struct ParentResult {
    pub molecule: Molecule,
    pub status: ParentComputationStatus,
    pub audit: ParentAudit,
}

/// One Parent function's audit trail, or (for [`super_parent`]) every
/// stage's, in order.
///
/// Does not implement `Serialize`/`Deserialize` even under the `serde`
/// feature: it embeds [`crate::tautomer::TautomerAuditRecord`], which does
/// not implement them (see that type's own doc comment).
#[derive(Debug, Clone)]
pub enum ParentAudit {
    /// fragment_parent / charge_parent / isotope_parent / stereo_parent.
    Transformation(TransformationRecord),
    /// tautomer_parent.
    Tautomer(crate::tautomer::TautomerAuditRecord),
    /// super_parent: one entry per stage, in the fixed order
    /// fragment_parent -> charge_parent -> isotope_parent -> stereo_parent
    /// -> tautomer_parent.
    Composed(Vec<ParentAudit>),
}

/// Compose all five Parent functions in one fixed order:
/// `fragment_parent` -> `charge_parent` -> `isotope_parent` ->
/// `stereo_parent` -> `tautomer_parent`, each stage's output feeding the
/// next stage's input.
///
/// This order is deliberate and does **not** follow
/// `StandardizationPipeline::run`'s stage order (which neutralizes charges
/// *before* fragment selection, for a different reason -- see that
/// pipeline's own comment): `super_parent` selects the representative
/// fragment *first*, because every subsequent Parent step should operate on
/// one fragment, not a not-yet-resolved multi-fragment molecule. See
/// `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` section 4.3.
///
/// `status`/the final `molecule` come from the `tautomer_parent` stage
/// (the only stage with a result-state question); `audit` is
/// `ParentAudit::Composed` with exactly 5 entries, one per stage, so a
/// caller can inspect every intermediate result, not just the final one.
pub fn super_parent(mol: &Molecule, limits: &TautomerLimits) -> ParentResult {
    if mol.atom_count() == 0 {
        return ParentResult {
            molecule: mol.clone(),
            status: ParentComputationStatus::InvalidInput(InvalidInputReason::EmptyMolecule),
            audit: ParentAudit::Composed(Vec::new()),
        };
    }

    let (after_fragment, fragment_record) = fragment_parent(mol);
    let (after_charge, charge_record) = charge_parent(&after_fragment);
    let (after_isotope, isotope_record) = isotope_parent(&after_charge);
    let (after_stereo, stereo_record) = stereo_parent(&after_isotope);
    let tautomer_result = tautomer_parent(&after_stereo, limits);

    ParentResult {
        molecule: tautomer_result.molecule,
        status: tautomer_result.status,
        audit: ParentAudit::Composed(vec![
            ParentAudit::Transformation(fragment_record),
            ParentAudit::Transformation(charge_record),
            ParentAudit::Transformation(isotope_record),
            ParentAudit::Transformation(stereo_record),
            tautomer_result.audit,
        ]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::{canonical_smiles, parse};

    // -- Phase 2 round-2B super_parent fixture test ---------------------------
    // Mirrors validation/tautomer_parent_identity_phase2_fixtures.jsonl's
    // tp2-23-super-parent-composed. See
    // docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 5.

    #[test]
    fn tp2_23_super_parent_composed_pins_every_intermediate_stage() {
        let input = "[NH3+][C@@H]([2H])C(=O)[O-].Cl";
        let mol = parse(input).unwrap();

        let (after_fragment, _) = fragment_parent(&mol);
        assert_eq!(
            canonical_smiles(&after_fragment),
            canonical_smiles(&parse("[2H][C@@H](C(=O)[O-])[NH3+]").unwrap()),
            "after fragment_parent: HCl salt dropped"
        );

        let (after_charge, _) = charge_parent(&after_fragment);
        assert_eq!(
            canonical_smiles(&after_charge),
            canonical_smiles(&parse("[2H][C@@H](C(O)=O)N").unwrap()),
            "after charge_parent: zwitterion neutralized"
        );

        let (after_isotope, _) = isotope_parent(&after_charge);
        assert_eq!(
            canonical_smiles(&after_isotope),
            canonical_smiles(&parse("[H][C@@H](C(O)=O)N").unwrap()),
            "after isotope_parent: 2H label dropped"
        );

        let (after_stereo, _) = stereo_parent(&after_isotope);
        assert_eq!(
            canonical_smiles(&after_stereo),
            canonical_smiles(&parse("[H]C(C(O)=O)N").unwrap()),
            "after stereo_parent: @@ dropped"
        );

        let result = super_parent(&mol, &TautomerLimits::default());
        assert_eq!(result.status, ParentComputationStatus::Completed);
        assert_eq!(
            canonical_smiles(&result.molecule),
            canonical_smiles(&parse("[H]C(C(O)=O)N").unwrap()),
            "super_parent's final output matches the pinned end-to-end chain"
        );
        match &result.audit {
            ParentAudit::Composed(stages) => assert_eq!(stages.len(), 5, "one entry per stage"),
            other => panic!("expected ParentAudit::Composed, got {other:?}"),
        }
    }

    #[test]
    fn super_parent_on_empty_molecule_is_invalid_input() {
        let empty = chematic_core::MoleculeBuilder::new().build();
        let result = super_parent(&empty, &TautomerLimits::default());
        assert_eq!(
            result.status,
            ParentComputationStatus::InvalidInput(InvalidInputReason::EmptyMolecule)
        );
    }

    #[test]
    fn super_parent_noop_on_toluene() {
        // tp2-holdout-04: an already-neutral/unlabeled/achiral/single-fragment
        // input must pass through all four mechanical stages unchanged.
        let mol = parse("Cc1ccccc1").unwrap();
        let result = super_parent(&mol, &TautomerLimits::default());
        assert_eq!(result.status, ParentComputationStatus::Completed);
        assert_eq!(canonical_smiles(&result.molecule), canonical_smiles(&mol));
    }
}
