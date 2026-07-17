//! Kekule-S0: verify `chematic_core::apply_kekule` preserving
//! `stereo_neighbor_order` (see `crates/chematic-core/src/kekulization.rs`)
//! actually fixes atom-mapped CIP R/S assignment for molecules routed
//! through an explicit kekulize+apply_kekule round trip before CIP
//! assignment -- exactly what `chematic-perception`'s
//! `apply_aromaticity_rdkit_parity_experimental` does internally.
//!
//! Before the fix, `mol.stereo_neighbor_order(idx)` was `None` on any
//! post-`apply_kekule` molecule, which both CIP engines treat as
//! "insufficient information" and *skip* (not mislabel) --
//! `assign_cip_accurate_experimental` pushes `SkipReason::NotFourSubstituents`
//! (a misleading name for this cause: the substituent count is fine, the
//! neighbor *order* is what's missing) and the legacy engine's assignment
//! for that atom is simply absent. Both engines are checked here, per atom,
//! comparing the RAW aromatic-form input against the same molecule after an
//! explicit kekulize+apply_kekule round trip.

use std::collections::HashMap;

use chematic_chem::assign_cip as assign_cip_legacy;
use chematic_cip::{CipBudget, assign_cip_accurate_experimental};
use chematic_core::{AtomIdx, CipCode, Molecule, apply_kekule, kekulize};

fn kekulized(mol: &Molecule) -> Molecule {
    let k = kekulize(mol).expect("kekulizable");
    apply_kekule(mol, &k)
}

fn legacy_map(mol: &Molecule) -> HashMap<u32, Option<CipCode>> {
    let assignment = assign_cip_legacy(mol);
    (0..mol.atom_count())
        .map(|i| {
            let idx = AtomIdx(i as u32);
            (i as u32, assignment.get(idx))
        })
        .collect()
}

fn accurate_map(mol: &Molecule) -> HashMap<u32, Option<CipCode>> {
    let assignment = assign_cip_accurate_experimental(mol, CipBudget::default_budget())
        .expect("CIP budget not exceeded for these small test molecules");
    let mut map: HashMap<u32, Option<CipCode>> =
        (0..mol.atom_count()).map(|i| (i as u32, None)).collect();
    for (idx, code) in &assignment.assignments {
        map.insert(idx.0, Some(*code));
    }
    map
}

#[track_caller]
fn assert_cip_invariant_under_kekulize(smi: &str) {
    let raw = chematic_smiles::parse(smi).expect("valid SMILES");
    let kek = kekulized(&raw);
    assert_eq!(
        raw.atom_count(),
        kek.atom_count(),
        "{smi}: apply_kekule must preserve atom count/index mapping"
    );

    let legacy_before = legacy_map(&raw);
    let legacy_after = legacy_map(&kek);
    assert_eq!(
        legacy_before, legacy_after,
        "{smi}: legacy chematic_chem::assign_cip differs before/after kekulize+apply_kekule"
    );

    let accurate_before = accurate_map(&raw);
    let accurate_after = accurate_map(&kek);
    assert_eq!(
        accurate_before, accurate_after,
        "{smi}: assign_cip_accurate_experimental differs before/after kekulize+apply_kekule"
    );
}

#[test]
fn simple_chiral_centers() {
    for smi in [
        "N[C@@H](C)C(=O)O",
        "N[C@H](C)C(=O)O",
        "[C@@H]1(N)CCCC1",
        "F[C@H]1CCCCC1",
    ] {
        assert_cip_invariant_under_kekulize(smi);
    }
}

#[test]
fn chiral_substituent_on_aromatic_ring() {
    // The exact shape that exposed the bug: a stereocenter whose
    // stereo_neighbor_order is defined relative to an aromatic ring atom
    // that itself needs kekulization.
    for smi in [
        "c1ccc(cc1)[C@H](F)Cl",
        "c1ccc(cc1)[C@@H](F)Cl",
        "Cc1cn([C@H]2CCCC[C@@H]2C)c(=O)n1C",
    ] {
        assert_cip_invariant_under_kekulize(smi);
    }
}

#[test]
fn real_corpus_case_multi_stereocenter() {
    // The Aromaticity-A1-1b-1 "experimental_only" canonical-round-trip
    // instability case that first surfaced this bug: 3 stereocenters, an
    // aromatic ring directly bonded to two of them.
    assert_cip_invariant_under_kekulize(
        "CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1",
    );
}

#[test]
fn pre_fix_regression_would_have_skipped_not_mislabeled() {
    // Documents the exact failure shape the fix closes: without
    // stereo_neighbor_order, both engines treat the atom as
    // "insufficient information" rather than emitting a wrong label. This
    // test would have failed pre-fix with `after[chiral_idx] == None` while
    // `before[chiral_idx] == Some(..)`.
    let raw = chematic_smiles::parse("c1ccc(cc1)[C@H](F)Cl").expect("valid SMILES");
    let chiral_idx = raw
        .atoms()
        .find(|(_, a)| a.chirality != chematic_core::Chirality::None)
        .map(|(idx, _)| idx.0)
        .expect("test setup sanity: exactly one chiral atom");
    let kek = kekulized(&raw);

    let before = accurate_map(&raw);
    let after = accurate_map(&kek);
    assert!(
        before[&chiral_idx].is_some(),
        "test setup sanity: atom {chiral_idx} should be a resolvable stereocenter"
    );
    assert_eq!(
        before[&chiral_idx], after[&chiral_idx],
        "atom {chiral_idx}'s CIP code must survive kekulize+apply_kekule, not silently become None"
    );
}
