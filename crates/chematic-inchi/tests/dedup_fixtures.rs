//! End-to-end fixtures for `chematic_inchi::dedup`: fast canonical-SMILES
//! candidate grouping, verified/reconciled by native standard InChI.
//!
//! Every "same molecule or not" call in this file that is not obvious from
//! the SMILES alone (tautomer/protonation-state/diastereomer pairs) was
//! cross-checked against RDKit 2026.03.3 in an independent venv
//! (`/tmp/chematic-dedupA-rework-venv`, not the repo's shared `.venv` and not
//! any other agent's venv) -- see the PR body for the venv's Python
//! executable path and the exact RDKit calls used. RDKit is a sanity
//! cross-check here; native InChI
//! ([`chematic_inchi::native::standard_inchi`]) is the actual shipped
//! verification mechanism.
//!
//! Requires the `native-inchi` feature (the whole point of this module is to
//! verify via the real oracle, not the pure-Rust approximation).

#![cfg(feature = "native-inchi")]

use chematic_inchi::dedup::{
    DedupRelation, IdentityPolicy, VerifiedGroup, compare, compare_molecules, deduplicate_verified,
};
use chematic_smiles::{parse, random_smiles};

fn mol(smiles: &str) -> chematic_core::Molecule {
    parse(smiles).unwrap_or_else(|e| panic!("parse {smiles:?}: {e}"))
}

// --- 1. Atom-renumbered same molecule --------------------------------------

#[test]
fn atom_renumbered_same_molecule_is_verified_duplicate() {
    // Aspirin, plus a reproducibly-renumbered respelling of the SAME parsed
    // molecule (chematic_smiles::random_smiles permutes atom order, not the
    // chemistry).
    let a = mol("CC(=O)Oc1ccccc1C(=O)O");
    let renumbered = random_smiles(&a, 42);
    let b = mol(&renumbered);
    assert_ne!(
        renumbered, "CC(=O)Oc1ccccc1C(=O)O",
        "renumbering should actually change the string, or this test proves nothing"
    );
    for policy in ALL_POLICIES {
        assert_eq!(
            compare_molecules(&a, &b, policy),
            DedupRelation::VerifiedDuplicate,
            "policy={policy:?}"
        );
    }
}

// --- 2. Alternate valid SMILES spelling of the same molecule ---------------

#[test]
fn alternate_spelling_same_molecule_is_verified_duplicate() {
    let a = mol("CCO");
    let b = mol("OCC");
    for policy in ALL_POLICIES {
        assert_eq!(
            compare_molecules(&a, &b, policy),
            DedupRelation::VerifiedDuplicate,
            "policy={policy:?}"
        );
    }
}

// --- 3. E/Z pair ------------------------------------------------------------

#[test]
fn ez_pair_distinct_under_stereo_aware_policies_merged_under_stereo_ignored() {
    let e = mol("C/C=C/C"); // (E)-but-2-ene
    let z = mol("C/C=C\\C"); // (Z)-but-2-ene

    // Stereo-preserving policies must never merge an E/Z pair.
    assert_eq!(
        compare_molecules(&e, &z, IdentityPolicy::StandardInchiString),
        DedupRelation::Distinct
    );
    assert_eq!(
        compare_molecules(&e, &z, IdentityPolicy::StandardInchiKey),
        DedupRelation::Distinct
    );
    assert_eq!(
        compare_molecules(&e, &z, IdentityPolicy::IsotopeIgnored),
        DedupRelation::Distinct
    );

    // StereoIgnored explicitly ignores stereo, so E and Z become the same
    // identity by design -- but the candidate keys (canonical SMILES) still
    // differ (they encode different directional bonds), so this must be
    // reported as CanonicalSplit, never silently folded into
    // VerifiedDuplicate.
    assert_eq!(
        compare_molecules(&e, &z, IdentityPolicy::StereoIgnored),
        DedupRelation::CanonicalSplit
    );
}

// --- 4. Enantiomer pair ------------------------------------------------------

#[test]
fn enantiomer_pair_distinct_unless_stereo_ignored() {
    let l_ala = mol("N[C@@H](C)C(=O)O"); // L-alanine
    let d_ala = mol("N[C@H](C)C(=O)O"); // D-alanine

    // Enantiomers must never be merged under any stereo-preserving policy.
    for policy in [
        IdentityPolicy::StandardInchiString,
        IdentityPolicy::StandardInchiKey,
        IdentityPolicy::IsotopeIgnored,
    ] {
        assert_eq!(
            compare_molecules(&l_ala, &d_ala, policy),
            DedupRelation::Distinct,
            "policy={policy:?}"
        );
    }
    // Under StereoIgnored, an enantiomeric pair -- which differs *only* in
    // stereo, exactly like the diastereomer case above -- correctly becomes
    // "same identity" by that policy's own definition (mirror-image /t
    // layers land on the same connectivity+H string once stripped). This is
    // NOT a false merge: it is the stated, documented behavior of a policy
    // whose entire purpose is to ignore stereo, and it is reported as
    // CanonicalSplit (verified-but-key-mismatched), never as a silent
    // VerifiedDuplicate. A caller that wants enantiomers kept apart must use
    // one of the three stereo-preserving policies above, never
    // StereoIgnored -- see the `IdentityPolicy::StereoIgnored` doc comment.
    assert_eq!(
        compare_molecules(&l_ala, &d_ala, IdentityPolicy::StereoIgnored),
        DedupRelation::CanonicalSplit
    );
}

// --- 5. Diastereomer pair ----------------------------------------------------

#[test]
fn diastereomer_pair_distinct_unless_stereo_ignored() {
    // L-(2R,3R)-tartaric acid vs meso-(2R,3S)-tartaric acid: same
    // connectivity, different (non-mirror-image) stereo -- genuine
    // diastereomers, never simple respellings of each other.
    let rr = mol("OC(=O)[C@H](O)[C@@H](O)C(=O)O");
    let meso = mol("OC(=O)[C@H](O)[C@H](O)C(=O)O");

    for policy in [
        IdentityPolicy::StandardInchiString,
        IdentityPolicy::StandardInchiKey,
        IdentityPolicy::IsotopeIgnored,
    ] {
        assert_eq!(
            compare_molecules(&rr, &meso, policy),
            DedupRelation::Distinct,
            "policy={policy:?}"
        );
    }
    // Under StereoIgnored, diastereomers -- which differ *only* in stereo --
    // correctly become "same identity" by that policy's own definition, but
    // the candidate keys differ, so this is CanonicalSplit, not a silent
    // VerifiedDuplicate. This is the intended, documented behavior of
    // StereoIgnored, not a bug: it demonstrates why the 4 policies are kept
    // separate rather than collapsed into one "identity" notion.
    assert_eq!(
        compare_molecules(&rr, &meso, IdentityPolicy::StereoIgnored),
        DedupRelation::CanonicalSplit
    );
}

// --- 6. Isotopologue ----------------------------------------------------------

#[test]
fn isotopologue_distinct_unless_isotope_ignored() {
    // Heavy-atom isotope (13C), not H-isotope: chosen deliberately. See the
    // PR body for the separate, independently-verified `chematic-inchi`
    // convert.rs fix required to make H-isotope labels (explicit D/T atoms)
    // round-trip through native InChI at all -- confirmed working here as a
    // regression check, not assumed.
    let plain = mol("CC");
    // NOTE: the bracket atom must carry its own explicit H count
    // (`[13CH3]`, not `[13C]`) -- bracket-atom syntax defaults to ZERO
    // implicit H, unlike the organic subset. `[13C]C` parses to a
    // valence-1 labeled carbon (still accepted by both chematic and RDKit,
    // formula `C2H3`, RDKit-cross-checked) rather than real ethane-13C1;
    // this fixture wants the latter.
    let labeled = mol("[13CH3]C");

    for policy in [
        IdentityPolicy::StandardInchiString,
        IdentityPolicy::StandardInchiKey,
        IdentityPolicy::StereoIgnored,
    ] {
        assert_eq!(
            compare_molecules(&plain, &labeled, policy),
            DedupRelation::Distinct,
            "policy={policy:?}"
        );
    }
    // IsotopeIgnored strips the isotope block, so the two collapse to the
    // same identity -- but different candidate keys (canonical SMILES
    // includes the isotope label), so CanonicalSplit, not a silent merge.
    assert_eq!(
        compare_molecules(&plain, &labeled, IdentityPolicy::IsotopeIgnored),
        DedupRelation::CanonicalSplit
    );
}

#[test]
fn h_isotope_isotopologue_also_distinct_and_reconciled() {
    // The explicit-H-isotope case that a convert.rs fix (merged to `main` as
    // PR #160, rebased onto by this branch) makes representable at all:
    // methane vs methane-d4. Before that fix, native InChI silently produced
    // byte-identical output for both (verified against RDKit's MolToInchi:
    // methane-d4 is `InChI=1S/CH4/h1H4/i1D4`), which would have been an
    // undetectable false merge under StandardInchiString/StandardInchiKey.
    // After the fix, they are correctly Distinct under stereo/isotope-
    // preserving policies and reconciled (CanonicalSplit) under
    // IsotopeIgnored.
    let plain = mol("C");
    let d4 = mol("[2H]C([2H])([2H])[2H]");
    for policy in [
        IdentityPolicy::StandardInchiString,
        IdentityPolicy::StandardInchiKey,
        IdentityPolicy::StereoIgnored,
    ] {
        assert_eq!(
            compare_molecules(&plain, &d4, policy),
            DedupRelation::Distinct,
            "policy={policy:?}"
        );
    }
    assert_eq!(
        compare_molecules(&plain, &d4, IdentityPolicy::IsotopeIgnored),
        DedupRelation::CanonicalSplit
    );
}

// --- 7. Protonation-state pair -----------------------------------------------

#[test]
fn protonation_state_pair_never_merged_under_any_policy() {
    // Acetic acid vs its conjugate base (acetate anion). Standard InChI
    // encodes this via its own /p layer (verified: acetic acid ->
    // `.../h1H3,(H,3,4)`, acetate -> the same plus `/p-1`) -- genuinely
    // different InChI, so this is Distinct under every policy without any
    // special-casing in this module. None of the 4 policies is
    // "charge-ignored", so this must never merge, under any policy.
    let acid = mol("CC(=O)O");
    let anion = mol("CC(=O)[O-]");
    for policy in ALL_POLICIES {
        assert_eq!(
            compare_molecules(&acid, &anion, policy),
            DedupRelation::Distinct,
            "policy={policy:?}"
        );
    }
}

// --- 8. Tautomer pair ---------------------------------------------------------

#[test]
fn tautomer_pair_reconciled_by_inchis_own_normalization_not_by_us() {
    // 2-pyridone / 2-hydroxypyridine: a textbook case where standard InChI's
    // own mobile-H normalization folds two tautomers to the SAME InChI
    // string (verified directly: both produce
    // `InChI=1S/C5H5NO/c7-5-3-1-2-4-6-5/h1-4H,(H,6,7)`). This module does not
    // add any tautomer-merging logic of its own -- the merge, when it
    // happens, is entirely a consequence of native InChI's own semantics.
    // Since the two forms are genuinely different molecular graphs (C=O + ring
    // N-H vs C-OH + aromatic ring N), their canonical SMILES differ, so the
    // correct, honest report is CanonicalSplit (verified same identity, but
    // the fast candidate key did not and could not agree) -- never a silent
    // VerifiedDuplicate that would look like an ordinary respelling.
    let pyridone = mol("O=c1cccc[nH]1");
    let hydroxypyridine = mol("Oc1ccccn1");

    assert_eq!(
        compare_molecules(
            &pyridone,
            &hydroxypyridine,
            IdentityPolicy::StandardInchiString
        ),
        DedupRelation::CanonicalSplit
    );
    assert_eq!(
        compare_molecules(
            &pyridone,
            &hydroxypyridine,
            IdentityPolicy::StandardInchiKey
        ),
        DedupRelation::CanonicalSplit
    );
}

// --- 9. Disconnected salt -----------------------------------------------------

#[test]
fn disconnected_salt_fragment_order_is_verified_duplicate() {
    // Sodium acetate written with the two fragments in opposite order.
    // Verified (native InChI probe) that both orders produce byte-identical
    // standard InChI:
    // `InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1`.
    let order1 = mol("CC(=O)[O-].[Na+]");
    let order2 = mol("[Na+].CC(=O)[O-]");
    for policy in ALL_POLICIES {
        assert_eq!(
            compare_molecules(&order1, &order2, policy),
            DedupRelation::VerifiedDuplicate,
            "policy={policy:?}"
        );
    }
}

// --- 10. Canonical-SMILES residual rows (re-derived at this branch's base
//         commit, NOT trusted from the RFC's own historical numbers) --------
//
// Freshly re-run at this PR's base commit (see PR body for the exact
// summary JSON): 96/5,000 real chematic permutation-invariance failures (all
// independently confirmed molecule-preserving), split 78
// idempotence-only + 18 random-relabeling-only + 0 both. The RFC's own
// historical "348/78" count predates PR #148's partial E/Z-carrier fix,
// which landed on `main` before this branch was cut -- re-deriving rather
// than trusting the doc's headline number was the right call: the true
// current count is materially smaller.
//
// Re-verified again after rebasing onto `main` past PR #160/#157 (neither
// touches canonical_smiles or the aromaticity/SMILES-writer code this
// residual depends on): re-running
// `scripts/canonical_residual_diagnosis.py` produced byte-identical
// `validation/results/canonical_residual_diagnosis*.{jsonl,json}` output --
// same 96/5,000 (78+18+0) split, not a stale carried-forward number.

#[test]
fn residual_row_relabeling_only_reconciled_via_native_inchi() {
    // One of the 18 random-relabeling-only residual cases at this branch's
    // base commit (see validation/results/canonical_residual_diagnosis_summary.json,
    // `permutation_invariance_failures_sample`, `detected_via_random_relabeling`
    // && !`detected_via_idempotence`). Both strings are chematic's own
    // canonical output for two RDKit-relabeled spellings of one input
    // molecule; native InChI confirms they encode the same molecule
    // (verified directly below via the InChI equality, not assumed from the
    // RFC's classification).
    let a = mol("OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC");
    let b = mol(r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c\1c(/c(c1O)O)=N/CCCCC");

    let key_a = chematic_smiles::canonical_smiles(&a);
    let key_b = chematic_smiles::canonical_smiles(&b);
    assert_ne!(
        key_a, key_b,
        "this fixture is only meaningful if the residual actually reproduces \
         (canonical SMILES must differ for this to be a CanonicalSplit case)"
    );

    assert_eq!(
        compare_molecules(&a, &b, IdentityPolicy::StandardInchiString),
        DedupRelation::CanonicalSplit,
        "same-molecule respelling must be reconciled (not silently merged, \
         not reported as Distinct) via native InChI"
    );
}

#[test]
fn residual_row_idempotence_only_reconciled_via_native_inchi() {
    // One of the 78 idempotence-only residual cases at this branch's base
    // commit (validation/results/canonical_residual_diagnosis_summary.json,
    // `detected_via_idempotence` && !`detected_via_random_relabeling`) --
    // Root Cause 3 in the RFC (aromaticity-perception round-trip
    // inconsistency), a different mechanism from the E/Z-carrier case above.
    //
    // Derived directly from the original SMILES by chaining
    // canonical_smiles twice (canonical(s) then canonical(canonical(s))),
    // rather than copying the diagnosis script's recorded output strings --
    // this way the "two different keys" property is reproduced fresh in
    // Rust, not assumed to carry over unchanged from the Python-side JSON.
    let original = mol("O=C(NCCc1c2n(c3ccccc13)CCc1ccccc1-2)C1CCC1");
    let out1 = chematic_smiles::canonical_smiles(&original);
    let a = mol(&out1); // canonical_smiles(&a) == out1's own re-canonicalization
    let out2 = chematic_smiles::canonical_smiles(&a);
    assert_ne!(
        out1, out2,
        "this fixture is only meaningful if the idempotence residual still \
         reproduces at this commit (canonical(canonical(s)) must differ from \
         canonical(s))"
    );

    // `original`'s own candidate key is out1; `a`'s own candidate key is
    // out2 -- different keys, same molecule (just two different points along
    // the same canonicalizer's non-fixed-point chain).
    assert_eq!(
        compare_molecules(&original, &a, IdentityPolicy::StandardInchiString),
        DedupRelation::CanonicalSplit,
        "same-molecule respelling (aromaticity round-trip residual) must be \
         reconciled via native InChI, not reported as Distinct"
    );
}

// --- 11. Synthetic positive control: force CanonicalCollision ---------------

#[test]
fn synthetic_collision_control_fails_closed() {
    // Real canonical-SMILES collisions between two DIFFERENT molecules are
    // not known to occur (0/4,992 in the corpus scan behind
    // docs/canonical_smiles_residual_rfc.md), so this proves
    // CanonicalCollision is reachable and fails closed by INJECTING an
    // identical candidate key for two genuinely different molecules
    // (benzene vs methanol), using the low-level `compare()` primitive
    // rather than `compare_molecules()`. Their real native InChI strings
    // differ, so the tool must report CanonicalCollision -- never merge.
    let benzene = mol("c1ccccc1");
    let methanol = mol("CO");
    let forced_key = "SYNTHETIC-SHARED-KEY";

    for policy in ALL_POLICIES {
        assert_eq!(
            compare(forced_key, &benzene, forced_key, &methanol, policy),
            DedupRelation::CanonicalCollision,
            "policy={policy:?}"
        );
    }
}

// --- VerificationUnavailable / InvalidMolecule -------------------------------

#[test]
fn zero_heavy_atom_molecule_is_invalid_not_unavailable() {
    // [H][H] has no heavy atoms; standard_inchi() returns InvalidInput for
    // it (see chematic_inchi::native's own guard). This is a property of
    // the molecule, not the oracle, so it must be InvalidMolecule, not
    // VerificationUnavailable.
    let h2 = mol("[H][H]");
    let water = mol("O");
    assert_eq!(
        compare_molecules(&h2, &water, IdentityPolicy::StandardInchiString),
        DedupRelation::InvalidMolecule
    );
}

// --- Worst-of-30 permutation robustness --------------------------------------

#[test]
fn worst_of_30_permutations_never_flip_to_distinct_or_collision() {
    // For a molecule that is NOT one of the canonicalizer's known residual
    // cases, every one of 30 reproducibly-generated valid respellings must
    // resolve the SAME way against a fixed reference molecule: never
    // Distinct, never CanonicalCollision, under StandardInchiString
    // identity. This is a robustness check on the dedup module itself (does
    // verification hold up
    // across many different, valid encodings of the same input?), not a
    // canonical-SMILES stability claim.
    let reference = mol("CC(=O)Oc1ccccc1C(=O)O"); // aspirin

    // Request 30; `random_smiles_vect` returns up to 30 *unique* respellings
    // within a bounded attempt budget, which may be fewer than 30 for a
    // molecule with limited DFS-order variety -- report the real count
    // rather than asserting an arbitrary target.
    let respellings = chematic_smiles::random_smiles_vect(&reference, 30, 1000);
    assert!(
        respellings.len() >= 10,
        "expected at least 10 unique respellings for a worst-of-30 test to be \
         meaningful, got only {}",
        respellings.len()
    );
    eprintln!(
        "worst-of-30 permutation test: {} unique respellings exercised",
        respellings.len()
    );

    for (i, s) in respellings.iter().enumerate() {
        let variant_mol = mol(s);
        let relation = compare_molecules(
            &reference,
            &variant_mol,
            IdentityPolicy::StandardInchiString,
        );
        assert!(
            matches!(
                relation,
                DedupRelation::VerifiedDuplicate | DedupRelation::CanonicalSplit
            ),
            "variant #{i} ({s:?}) resolved to {relation:?}, expected VerifiedDuplicate or CanonicalSplit"
        );
    }
}

// --- Regression: explicit-H stereocenter (PR #160 general-case fix) --------

#[test]
fn explicit_h_stereocenter_enantiomers_no_longer_false_duplicate() {
    // Historical known-limitation example from this module's own docs
    // (resolved on `main`, ahead of this rebase, by PR #160's Stereo0D fix
    // for the general single-H-like-substituent case): an explicit
    // graph-H stereocenter substituent used to silently drop its /t layer
    // entirely, making this enantiomer pair a false VerifiedDuplicate.
    // Re-verified here as a regression check, not assumed from the commit
    // message alone.
    let r = mol("[C@](Br)(Cl)(F)[2H]");
    let s = mol("[C@@](Br)(Cl)(F)[2H]");
    assert_eq!(
        compare_molecules(&r, &s, IdentityPolicy::StandardInchiString),
        DedupRelation::Distinct,
        "must be Distinct, not a false VerifiedDuplicate"
    );
}

// --- Guard: stereocenter with 2+ H-like substituents ------------------------
//
// PR #160 fixed the general single-H-like-substituent case (previous
// section) but explicitly left this narrower shape unfixed at the
// `crate::native::convert` level: a stereocentre with TWO H-like
// substituents (e.g. D+T on one carbon, or bracket-H + explicit D) still
// silently drops its `/t` layer. Left unguarded, this module could report
// two genuinely different stereoisomers as a false VerifiedDuplicate, since
// `standard_inchi` itself still returns a (stereo-incomplete) string rather
// than an error. `dedup::identity_verify` detects this input shape
// structurally (`crate::native::has_unrepresentable_multi_h_stereocenter`)
// and fails closed to VerificationUnavailable, for every policy, before
// ever trusting the string.

#[test]
fn two_explicit_h_isotopes_guard_never_false_duplicate() {
    // D+T on the same carbon (matches convert.rs's own
    // `two_h_like_substituents_on_one_centre_drops_stereo_not_corrupts_it`
    // fixture exactly).
    let r = mol("[C@](Br)(F)([2H])[3H]");
    let s = mol("[C@@](Br)(F)([2H])[3H]");
    for policy in ALL_POLICIES {
        let rel = compare_molecules(&r, &s, policy);
        assert_ne!(
            rel,
            DedupRelation::VerifiedDuplicate,
            "policy={policy:?}: must never be a false VerifiedDuplicate"
        );
        assert_eq!(
            rel,
            DedupRelation::VerificationUnavailable,
            "policy={policy:?}: must be explicitly VerificationUnavailable (not Distinct, \
             not silently dropped), got {rel:?}"
        );
    }
}

#[test]
fn bracket_h_plus_explicit_isotope_h_guard_never_false_duplicate() {
    // Bracket-H (implicit) + explicit [2H] on the same carbon -- same guard,
    // a different input shape (sentinel + explicit atom, rather than two
    // explicit atoms). Matches convert.rs's own
    // `bracket_h_plus_explicit_isotope_h_drops_stereo_not_corrupts_it`
    // fixture exactly.
    let r = mol("[C@H](Br)([2H])F");
    let s = mol("[C@@H](Br)([2H])F");
    for policy in ALL_POLICIES {
        let rel = compare_molecules(&r, &s, policy);
        assert_ne!(
            rel,
            DedupRelation::VerifiedDuplicate,
            "policy={policy:?}: must never be a false VerifiedDuplicate"
        );
        assert_eq!(
            rel,
            DedupRelation::VerificationUnavailable,
            "policy={policy:?}: must be explicitly VerificationUnavailable (not Distinct, \
             not silently dropped), got {rel:?}"
        );
    }
}

#[test]
fn two_h_like_substituents_guard_never_lands_in_verified_group() {
    // Same guarded molecules, run through the batch API: none may be
    // grouped with anything, and all 4 must show up in
    // `verification_unavailable`, never in `groups`.
    let mols = [
        mol("[C@](Br)(F)([2H])[3H]"),  // 0
        mol("[C@@](Br)(F)([2H])[3H]"), // 1
        mol("[C@H](Br)([2H])F"),       // 2
        mol("[C@@H](Br)([2H])F"),      // 3
    ];
    for policy in ALL_POLICIES {
        let report = deduplicate_verified(&mols, policy);
        assert!(
            report.groups.is_empty(),
            "policy={policy:?}: none of these guarded molecules may land in a \
             VerifiedGroup: {:?}",
            report.groups
        );
        let mut unavailable = report.verification_unavailable.clone();
        unavailable.sort_unstable();
        assert_eq!(
            unavailable,
            vec![0, 1, 2, 3],
            "policy={policy:?}: all 4 guarded molecules must be VerificationUnavailable"
        );
        assert!(report.invalid_molecules.is_empty(), "policy={policy:?}");
    }
}

#[test]
fn ordinary_single_h_like_enantiomers_still_distinct_not_overcorrected() {
    // The guard must be specific to 2+ H-like substituents -- it must NOT
    // start flagging the ordinary, already-fixed single-H-like cases (plain
    // H, D-only, T-only) as unavailable. Re-checks all three variants, not
    // just the D one already covered above.
    for (r_smi, s_smi) in [
        ("[C@](Br)(Cl)(F)[H]", "[C@@](Br)(Cl)(F)[H]"),
        ("[C@](Br)(Cl)(F)[2H]", "[C@@](Br)(Cl)(F)[2H]"),
        ("[C@](Br)(Cl)(F)[3H]", "[C@@](Br)(Cl)(F)[3H]"),
        ("[C@H](Br)(Cl)F", "[C@@H](Br)(Cl)F"),
    ] {
        let r = mol(r_smi);
        let s = mol(s_smi);
        assert_eq!(
            compare_molecules(&r, &s, IdentityPolicy::StandardInchiString),
            DedupRelation::Distinct,
            "{r_smi:?}/{s_smi:?}: ordinary single-H-like enantiomers must still be Distinct, \
             not over-guarded into VerificationUnavailable"
        );
    }
}

// --- 12. IsotopeIgnored must never mutate the caller's molecule ------------

#[test]
fn isotope_ignored_never_mutates_original_molecule() {
    let original = mol("[2H]C([2H])([2H])[2H]"); // methane-d4
    let atom_count = original.atom_count();
    let isotopes_before: Vec<Option<u16>> = (0..atom_count)
        .map(|i| original.atom(chematic_core::AtomIdx(i as u32)).isotope)
        .collect();
    assert!(
        isotopes_before.iter().any(Option::is_some),
        "fixture only meaningful if the original actually carries an isotope label"
    );

    let plain = mol("C");
    // IsotopeIgnored clones internally and clears isotopes on the clone --
    // `original` (passed by shared reference) must come out unchanged.
    let _ = compare_molecules(&original, &plain, IdentityPolicy::IsotopeIgnored);
    let _ = compare(
        "unused-key-a",
        &original,
        "unused-key-b",
        &plain,
        IdentityPolicy::IsotopeIgnored,
    );

    let isotopes_after: Vec<Option<u16>> = (0..atom_count)
        .map(|i| original.atom(chematic_core::AtomIdx(i as u32)).isotope)
        .collect();
    assert_eq!(
        isotopes_before, isotopes_after,
        "IdentityPolicy::IsotopeIgnored must never mutate the caller's molecule"
    );
}

// --- 13. Batch reconciliation (`deduplicate_verified`) ---------------------

#[test]
fn deduplicate_verified_unifies_residual_row_pair_into_one_group() {
    // Required fixture: two of this project's own canonical-SMILES residual
    // rows -- same molecule, different canonical string (see section 10
    // above; re-verified at this branch's post-rebase HEAD, same 96/5,000
    // split as before the rebase). A different canonical-SMILES key must
    // not stop `deduplicate_verified` from unifying them into ONE
    // `VerifiedGroup`, reconciled across the whole collection (not just
    // detectable via a manual pairwise `compare` call).
    let a = mol("OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC");
    let b = mol(r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c\1c(/c(c1O)O)=N/CCCCC");
    let key_a = chematic_smiles::canonical_smiles(&a);
    let key_b = chematic_smiles::canonical_smiles(&b);
    assert_ne!(
        key_a, key_b,
        "this fixture is only meaningful if the residual actually reproduces"
    );

    let mols = [a, b];
    let report = deduplicate_verified(&mols, IdentityPolicy::StandardInchiString);

    assert_eq!(report.groups.len(), 1, "{:?}", report.groups);
    assert_eq!(report.groups[0].members, vec![0, 1]);
    assert_eq!(
        report.canonical_splits.len(),
        1,
        "{:?}",
        report.canonical_splits
    );
    assert_eq!(report.canonical_splits[0].members, vec![0, 1]);
    assert_eq!(report.canonical_splits[0].canonical_subgroups.len(), 2);
    assert!(report.canonical_collisions.is_empty());
    assert!(report.verification_unavailable.is_empty());
    assert!(report.invalid_molecules.is_empty());
}

#[test]
fn deduplicate_verified_on_synthetic_multi_molecule_corpus() {
    // A small corpus mixing THREE separate true duplicate pairs, several
    // genuinely distinct molecules (including a protonation-state pair,
    // which must never merge under any policy), and one zero-heavy-atom
    // invalid molecule -- not just the single required residual-row
    // fixture above.
    let mols: Vec<chematic_core::Molecule> = vec![
        mol("CCO"),                   // 0: ethanol
        mol("OCC"),                   // 1: ethanol, respelled -> dup of 0
        mol("CCN"),                   // 2: ethylamine, distinct singleton
        mol("CC(=O)O"),               // 3: acetic acid
        mol("OC(=O)C"),               // 4: acetic acid, respelled -> dup of 3
        mol("CC(=O)[O-]"),            // 5: acetate anion, distinct (protonation state)
        mol("c1ccccc1"),              // 6: benzene, aromatic form
        mol("C1=CC=CC=C1"),           // 7: benzene, Kekule form -> dup of 6
        mol("[H][H]"),                // 8: invalid (zero heavy atoms)
        mol("CC(=O)Oc1ccccc1C(=O)O"), // 9: aspirin, distinct singleton
    ];

    let report = deduplicate_verified(&mols, IdentityPolicy::StandardInchiString);

    assert!(
        report.verification_unavailable.is_empty(),
        "{:?}",
        report.verification_unavailable
    );
    assert_eq!(report.invalid_molecules, vec![8]);
    assert!(
        report.canonical_collisions.is_empty(),
        "no real canonical-SMILES collision is expected in this corpus: {:?}",
        report.canonical_collisions
    );

    assert_eq!(
        report.groups.len(),
        3,
        "expected exactly 3 duplicate pairs: {:?}",
        report.groups
    );
    let find_group = |idx: usize| -> Option<&VerifiedGroup> {
        report.groups.iter().find(|g| g.members.contains(&idx))
    };
    assert_eq!(find_group(0).unwrap().members, vec![0, 1]);
    assert_eq!(find_group(3).unwrap().members, vec![3, 4]);
    assert_eq!(find_group(6).unwrap().members, vec![6, 7]);

    // Distinct singletons must never appear in any group.
    for idx in [2usize, 5, 9] {
        assert!(
            find_group(idx).is_none(),
            "index {idx} must not be grouped with anything"
        );
    }

    // Deterministic order: groups sorted by lowest member index.
    let starts: Vec<usize> = report.groups.iter().map(|g| g.members[0]).collect();
    let mut sorted_starts = starts.clone();
    sorted_starts.sort_unstable();
    assert_eq!(starts, sorted_starts);
}

const ALL_POLICIES: [IdentityPolicy; 4] = [
    IdentityPolicy::StandardInchiString,
    IdentityPolicy::StandardInchiKey,
    IdentityPolicy::StereoIgnored,
    IdentityPolicy::IsotopeIgnored,
];
