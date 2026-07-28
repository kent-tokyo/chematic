//! Negative controls for the torsion-knowledge v2 layer (Wave 2 spec §13).
//!
//! Every test in this file demonstrates that a specific bad behavior does
//! **not** occur, and (for the "harness itself" requirement) that the
//! checks proving so are not vacuous -- a deliberately corrupted expected
//! value is shown to actually make the check fail, via
//! `std::panic::catch_unwind` around an intentionally-wrong assertion
//! (the only way to demonstrate "this check can fail" inside a test suite
//! that itself must pass).
//!
//! Additional negative controls needing access to private internals
//! (`resolve_same_tier`, `Candidate`) live in
//! `etkdg_knowledge::matcher::tests` and `etkdg_knowledge::energy::tests`
//! instead (not reachable from this integration-test crate) -- see those
//! modules for: rule-order invariance, silent same-tier-conflict
//! acceptance, artificially-corrupted-geometry detection (ring tear,
//! chirality flip), and non-finite-energy exclusion.

use chematic_3d::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};
use chematic_3d::etkdg_knowledge::{
    RingMembership, RingMembershipIndex, TorsionKnowledgeConfig, TorsionKnowledgeSource,
    build_torsion_knowledge, classify_bond, macrocycle_14_bound_adjustments,
};
use chematic_core::AtomIdx;
use chematic_smiles::parse;

/// Spec §13: "including a non-1,4 pair in macrocycle 1-4 adjustments" must
/// FAIL to occur. Every returned pair must be a genuine 1-4 (not directly
/// bonded, i.e. not a 1-2 or a short-circuited 1-3).
#[test]
fn macrocycle_14_adjustments_never_include_a_non_14_pair() {
    for smiles in ["C1CCCCCCCCCCC1", "O=C1CCCCCCCCCCN1", "O1CCOCCOCCOCC1"] {
        let mol = parse(smiles).unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let adjustments = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(
            !adjustments.is_empty(),
            "{smiles}: expected at least one 1-4 pair"
        );
        for adj in &adjustments {
            let (a, b) = adj.atom_pair;
            assert_ne!(a, b, "{smiles}: pair must be two distinct atoms");
            assert!(
                mol.bond_between(a, b).is_none(),
                "{smiles}: pair ({},{}) is directly bonded -- not a genuine 1-4",
                a.0,
                b.0
            );
        }
    }
}

/// Spec §13: "applying an acyclic torsion rule to an aromatic ring bond"
/// must FAIL to occur, checked across a broader fixture set than the
/// single-molecule unit test in `matcher.rs`.
#[test]
fn standard_tier_never_fires_on_aromatic_ring_bonds_across_fixtures() {
    for smiles in ["c1ccccc1", "c1ccncc1", "c1ccc2ccccc2c1", "c1ccoc1"] {
        let mol = parse(smiles).unwrap();
        let config = TorsionKnowledgeConfig {
            use_exp_torsions: true,
            ..TorsionKnowledgeConfig::default()
        };
        let report = build_torsion_knowledge(&mol, &config);
        for pot in &report.potentials {
            assert_ne!(
                pot.source,
                TorsionKnowledgeSource::StandardExperimental,
                "{smiles}: StandardExperimental must never fire on an aromatic-only molecule"
            );
        }
    }
}

/// Spec §13: "coordinates changing when all flags are false" must FAIL to
/// occur. Since an empty `TorsionKnowledgeReport` has no potentials, there
/// is nothing for `optimize_torsions` to act on -- the coordinate-level
/// no-op follows directly from the report-level no-op, verified here on a
/// realistic multi-functional-group molecule (not just a toy case).
#[test]
fn disabled_flags_report_is_empty_on_a_realistic_molecule() {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
    let report = build_torsion_knowledge(&mol, &TorsionKnowledgeConfig::default());
    assert!(report.potentials.is_empty());
    let adjustments =
        macrocycle_14_bound_adjustments(&mol, &TorsionKnowledgeConfig::default()).unwrap();
    assert!(adjustments.is_empty());

    // Coordinates a Coordinator integration would produce with no
    // potentials to optimize are, trivially, the raw embedder's own output
    // -- confirmed unchanged by construction (there is no optimize_torsions
    // call to make when `report.potentials` is empty).
    let params = EmbedParameters::default();
    let coords_a = embed_distance_geometry_v2(&mol, &params).unwrap();
    let coords_b = embed_distance_geometry_v2(&mol, &params).unwrap();
    assert_eq!(
        coords_a, coords_b,
        "same-input reproducibility (sanity check)"
    );
}

// ---------------------------------------------------------------------------
// "Negative controls for the harness itself": deliberately corrupt an
// expected torsion minimum, ring size, or rule id and confirm the relevant
// check actually fails. Demonstrated via `catch_unwind` around an
// intentionally-wrong assertion -- proof the check has teeth, executed as
// part of a test suite that itself must still pass overall.
// ---------------------------------------------------------------------------

/// A corrupted ring-size expectation must make the check fail.
#[test]
fn harness_self_test_corrupted_ring_size_is_caught() {
    let mol = parse("C1CCCCC1").unwrap(); // cyclohexane: real ring size is 6
    let rings = RingMembershipIndex::build(&mol);
    let classification = classify_bond(&mol, &rings, AtomIdx(0), AtomIdx(1));

    // The TRUE assertion (must not panic): confirms the real behavior.
    assert_eq!(classification.ring, RingMembership::SmallRing(6));

    // The DELIBERATELY WRONG assertion: must panic, proving the check
    // mechanism can actually fail when the expectation is wrong (a check
    // that can't fail isn't a check).
    let result = std::panic::catch_unwind(|| {
        assert_eq!(classification.ring.clone(), RingMembership::SmallRing(7));
    });
    assert!(
        result.is_err(),
        "a deliberately wrong ring-size expectation (7 instead of 6) must panic, not pass silently"
    );
}

/// A corrupted rule-id expectation must make the check fail.
#[test]
fn harness_self_test_corrupted_rule_id_is_caught() {
    let mol = parse("CC(=O)NC").unwrap(); // N-methylacetamide
    let config = TorsionKnowledgeConfig {
        use_exp_torsions: true,
        ..TorsionKnowledgeConfig::default()
    };
    let report = build_torsion_knowledge(&mol, &config);
    assert!(
        !report.matched_rule_ids.is_empty(),
        "expected at least one matched rule"
    );

    let real_id = report.matched_rule_ids[0].clone();
    // TRUE assertion: the real id is present.
    assert!(report.matched_rule_ids.contains(&real_id));

    // DELIBERATELY WRONG assertion: a rule id that cannot possibly be
    // present must fail the check, not pass silently.
    let result = std::panic::catch_unwind(|| {
        assert!(
            report
                .matched_rule_ids
                .contains(&"standard:this_rule_id_does_not_exist".to_string()),
            "corrupted rule id lookup"
        );
    });
    assert!(
        result.is_err(),
        "a fabricated rule id must fail the containment check"
    );
}

/// A corrupted expected torsion-preference minimum must make the check
/// fail. Uses the legacy `score_torsion` API (angle + penalty model) since
/// it directly exposes a single numeric "minimum" to corrupt.
#[test]
fn harness_self_test_corrupted_torsion_minimum_is_caught() {
    use chematic_3d::etkdg_knowledge::{TorsionPreference, score_torsion};

    let pref = TorsionPreference {
        angle_deg: 180.0,
        penalty_per_degree: 0.1,
    };
    // TRUE assertion: scoring at the real minimum gives ~0 penalty.
    assert!(score_torsion(180.0, &pref).abs() < 1e-9);

    // DELIBERATELY WRONG assertion: claiming 90 degrees (clearly off the
    // real 180-degree minimum) also scores ~0 must fail.
    let result = std::panic::catch_unwind(|| {
        assert!(
            score_torsion(90.0, &pref).abs() < 1e-9,
            "wrong claimed minimum"
        );
    });
    assert!(
        result.is_err(),
        "a deliberately wrong torsion-minimum expectation (90 instead of 180) must fail the score check"
    );
}
