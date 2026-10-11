//! Public reaction APIs exercised over the checked-in compatibility fixtures.

use chematic_rxn::{
    RdkitDetailedProfileOutcome, ReactionTransformLimits,
    run_reactants_traced_rdkit_2026_03_6_detailed, run_reactants_traced_with_rejection_diagnostics,
};
use chematic_smiles::parse;

fn fixture_cases() -> Vec<serde_json::Value> {
    let base: serde_json::Value = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/reaction_product_parity_cases.json"
    )))
    .unwrap();
    let strata: serde_json::Value = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/reaction_product_parity_strata_v2.json"
    )))
    .unwrap();
    let supplement: serde_json::Value = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/reaction_product_new_maps_v3.json"
    )))
    .unwrap();

    base.as_array()
        .unwrap()
        .iter()
        .chain(strata["cases"].as_array().unwrap())
        .chain(supplement["cases"].as_array().unwrap())
        .cloned()
        .collect()
}

fn assert_accounting(
    accepted: usize,
    applied: usize,
    rejected: usize,
    accepted_sets: usize,
    rejected_sets: usize,
) {
    assert_eq!(accepted, applied + rejected);
    assert_eq!(applied, accepted_sets);
    assert_eq!(rejected, rejected_sets);
}

#[test]
fn detailed_native_and_rdkit_profiles_account_for_every_fixture_match() {
    let limits = ReactionTransformLimits {
        max_matches: 20_000,
    };
    let mut exercised = 0usize;

    for case in fixture_cases() {
        let smirks = case["smirks"].as_str().unwrap();
        let Some(reactant_texts) = case["reactants"].as_array() else {
            continue;
        };
        let Ok(reactants) = reactant_texts
            .iter()
            .map(|value| parse(value.as_str().unwrap()))
            .collect::<Result<Vec<_>, _>>()
        else {
            continue;
        };
        let refs = reactants.iter().collect::<Vec<_>>();

        if let Ok(report) = run_reactants_traced_with_rejection_diagnostics(smirks, &refs, &limits)
        {
            let d = report.report.diagnostics;
            assert_accounting(
                d.accepted_matches,
                d.applied_products,
                d.valence_rejected_matches,
                report.report.products.len(),
                report.rejection_diagnostics.len(),
            );
            for rejected_set in report.rejection_diagnostics {
                assert!(rejected_set.iter().any(|product| !product.accepted));
                assert!(
                    rejected_set
                        .iter()
                        .all(|product| { product.accepted == product.reason.is_none() })
                );
            }
        }

        if let Ok(RdkitDetailedProfileOutcome::Report(report)) =
            run_reactants_traced_rdkit_2026_03_6_detailed(smirks, &refs, &limits)
        {
            let d = report.report.diagnostics;
            assert_accounting(
                d.accepted_matches,
                d.applied_products,
                d.valence_rejected_matches,
                report.report.products.len(),
                report.rejection_diagnostics.len(),
            );
            assert_eq!(
                report.rejected_products.len(),
                report.rejection_diagnostics.len()
            );
        }
        exercised += 1;
    }

    assert!(exercised >= 83, "only exercised {exercised} fixtures");
}

#[test]
fn rdkit_profile_stereo_matrix_is_deterministic_and_bounded() {
    let reactants = [
        "N[C@@H](C)C(=O)O",
        "N[C@H](C)C(=O)O",
        "F[C@H](Cl)Br",
        "F[C@@H](Cl)Br",
        "C1[C@H](F)CCC1",
        "C1[C@@H](F)CCC1",
    ];
    let templates = [
        "[N:1][C@@H:2]([C:3])[C:4](=[O:5])[O:6]>>[N:1][C@@H:2]([C:3])[C:4](=[O:5])[O:6]C",
        "[N:1][C@@H:2]([C:3])[C:4](=[O:5])[O:6]>>[N:1][C@H:2]([C:3])[C:4](=[O:5])[O:6]",
        "[F:1][C@H:2]([Cl:3])[Br:4]>>[F:1][C@@H:2]([Cl:3])[Br:4]",
        "[F:1][C@@H:2]([Cl:3])[Br:4]>>[F:1][C@@H:2]([Cl:3])[Br:4]",
        "[C:1]1[C@H:2]([F:3])[C:4][C:5][C:6]1>>[C:1]1[C@@H:2]([F:3])[C:4][C:5][C:6]1",
    ];
    let limits = ReactionTransformLimits { max_matches: 1_000 };

    for smiles in reactants {
        let mol = parse(smiles).unwrap();
        for smirks in templates {
            let first = run_reactants_traced_rdkit_2026_03_6_detailed(smirks, &[&mol], &limits);
            let second = run_reactants_traced_rdkit_2026_03_6_detailed(smirks, &[&mol], &limits);
            assert_eq!(first.is_ok(), second.is_ok());
            if let (
                Ok(RdkitDetailedProfileOutcome::Report(a)),
                Ok(RdkitDetailedProfileOutcome::Report(b)),
            ) = (first, second)
            {
                assert_eq!(a.report.products.len(), b.report.products.len());
                assert_eq!(a.rejected_products.len(), b.rejected_products.len());
                assert_eq!(a.report.diagnostics, b.report.diagnostics);
            }
        }
    }
}
