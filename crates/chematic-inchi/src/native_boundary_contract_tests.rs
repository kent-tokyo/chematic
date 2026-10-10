use crate::{dedup::*, rdkit_mol_from_inchi};
use chematic_core::MoleculeBuilder;
use chematic_smiles::{parse, rdkit_canonical_smiles};
use serde_json::Value;
#[test]
fn native_inchi_boundary_reconstruction_matches_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-inchi-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let input = row["inchi"].as_str().unwrap();
        let actual = rdkit_mol_from_inchi(input).and_then(|m| {
            rdkit_canonical_smiles(&m)
                .map_err(|e| crate::InchiError::InvalidInput(format!("{e:?}")))
        });
        if let Some(expected) = row["smiles"].as_str() {
            match actual {
                Ok(got) if got == expected => (),
                other => failures.push(format!(
                    "{}: {other:?}; expected {expected}",
                    row["source_smiles"]
                )),
            }
        } else if actual.is_ok() {
            failures.push(format!(
                "{}: expected an invalid reconstruction",
                row["source_smiles"]
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn accurate_preflight_batch_reconciles_stereo_and_isotope_policies() {
    let mut mols = [
        "CCO",
        "OCC",
        "CCN",
        "C[C@H](O)F",
        "C[C@@H](O)F",
        "[13CH3]CO",
    ]
    .iter()
    .map(|s| parse(s).unwrap())
    .collect::<Vec<_>>();
    mols.push(MoleculeBuilder::new().build());
    for policy in [
        IdentityPolicy::StandardInchiString,
        IdentityPolicy::StandardInchiKey,
        IdentityPolicy::StereoIgnored,
        IdentityPolicy::IsotopeIgnored,
    ] {
        let report = deduplicate_verified_with_accurate_cip_preflight(&mols, policy);
        assert_eq!(report.invalid_molecules, vec![6]);
        assert!(report.verification_unavailable.is_empty());
        let groups = report
            .groups
            .iter()
            .map(|g| g.members.clone())
            .collect::<Vec<_>>();
        let expected = match policy {
            IdentityPolicy::StereoIgnored => vec![vec![0, 1], vec![3, 4]],
            IdentityPolicy::IsotopeIgnored => vec![vec![0, 1, 5]],
            _ => vec![vec![0, 1]],
        };
        assert_eq!(groups, expected, "{policy:?}");
        assert!(report.canonical_collisions.is_empty());
        let split_count = if matches!(
            policy,
            IdentityPolicy::StereoIgnored | IdentityPolicy::IsotopeIgnored
        ) {
            1
        } else {
            0
        };
        assert_eq!(report.canonical_splits.len(), split_count);
        assert_eq!(
            compare_with_accurate_cip_preflight("same", &mols[0], "same", &mols[2], policy),
            DedupRelation::CanonicalCollision
        );
        assert_eq!(
            compare_with_accurate_cip_preflight("first", &mols[0], "second", &mols[1], policy),
            DedupRelation::CanonicalSplit
        );
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(&mols[0], &mols[1], policy),
            DedupRelation::VerifiedDuplicate
        );
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(&mols[0], &mols[2], policy),
            DedupRelation::Distinct
        );
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(&mols[0], &mols[6], policy),
            DedupRelation::InvalidMolecule
        );
    }
}
