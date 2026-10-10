use crate::{canonical_smiles, parse, rdkit_canonical_smiles};
use chematic_core::Chirality;

#[test]
fn tetrahedral_output_preserves_each_recorded_hydrogen_slot() {
    for source in [
        "F[C@H](Cl)Br",
        "[C@@H](F)(Cl)Br",
        "N[C@H](C)O",
        "F[C@H]1CCCC1Cl",
        "[2H][C@](F)(Cl)Br",
    ] {
        let mol = parse(source).unwrap();
        let expected = rdkit_canonical_smiles(&mol).unwrap();
        let expected_native = canonical_smiles(&mol);
        for (center, atom) in mol.atoms().filter(|(_, a)| a.chirality.is_tetrahedral()) {
            let original = mol.stereo_neighbor_order(center).unwrap();
            assert_eq!(original.len(), 4);
            for a in 0..4 {
                for b in 0..4 {
                    if a == b {
                        continue;
                    }
                    for c in 0..4 {
                        if c == a || c == b {
                            continue;
                        }
                        let d = (0..4).find(|&d| d != a && d != b && d != c).unwrap();
                        let indices = [a, b, c, d];
                        let inversions = (0..4)
                            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                            .filter(|&(i, j)| indices[i] > indices[j])
                            .count();
                        let mut changed = mol.clone();
                        changed.set_stereo_neighbor_order(
                            center,
                            indices.iter().map(|&i| original[i]).collect(),
                        );
                        let clockwise =
                            (atom.chirality == Chirality::Clockwise) != (inversions % 2 == 1);
                        changed.set_chirality(
                            center,
                            if clockwise {
                                Chirality::Clockwise
                            } else {
                                Chirality::CounterClockwise
                            },
                        );
                        assert_eq!(
                            rdkit_canonical_smiles(&changed).unwrap(),
                            expected,
                            "{source}: {indices:?}"
                        );
                        assert_eq!(
                            canonical_smiles(&changed),
                            expected_native,
                            "{source}: {indices:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn query_logical_combinations_serialize_like_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-serialization-geometry-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["queries"].as_array().unwrap() {
        let query = row["query"].as_str().unwrap();
        let result = crate::rdkit_smarts_to_smarts(query);
        match result {
            Ok(actual) if actual == row["canonical"].as_str().unwrap() => {}
            other => failures.push(format!("{query}: {other:?} != {}", row["canonical"])),
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
#[test]
fn census_depictions_match_rdkit_or_the_documented_coordinate_residuals() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-serialization-geometry-boundary.json"
    ))
    .unwrap();
    let mut failures = Vec::new();
    for row in fixture["geometries"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = crate::parse(source).unwrap();
        let actual = crate::rdkit_2d_coords(&mol).unwrap();
        let expected = row["coords"].as_array().unwrap();
        assert_eq!(actual.len(), expected.len(), "{source}");
        assert!(actual.iter().flatten().all(|v| v.is_finite()), "{source}");
        assert_eq!(actual, crate::rdkit_2d_coords(&mol).unwrap(), "{source}");
        if actual.iter().enumerate().any(|(i, point)| {
            (0..2).any(|axis| (point[axis] - expected[i][axis].as_f64().unwrap()).abs() > 1e-9)
        }) {
            failures.push(source.to_owned());
        }
    }
    // These five complex layouts differ from the pinned reference. Keep the
    // complete sample and an explicit residual gate so new differences fail.
    let known_residuals = [
        "CC(C)C[C@H](O)[C@H](O)[C@H](CC1CCCCC1)NC(=O)C(NC(=O)[C@H](Cc1ccccc1)NS(=O)(=O)N1CCOCC1)OCC(F)(F)F",
        "O=C(OCc1ccccc1)C(Cc1ccc(C(F)(F)P(=O)(O)O)cc1)(Cc1ccc(C(F)(F)P(=O)(O)O)cc1)C(=O)OCc1ccccc1",
        r"COC1/C=C/OC2(C)Oc3c(C)c(O)c4c(O)c(c(/C=N/N5CCN(C)CC5)c(O)c4c3C2=O)NC(=O)/C(C)=C\C=C\C(C)C(O)C(C)C(O)C(C)C(OC(C)=O)C1C",
        "COC(=O)c1ccc(C(Cc2ccc(C(F)(F)P(=O)(O)O)cc2)(Cc2ccc(C(F)(F)P(=O)(O)O)cc2)n2nnc3ccccc32)cc1",
        "CN(Cc1ccc(S(=O)(=O)c2ccccc2)cc1)c1ccc2c3c(cccc13)C(N)=N2",
    ];
    assert_eq!(failures, known_residuals);
}

#[test]
fn conjugated_directional_smiles_respelling_matches_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-canonical-direction-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let actual = crate::parse(source)
            .map_err(|e| e.to_string())
            .and_then(|mol| crate::rdkit_canonical_smiles(&mol).map_err(|e| e.to_string()));
        if actual.as_ref().ok().map(String::as_str) != row["canonical"].as_str() {
            failures.push(format!("{source}: {actual:?} != {}", row["canonical"]));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
