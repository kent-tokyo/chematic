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
fn census_depictions_preserve_reference_bond_lengths_and_are_deterministic() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-serialization-geometry-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    for (case, row) in fixture["geometries"].as_array().unwrap().iter().enumerate() {
        let source = row["smiles"].as_str().unwrap();
        let mol = crate::parse(source).unwrap();
        let actual = crate::rdkit_2d_coords(&mol).unwrap();
        let expected: Vec<[f64; 2]> = serde_json::from_value(row["coords"].clone()).unwrap();
        assert_eq!(actual.len(), expected.len(), "{source}");
        assert!(actual.iter().flatten().all(|v| v.is_finite()), "{source}");
        assert_eq!(actual, crate::rdkit_2d_coords(&mol).unwrap(), "{source}");
        let parsed = crate::rdkit_parsed_molecule(&mol).unwrap();
        assert_eq!(parsed.atom_count(), actual.len(), "{source}");
        // Complete layouts can differ across CPU architectures at floating-
        // point tie breaks. Bond lengths are invariant under rigid motion and
        // independent of alternative fragment placement.
        for bi in 0..parsed.bond_count() {
            let bond = parsed.bond(chematic_core::BondIdx(bi as u32));
            let a = bond.atom1.0 as usize;
            let b = bond.atom2.0 as usize;
            let length = |points: &[[f64; 2]]| {
                (points[a][0] - points[b][0]).hypot(points[a][1] - points[b][1])
            };
            assert!(
                (length(&actual) - length(&expected)).abs() < 1e-6,
                "bond {bi} in {source}: {} != {}",
                length(&actual),
                length(&expected)
            );
        }
        // The eight small hand-picked layouts have stable exact references.
        if case < 8 {
            for (point, reference) in actual.iter().zip(&expected) {
                for (value, expected) in point.iter().zip(reference) {
                    assert!((value - expected).abs() < 1e-9, "{source}");
                }
            }
        }
    }
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
