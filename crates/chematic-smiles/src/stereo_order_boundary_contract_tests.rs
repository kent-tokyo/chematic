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
            let actual_length = length(&actual);
            let expected_length = length(&expected);
            // Collision repair may shrink a terminal bond by 0.9 (depict.rs). This
            // one observed platform tie selects either the 1.35 or 1.5 layout.
            let known_collision_scale = source
                == "CCCCc1cn(-c2c(C(C)C)cccc2C(C)C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1"
                && [actual_length, expected_length]
                    .iter()
                    .all(|v| (*v - 1.35).abs() < 1e-6 || (*v - 1.5).abs() < 1e-6);
            assert!(
                (actual_length - expected_length).abs() < 1e-6 || known_collision_scale,
                "bond {bi} in {source}: {} != {}",
                actual_length,
                expected_length
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

#[test]
fn disconnected_rooted_smarts_and_cx_extensions_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-smarts-components-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = if source.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            parse(source).unwrap()
        };
        let root = row["root"].as_u64().map(|v| v as usize);
        let result = crate::rdkit_smarts(&mol, row["isomeric"].as_bool().unwrap(), root);
        match result {
            Ok(actual) if actual == row["smarts"].as_str().unwrap() => {}
            other => failures.push(format!(
                "{source},root={root:?},isomeric={}: {other:?} != {}",
                row["isomeric"], row["smarts"]
            )),
        }
    }
    for row in fixture["cx"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = if source.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            parse(source).unwrap()
        };
        let result = crate::rdkit_cx_smarts(&mol);
        match result {
            Ok(actual) if actual == row["cx_smarts"].as_str().unwrap() => {}
            other => failures.push(format!("CX {source}: {other:?} != {}", row["cx_smarts"])),
        }
        if mol.atom_count() > 0 {
            let error = crate::rdkit_smarts(&mol, true, Some(mol.atom_count())).unwrap_err();
            assert!(error.to_string().contains("bad atom index"));
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 196);
    assert_eq!(fixture["cx"].as_array().unwrap().len(), 17);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn branched_polyene_direction_respelling_matches_native_reader_cleanup() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-branched-directions-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = parse(source).unwrap();
        let result = rdkit_canonical_smiles(&mol);
        match result {
            Ok(actual) if actual == row["canonical"].as_str().unwrap() => {}
            other => failures.push(format!("{source}: {other:?} != {}", row["canonical"])),
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 4698);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
