#[cfg(test)]
mod invert_stereocenter_tests {
    use super::*;
    use chematic_core::{Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder};

    #[test]
    fn wasm_api_synchronizes_chirality_and_wedge() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(Atom::new(Element::C));
        let f = builder.add_atom(Atom::new(Element::F));
        let cl = builder.add_atom(Atom::new(Element::CL));
        let br = builder.add_atom(Atom::new(Element::BR));
        let i = builder.add_atom(Atom::new(Element::I));
        let wedge = builder.add_bond(center, f, BondOrder::Up).unwrap();
        builder.add_bond(center, cl, BondOrder::Single).unwrap();
        builder.add_bond(center, br, BondOrder::Single).unwrap();
        builder.add_bond(center, i, BondOrder::Single).unwrap();
        builder.set_stereo_neighbor_order(center, vec![f.0, cl.0, br.0, i.0]);
        let mut molecule = builder.build();
        molecule.set_chirality(center, Chirality::Clockwise);
        let handle = MolHandle {
            inner: std::rc::Rc::new(molecule),
        };

        let inverted = invert_stereocenter_at(&handle, center.0).unwrap();
        assert_eq!(
            inverted.inner.atom(AtomIdx(center.0)).chirality,
            Chirality::CounterClockwise
        );
        assert_eq!(inverted.inner.bond(wedge).order, BondOrder::Down);
    }
}

// ---------------------------------------------------------------------------
// mol_transforms (3D geometry manipulation)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod mcs_config_tests {
    use super::*;

    #[test]
    fn default_config_matches_bare_mcs_smiles_json() {
        let smiles = r#"["CC(=O)Oc1ccccc1C(=O)O","CC(=O)Nc1ccc(O)cc1"]"#;
        let full = mcs_smiles_json_with_config(smiles, "{}").unwrap();
        let bare = mcs_smiles_json(smiles).unwrap();
        let value: serde_json::Value = serde_json::from_str(&full).unwrap();
        assert_eq!(value["smiles"], serde_json::Value::String(bare));
        assert_eq!(value["wasTimedOut"], false);
        assert_eq!(value["wasResourceLimited"], false);
    }

    #[test]
    fn concrete_mcs_projection_uses_first_input_alternative() {
        for (first, second, config) in [
            ("NO", "CO", r#"{"atomCompare":"any"}"#),
            ("CO", "NO", r#"{"atomCompare":"any"}"#),
            ("C=N", "CN", r#"{"bondCompare":"any"}"#),
            ("CN", "C=N", r#"{"bondCompare":"any"}"#),
        ] {
            let inputs = serde_json::to_string(&[first, second]).unwrap();
            let result = mcs_smiles_json_with_config(&inputs, config).unwrap();
            let value: serde_json::Value = serde_json::from_str(&result).unwrap();
            let expected = chematic_smiles::canonical_smiles(&chematic_smiles::parse(first).unwrap());
            assert_eq!(value["smiles"], expected, "{first} versus {second}");
        }
    }

    #[test]
    fn no_common_substructure_returns_null_smiles() {
        // Two disconnected, chemically unrelated single-heavy-atom molecules.
        let smiles = r#"["[He]","[Ne]"]"#;
        let result = mcs_smiles_json_with_config(smiles, "{}").unwrap();
        let value: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(value["smiles"], serde_json::Value::Null);
        assert_eq!(value["wasTimedOut"], false);
        assert_eq!(value["wasResourceLimited"], false);
    }

    #[test]
    fn match_charge_true_prevents_charged_neutral_match() {
        // Acetate vs acetic acid: MCS with match_charge=false matches the
        // full 4-heavy-atom core; match_charge=true must reject the
        // charge-differing carboxylate oxygen, shrinking the MCS.
        let smiles = r#"["CC(=O)[O-]","CC(=O)O"]"#;
        let without = mcs_smiles_json_with_config(smiles, r#"{"matchCharge":false}"#).unwrap();
        let with = mcs_smiles_json_with_config(smiles, r#"{"matchCharge":true}"#).unwrap();
        assert_ne!(
            without, with,
            "match_charge=true must change the MCS result"
        );
    }

    #[test]
    fn atom_compare_any_heavy_atom_widens_match() {
        // Benzene vs pyridine: default (element-exact) MCS excludes the N
        // position; any_heavy_atom compare should include it (6-atom ring).
        let smiles = r#"["c1ccccc1","c1ccncc1"]"#;
        let elements =
            mcs_smiles_json_with_config(smiles, r#"{"atomCompare":"elements"}"#).unwrap();
        let any_heavy =
            mcs_smiles_json_with_config(smiles, r#"{"atomCompare":"any_heavy_atom"}"#).unwrap();
        assert_ne!(
            elements, any_heavy,
            "any_heavy_atom compare must change the MCS result vs elements"
        );
    }

    #[test]
    fn rdkit_heavy_comparator_keeps_explicit_hydrogen_pairs() {
        let smiles = r#"["[2H]O[2H]","[H]OC"]"#;
        let legacy = mcs_smiles_json_with_config(
            smiles,
            r#"{"atomCompare":"any_heavy_atom"}"#,
        )
        .unwrap();
        let rdkit = mcs_smiles_json_with_config(
            smiles,
            r#"{"atomCompare":"rdkit_any_heavy_atom"}"#,
        )
        .unwrap();
        let legacy: serde_json::Value = serde_json::from_str(&legacy).unwrap();
        let rdkit: serde_json::Value = serde_json::from_str(&rdkit).unwrap();
        let legacy_smiles = legacy["smiles"].as_str().unwrap();
        let rdkit_smiles = rdkit["smiles"].as_str().unwrap();
        assert_eq!(chematic_smiles::parse(legacy_smiles).unwrap().atom_count(), 1);
        assert_eq!(chematic_smiles::parse(rdkit_smiles).unwrap().atom_count(), 2);
        assert_eq!(rdkit["wasTimedOut"], false);
        assert_eq!(rdkit["wasResourceLimited"], false);
    }

    #[test]
    fn timeout_zero_is_reported_as_timed_out() {
        let smiles = r#"["CC(=O)Oc1ccccc1C(=O)O","CC(=O)Nc1ccc(O)cc1"]"#;
        let result = mcs_smiles_json_with_config(smiles, r#"{"timeoutMs":0}"#).unwrap();
        let value: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(value["wasTimedOut"], true);
        assert_eq!(value["wasResourceLimited"], false);
    }

    #[test]
    fn search_node_limit_is_reported_separately_from_timeout() {
        let smiles = r#"["CCC","CCC"]"#;
        let result =
            mcs_smiles_json_with_config(smiles, r#"{"maxSearchNodes":0}"#).unwrap();
        let value: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert_eq!(value["wasTimedOut"], false);
        assert_eq!(value["wasResourceLimited"], true);
    }

    // Error paths that surface as a real thrown `JsValue` are only safely
    // constructible on the wasm32 target (constructing a `JsValue` natively
    // aborts the process, per wasm-bindgen's own design) -- covered instead
    // by `tests/mcs.test.mjs`, run under the actual wasm runtime via CI's
    // `Test (WASM)` job. What's safe and worth testing natively here is the
    // pure-serde config-parsing layer, which never touches `JsValue`.

    #[test]
    fn invalid_atom_compare_string_rejected_at_parse_time() {
        let err = serde_json::from_str::<McsConfigJson>(r#"{"atomCompare":"bogus"}"#).unwrap_err();
        assert!(err.to_string().contains("atomCompare") || err.to_string().contains("unknown"));
    }

    #[test]
    fn unknown_config_field_rejected_at_parse_time() {
        let err = serde_json::from_str::<McsConfigJson>(r#"{"notAField":true}"#).unwrap_err();
        assert!(err.to_string().contains("notAField") || err.to_string().contains("unknown"));
    }
}
