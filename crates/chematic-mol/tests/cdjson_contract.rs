use chematic_core::{AtomIdx, implicit_hcount};
use chematic_mol::{
    CdjsonError, CdjsonParseLimits, parse_cdjson, parse_cdjson_with_limits, write_cdjson_graph,
};
use serde_json::json;

#[test]
fn omitted_carbon_label_and_hydrogen_use_native_inference() {
    let d = diagnostic(parse_cdjson(r#"{"m":[{"a":[{"x":0,"y":0}]}]}"#)).unwrap();
    let mol = &d.molecules()[0].molecule;
    assert_eq!(mol.atom(AtomIdx(0)).element.atomic_number(), 6);
    assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, None);
    assert_eq!(implicit_hcount(mol, AtomIdx(0)), 4);
}

#[test]
fn explicit_zero_hydrogens_are_distinct_from_default_and_minus_one() {
    for (h, expected) in [(0, 0), (2, 2), (4, 4), (-1, 4)] {
        let text = json!({"a":[{"x":0,"y":0,"h":h}]}).to_string();
        let d = diagnostic(parse_cdjson(&text)).unwrap();
        assert_eq!(
            implicit_hcount(&d.molecules()[0].molecule, AtomIdx(0)),
            expected
        );
    }
}

#[test]
fn source_ids_charge_isotope_and_coordinates_are_preserved() {
    let source = "{\r\n\"m\":[{\"i\":\"m1\",\"a\":[{\"i\":\"a1\",\"l\":\"N\",\"c\":1,\"m\":15,\"h\":4,\"x\":1.2345678901234567,\"y\":-2,\"z\":3}]}]}\r\n";
    let d = diagnostic(parse_cdjson(source)).unwrap();
    assert_eq!(d.source_text(), source);
    let m = &d.molecules()[0];
    assert_eq!(m.id.as_deref(), Some("m1"));
    assert_eq!(m.atom_ids, [Some("a1".into())]);
    assert_eq!(m.coordinates[0], [1.2345678901234567, -2., 3.]);
    let atom = m.molecule.atom(AtomIdx(0));
    assert_eq!(
        (atom.charge, atom.isotope, atom.hydrogen_count),
        (1, Some(15), Some(4))
    );
}

#[test]
fn all_content_molecules_are_returned_and_empty_content_is_valid() {
    let d = diagnostic(parse_cdjson(
        r#"{"m":[{"a":[]},{"a":[{"x":0,"y":0,"l":"Cl"}]}],"s":[]}"#,
    ))
    .unwrap();
    assert_eq!(d.molecules().len(), 2);
    assert_eq!(d.molecules()[0].molecule.atom_count(), 0);
    assert_eq!(d.molecules()[1].molecule.atom_count(), 1);
    assert!(
        diagnostic(parse_cdjson("{}"))
            .unwrap()
            .molecules()
            .is_empty()
    );
}

#[test]
fn bond_indices_defaults_and_ids_are_checked() {
    let d = diagnostic(parse_cdjson(
        r#"{"a":[{"x":0,"y":0},{"x":20,"y":0,"l":"O"}],"b":[{"i":"b7","b":0,"e":1,"s":"none"}]}"#,
    ))
    .unwrap();
    let m = &d.molecules()[0];
    assert_eq!(m.bond_ids, [Some("b7".into())]);
    assert_eq!(implicit_hcount(&m.molecule, AtomIdx(0)), 3);
    assert_eq!(implicit_hcount(&m.molecule, AtomIdx(1)), 1);
    for bond in [
        json!({"b":0,"e":0}),
        json!({"b":0,"e":2}),
        json!({"b":-1,"e":1}),
        json!({"b":0.0,"e":1}),
    ] {
        let input = json!({"a":[{"x":0,"y":0},{"x":1,"y":1}],"b":[bond]}).to_string();
        assert!(diagnostic(parse_cdjson(&input)).is_err());
    }
    assert!(
        diagnostic(parse_cdjson(
            r#"{"a":[{"x":0,"y":0},{"x":1,"y":1}],"b":[{"b":0,"e":1},{"b":1,"e":0}]}"#
        ))
        .is_err()
    );
}

#[test]
fn unsupported_chemistry_and_unknown_fields_never_disappear_silently() {
    for field in [
        json!({"q":{}}),
        json!({"s2":{"t":"abs"}}),
        json!({"r":1}),
        json!({"p":1}),
        json!({"custom":7}),
        json!({"l":"8"}),
    ] {
        let mut atom = json!({"x":0,"y":0});
        atom.as_object_mut()
            .unwrap()
            .extend(field.as_object().unwrap().clone());
        assert!(matches!(
            diagnostic(parse_cdjson(&json!({"a":[atom]}).to_string())),
            Err(CdjsonError::Unsupported(_))
        ));
    }
    assert!(diagnostic(parse_cdjson(r#"{"m":[],"s":[{"t":"Line"}]}"#)).is_err());
    for order in [0., 0.5, 1.5, 2.5, 4.] {
        assert!(
            diagnostic(parse_cdjson(
                &json!({"a":[{"x":0,"y":0},{"x":1,"y":1}],"b":[{"b":0,"e":1,"o":order}]})
                    .to_string()
            ))
            .is_err()
        );
    }
}

#[test]
fn malformed_atoms_and_duplicate_decoded_json_keys_are_rejected() {
    for source in [
        r#"{"a":[{"x":0,"x":1,"y":0}]}"#,
        r#"{"a":[{"x":0,"\u0078":1,"y":0}]}"#,
        r#"{"a":[{"x":0,"y":0,"h":null}]}"#,
        r#"{"a":[{"x":0}]}"#,
        r#"{"a":[{"x":1e999,"y":0}]}"#,
    ] {
        assert!(diagnostic(parse_cdjson(source)).is_err(), "{source}");
    }
    for (key, value) in [
        ("c", json!(128)),
        ("c", json!(1.0)),
        ("h", json!(-2)),
        ("h", json!(256)),
        ("m", json!(0)),
        ("m", json!(65536)),
    ] {
        let mut atom = json!({"x":0,"y":0});
        atom[key] = value;
        assert!(diagnostic(parse_cdjson(&json!({"a":[atom]}).to_string())).is_err());
    }
    assert!(
        diagnostic(parse_cdjson(
            r#"{"m":[{"i":"same","a":[{"i":"same","x":0,"y":0}]}]}"#
        ))
        .is_err()
    );
}

#[test]
fn caller_limits_cover_bytes_and_aggregate_counts() {
    let input =
        r#"{"m":[{"a":[{"x":0,"y":0},{"x":1,"y":0}],"b":[{"b":0,"e":1}]},{"a":[{"x":0,"y":0}]}]}"#;
    for limits in [
        CdjsonParseLimits {
            max_input_bytes: input.len() - 1,
            ..Default::default()
        },
        CdjsonParseLimits {
            max_molecules: 1,
            ..Default::default()
        },
        CdjsonParseLimits {
            max_atoms: 2,
            ..Default::default()
        },
        CdjsonParseLimits {
            max_bonds: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            diagnostic(parse_cdjson_with_limits(input, &limits)),
            Err(CdjsonError::ResourceLimit { .. })
        ));
    }
    let limits = CdjsonParseLimits {
        max_input_bytes: input.len(),
        max_molecules: 2,
        max_atoms: 3,
        max_bonds: 1,
    };
    assert_eq!(
        diagnostic(parse_cdjson_with_limits(input, &limits))
            .unwrap()
            .molecules()
            .len(),
        2
    );
}

#[test]
fn hard_caps_cannot_be_raised_and_are_checked_before_atom_projection() {
    let source = json!({"a":vec![json!({"x":0,"y":0});65537]}).to_string();
    let limits = CdjsonParseLimits {
        max_atoms: usize::MAX,
        ..Default::default()
    };
    assert!(matches!(
        diagnostic(parse_cdjson_with_limits(&source, &limits)),
        Err(CdjsonError::ResourceLimit {
            resource: "atoms",
            actual: 65537,
            limit: 65536
        })
    ));
    let source = json!({"m":vec![json!({"a":[]});1025]}).to_string();
    let limits = CdjsonParseLimits {
        max_molecules: usize::MAX,
        ..Default::default()
    };
    assert!(matches!(
        diagnostic(parse_cdjson_with_limits(&source, &limits)),
        Err(CdjsonError::ResourceLimit {
            resource: "molecules",
            limit: 1024,
            ..
        })
    ));
}

#[test]
fn exact_atom_and_molecule_caps_are_accepted() {
    let source = json!({"a":vec![json!({"x":0,"y":0});65536]}).to_string();
    let d = diagnostic(parse_cdjson(&source)).unwrap();
    assert_eq!(d.molecules()[0].molecule.atom_count(), 65536);
    let source = json!({"m":vec![json!({"a":[]});1024]}).to_string();
    assert_eq!(
        diagnostic(parse_cdjson(&source)).unwrap().molecules().len(),
        1024
    );
}

#[test]
fn graph_writer_roundtrips_elements_isotopes_charge_and_hydrogens() {
    for smiles in [
        "C",
        "CC",
        "C(=O)CO",
        "[C]",
        "[CH2]",
        "[NH4+]",
        "[Cl-]",
        "[13CH3]O",
        "[2H]O[2H]",
        "C1=CC=CC=C1",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let coordinates: Vec<_> = (0..mol.atom_count())
            .map(|i| [i as f64 * 1.2345678901234567, -2., 3.])
            .collect();
        let text = diagnostic(write_cdjson_graph(&mol, &coordinates)).unwrap();
        let d = diagnostic(parse_cdjson(&text)).unwrap();
        let recovered = &d.molecules()[0];
        assert_eq!(recovered.coordinates, coordinates);
        for ((idx, a), (_, b)) in mol.atoms().zip(recovered.molecule.atoms()) {
            assert_eq!(
                (a.element, a.charge, a.isotope),
                (b.element, b.charge, b.isotope)
            );
            assert_eq!(
                implicit_hcount(&mol, idx),
                implicit_hcount(&recovered.molecule, idx)
            );
        }
    }
}

#[test]
fn graph_writer_refuses_stereo_aromaticity_maps_and_bad_coordinates() {
    for smiles in [
        "c1ccccc1",
        "N->[Cu+2]",
        "F/C=C\\F",
        "[C@H](F)(Cl)Br",
        "[CH3:7]",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        assert!(
            diagnostic(write_cdjson_graph(&mol, &vec![[0.; 3]; mol.atom_count()])).is_err(),
            "{smiles}"
        );
    }
    let mol = chematic_smiles::parse("C").unwrap();
    assert!(diagnostic(write_cdjson_graph(&mol, &[])).is_err());
    assert!(diagnostic(write_cdjson_graph(&mol, &[[f64::NAN, 0., 0.]])).is_err());
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
