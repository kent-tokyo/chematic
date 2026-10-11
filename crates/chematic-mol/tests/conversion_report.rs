use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::{
    ConversionLimits, ConversionPolicy as Policy, PreservationStatus as Status,
    convert_with_report, convert_with_report_json, convert_with_report_with_limits,
};

const SOURCE: &str = include_str!("../../../validation/mol2_preservation/amide.mol2");

#[test]
fn strict_refuses_without_molecule_output_but_supplies_loss_reasons() {
    let report = convert_with_report(SOURCE, "mol2", "smiles", Policy::Strict).unwrap();
    assert_eq!(report.status, "refused");
    assert!(report.output.is_none());
    assert!(report.candidate_graph_verified);
    for field in [
        "coordinates",
        "partial_charges",
        "atom_types",
        "atom_identity",
        "substructures",
        "opaque_sections",
        "bond_records",
        "molecule_metadata",
        "source_text",
    ] {
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == field && o.status == Status::Dropped),
            "missing {field}"
        );
    }
    for field in ["formal_charge_interpretation", "stereochemistry"] {
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == field && o.status == Status::Unverified)
        );
    }
}

#[test]
fn allow_loss_keeps_neutral_amide_and_exposes_source_correspondence() {
    let report = convert_with_report(SOURCE, "mol2", "smi", Policy::AllowLoss).unwrap();
    assert_eq!(report.status, "completed_with_loss");
    assert!(report.candidate_graph_verified);
    let mol = chematic_smiles::parse(report.output.as_ref().unwrap()).unwrap();
    assert!(mol.atoms().all(|(_, atom)| atom.charge == 0));
    let index = |id| {
        AtomIdx(
            report
                .candidate_atom_mapping
                .iter()
                .find(|m| m.source_id == id)
                .unwrap()
                .candidate_index as u32,
        )
    };
    assert_eq!(
        mol.bond_between(index(7), index(90)).unwrap().1.order,
        BondOrder::Double
    );
    assert_eq!(
        mol.bond_between(index(7), index(3)).unwrap().1.order,
        BondOrder::Single
    );
    assert_eq!(report.candidate_bond_mapping.len(), 3);
}

#[test]
fn same_format_all_policies_retain_exact_bytes_without_claiming_interpretation() {
    for policy in [Policy::Strict, Policy::AllowLoss, Policy::PreserveSource] {
        let input = SOURCE.replace('\n', "\r\n");
        let report = convert_with_report(&input, " .MOL2 ", ".mol2", policy).unwrap();
        assert_eq!(report.output.as_deref(), Some(input.as_str()));
        assert_eq!(report.status, "completed");
        assert_eq!(report.verification_scope, "exact_source_bytes");
        assert!(!report.candidate_graph_verified);
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.status == Status::OpaqueRetained)
        );
        assert!(report.candidate_atom_mapping.is_empty());
    }
}

#[test]
fn unsupported_routes_policies_and_invalid_source_are_errors() {
    let limits = ConversionLimits::default();
    for (input, from, to, policy) in [
        (SOURCE, "sdf", "smiles", "strict"),
        (SOURCE, "mol2", "cml", "allow_loss"),
        (SOURCE, "mol2", "smiles", "preserve_source"),
        (SOURCE, "mol2", "smiles", "typo"),
        ("invalid", "mol2", "mol2", "strict"),
    ] {
        assert!(convert_with_report_json(input, from, to, policy, &limits).is_err());
    }
    assert!(
        convert_with_report(
            &SOURCE.replace("ACM -0.70", "ACM NaN"),
            "mol2",
            "mol2",
            Policy::PreserveSource
        )
        .is_err()
    );
}

#[test]
fn both_input_and_escaped_report_sizes_are_bounded() {
    let mut limits = ConversionLimits::default();
    limits.mol2.max_input_bytes = SOURCE.len() - 1;
    assert!(
        convert_with_report_with_limits(SOURCE, "mol2", "mol2", Policy::Strict, &limits).is_err()
    );
    limits.mol2.max_input_bytes = SOURCE.len();
    limits.max_output_bytes = SOURCE.len();
    // Source output fits exactly, JSON envelope/escaping does not.
    assert!(
        convert_with_report_with_limits(SOURCE, "mol2", "mol2", Policy::Strict, &limits).is_ok()
    );
    assert!(convert_with_report_json(SOURCE, "mol2", "mol2", "strict", &limits).is_err());
    limits.max_output_bytes = SOURCE.len() - 1;
    assert!(
        convert_with_report_with_limits(SOURCE, "mol2", "mol2", Policy::Strict, &limits).is_err()
    );
}

#[test]
fn missing_charge_column_does_not_invent_partial_charge_loss() {
    let source = "@<TRIPOS>MOLECULE\natom\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n8 C1 0 0 0 C.3\n";
    let report = convert_with_report(source, "mol2", "smiles", Policy::AllowLoss).unwrap();
    assert!(
        !report
            .observations
            .iter()
            .any(|o| o.field == "partial_charges"
                || o.field == "opaque_sections"
                || o.field == "substructures")
    );
}

#[test]
fn zero_bond_projection_change_is_detected_not_reported_preserved() {
    let source = SOURCE.replace("12 7 3 am", "12 7 3 nc");
    let report = convert_with_report(&source, "mol2", "smiles", Policy::AllowLoss).unwrap();
    assert!(!report.candidate_graph_verified);
    assert!(report.candidate_atom_mapping.is_empty());
    assert!(
        report
            .observations
            .iter()
            .any(|o| o.field == "projected_graph" && o.status == Status::Unverified)
    );
}

#[test]
fn branch_ring_and_disconnected_correspondence_covers_every_atom_and_bond() {
    for smiles in [
        "CC(O)(F)N",
        "C1CCOC1.CN",
        "[H]C([H])([H])[H]",
        "c1ccncc1",
        "[N+](C)(C)(C)C",
    ] {
        let graph = chematic_smiles::parse(smiles).unwrap();
        let source = chematic_mol::write_mol2(&graph, &[]);
        let report = convert_with_report(&source, "mol2", "smiles", Policy::AllowLoss).unwrap();
        assert!(
            report.candidate_graph_verified,
            "{smiles}: {:?}",
            report.output
        );
        assert_eq!(report.candidate_atom_mapping.len(), graph.atom_count());
        assert_eq!(report.candidate_bond_mapping.len(), graph.bond_count());
    }
}

#[test]
fn recursive_projection_has_an_independent_atom_budget() {
    let mut source =
        String::from("@<TRIPOS>MOLECULE\nlarge\n513 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n");
    for id in 1..=513 {
        source.push_str(&format!("{id} C{id} 0 0 0 C.3\n"));
    }
    assert!(
        convert_with_report(&source, "mol2", "smiles", Policy::AllowLoss)
            .unwrap_err()
            .to_string()
            .contains("512")
    );
    assert!(convert_with_report(&source, "mol2", "mol2", Policy::Strict).is_ok());
}

#[test]
fn cjson_retains_xyz_geometry_and_projected_formal_charges_without_partial_charge_rounding() {
    let source = SOURCE.replace("0.123456789 0 0 C.3", "0.123456789 2.5 -3.75 N.4");
    let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
    assert_eq!(report.status, "completed_with_loss");
    assert!(report.candidate_graph_verified);
    let (mol, coords) = chematic_mol::parse_cjson(report.output.as_ref().unwrap()).unwrap();
    assert_eq!(coords[0], (0.123456789, 2.5, -3.75));
    assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
    assert_eq!(mol.atom(AtomIdx(2)).charge, 0); // source partial charge -0.70
    for field in ["coordinates", "projected_formal_charges"] {
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == field && o.status == Status::Preserved),
            "{field}"
        );
    }
    assert!(
        report
            .observations
            .iter()
            .any(|o| o.field == "partial_charges" && o.status == Status::Preserved)
    );
    assert!(
        !report
            .observations
            .iter()
            .any(|o| o.field == "implicit_hydrogens" && o.status == Status::Inferred)
    );
    assert_eq!(report.candidate_atom_mapping[0].source_id, 50);
    assert_eq!(report.candidate_atom_mapping[0].candidate_index, 0);
    assert_eq!(report.candidate_bond_mapping.len(), 3);
}

#[test]
fn cjson_strict_still_refuses_other_source_losses_and_preserve_source_is_same_format_only() {
    let report = convert_with_report(SOURCE, "mol2", "cjson", Policy::Strict).unwrap();
    assert_eq!(report.status, "refused");
    assert!(report.output.is_none());
    assert!(
        report
            .observations
            .iter()
            .any(|o| o.field == "coordinates" && o.status == Status::Preserved)
    );
    assert!(convert_with_report(SOURCE, "mol2", "cjson", Policy::PreserveSource).is_err());
}

#[test]
fn cjson_coordinate_check_covers_f64_precision_and_signed_zero() {
    for literal in [
        "-0.0",
        "-1.4175225513503917",
        "1.2345678901234567",
        "0.12345678901234567",
        "2.2250738585072014e-308",
        "5e-324",
        "1.7976931348623157e308",
    ] {
        let source = SOURCE.replace("0.123456789", literal);
        let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
        let (_, coordinates) = chematic_mol::parse_cjson(report.output.as_ref().unwrap()).unwrap();
        let actual = coordinates[0].0;
        let expected: f64 = literal.parse().unwrap();
        assert_eq!(actual.to_bits(), expected.to_bits(), "{literal}");
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == "coordinates" && o.status == Status::Preserved),
            "{literal}"
        );
    }
}

#[test]
fn cjson_zero_and_query_bond_changes_do_not_claim_graph_preservation() {
    for token in ["nc", "du", "un"] {
        let source = SOURCE.replace("12 7 3 am", &format!("12 7 3 {token}"));
        let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
        assert!(!report.candidate_graph_verified, "{token}");
        assert!(report.candidate_atom_mapping.is_empty());
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == "projected_graph" && o.status == Status::Unverified)
        );
    }
}

#[test]
fn cjson_uses_reader_atom_budget_instead_of_recursive_smiles_budget() {
    let mut source =
        String::from("@<TRIPOS>MOLECULE\nlarge\n513 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n");
    for id in 1..=513 {
        source.push_str(&format!("{id} C{id} 0 0 1 C.3\n"));
    }
    let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
    assert_eq!(report.candidate_atom_mapping.len(), 513);
    assert!(report.candidate_graph_verified);
    let mut limits = ConversionLimits::default();
    limits.mol2.max_atoms = 512;
    assert!(
        convert_with_report_with_limits(&source, "mol2", "cjson", Policy::AllowLoss, &limits)
            .is_err()
    );
}

#[test]
fn cjson_partial_charges_keep_method_name_and_f64_bits_separate_from_formal_charges() {
    for literal in [
        "-0.0",
        "-1.4175225513503917",
        "5e-324",
        "1.7976931348623157e308",
    ] {
        let source = SOURCE.replace("ACM 0.70", &format!("ACM {literal}"));
        let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
        let json: serde_json::Value =
            serde_json::from_str(report.output.as_ref().unwrap()).unwrap();
        assert_eq!(
            json["partialCharges"]["USER_CHARGES"][0]
                .as_f64()
                .unwrap()
                .to_bits(),
            literal.parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(json["properties"]["name"], "acetamide");
        let (mol, _) = chematic_mol::parse_cjson(report.output.as_ref().unwrap()).unwrap();
        assert_eq!(mol.atom(AtomIdx(0)).charge, 0);
        for field in ["partial_charges", "charge_method", "molecule_name"] {
            assert!(
                report
                    .observations
                    .iter()
                    .any(|o| o.field == field && o.status == Status::Preserved)
            );
        }
    }
}

#[test]
fn cjson_does_not_fill_partial_charges_or_claim_undeclared_charge_methods() {
    for source in [
        SOURCE.replace("ACM 0.70", "ACM"),
        SOURCE.replace("USER_CHARGES", "NO_CHARGES"),
        SOURCE.replace("USER_CHARGES", ""),
    ] {
        let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
        let json: serde_json::Value =
            serde_json::from_str(report.output.as_ref().unwrap()).unwrap();
        assert!(json.get("partialCharges").is_none());
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == "partial_charges" && o.status == Status::Dropped)
        );
        assert!(
            !report
                .observations
                .iter()
                .any(|o| o.field == "charge_method" && o.status == Status::Preserved)
        );
    }
    let source = SOURCE.replace("USER_CHARGES", "DECLARED/custom~method");
    let report = convert_with_report(&source, "mol2", "cjson", Policy::AllowLoss).unwrap();
    let json: serde_json::Value = serde_json::from_str(report.output.as_ref().unwrap()).unwrap();
    assert_eq!(json["partialCharges"]["DECLARED/custom~method"][0], 0.70);
}
