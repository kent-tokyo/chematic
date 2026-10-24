use super::{
    BatchLimits, MAX_OUTPUT_BYTES, MAX_SMILES_ATOMS, MAX_SMILES_INPUT_BYTES,
    batch_canonicalize_json, batch_descriptors_json, batch_fingerprints_json, batch_reactions_json,
    batch_report_json, batch_similarity_json, batch_standardize_json, batch_substructure_json,
    convert_text, cxsmiles_json, descriptors_json, fingerprint_json, parse_cli_smiles, parse_json,
    reaction_balance_json, reaction_fingerprint_json, reaction_json, reaction_match_json,
    reaction_similarity_json, report_json, similarity_json, standardize_json, substructure_json,
};

fn default_batch_limits() -> BatchLimits {
    BatchLimits::default()
}

#[test]
fn converts_smiles_to_mol2_and_back() {
    let mol2 = convert_text("CCO", "smiles", "mol2").unwrap();
    assert!(mol2.contains("@<TRIPOS>MOLECULE"));
    let smiles = convert_text(&mol2, ".mol2", "smi").unwrap();
    let reparsed = chematic_smiles::parse(&smiles).unwrap();
    let expected = chematic_smiles::parse("CCO").unwrap();
    assert_eq!(reparsed.atom_count(), expected.atom_count());
    assert_eq!(reparsed.bond_count(), expected.bond_count());
}

#[test]
fn converts_mol2_to_mol2_without_dropping_interchange_fields() {
    let input = "@<TRIPOS>MOLECULE\nethanol\n 3 2 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG1 0.125\n2 C2 1.5 0 0 C.3 1 LIG1 0.250\n3 O1 3 0 0 O.3 1 LIG1 -0.375\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n";
    let rewritten = convert_text(input, "mol2", "mol2").unwrap();
    let record = chematic_mol::parse_mol2_record(&rewritten).unwrap();
    assert_eq!(record.atoms[0].atom_type, "C.3");
    assert_eq!(record.atoms[2].partial_charge, Some(-0.375));
    assert_eq!(record.atoms[0].subst_name.as_deref(), Some("LIG1"));
}

#[test]
fn mol2_self_conversion_retains_source_metadata_and_geometry() {
    let source = include_str!("../../../validation/mol2_preservation/amide.mol2");
    assert_eq!(convert_text(source, ".mol2", "mol2").unwrap(), source);
    let invalid = source.replace("ACM -0.70", "ACM NaN");
    assert!(convert_text(&invalid, "mol2", "mol2").is_err());
}

#[test]
fn cxsmiles_json_preserves_mapped_wildcard_attachment_identity() {
    let value: serde_json::Value =
        serde_json::from_str(&cxsmiles_json("[*:11][C:12] |$_AP7;$|").unwrap()).unwrap();
    assert_eq!(value["atom_maps"], serde_json::json!([11, 12]));
    assert_eq!(value["atom_labels"], serde_json::json!(["_AP7", null]));
    assert_eq!(
        value["marked_attachment_points"],
        serde_json::json!([{"atom_index": 0, "atom_map": 11, "identifier": 7}])
    );
    assert_eq!(value["cxsmiles"], "[*:11][C:12] |$_AP7;$|");
}

#[test]
fn converts_v3000_fixture_to_a_structurally_equivalent_smiles() {
    let block = include_str!("../../../benchmarks/fixtures/ethanol.v3000");
    let smiles = convert_text(block, "mol_v3000", "smiles").unwrap();
    let converted = chematic_smiles::parse(&smiles).unwrap();
    assert_eq!(converted.atom_count(), 3);
    assert_eq!(converted.bond_count(), 2);
    assert_eq!(converted.total_formula(), "C2H6O");
}

#[test]
fn v3000_to_v3000_keeps_layout_coordinates() {
    let block = include_str!("../../../benchmarks/fixtures/ethanol.v3000");
    let rewritten = convert_text(block, "mol_v3000", "mol_v3000").unwrap();
    assert!(
        rewritten.contains("M  V30 2 C 1.5000 0.0000 0.0000 0"),
        "same-format V3000 conversion must not collapse coordinates: {rewritten}"
    );
}

#[test]
fn refuses_to_flatten_v3000_query_attributes_into_smiles() {
    let block = "query\n  chematic\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 SUBST=3\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n";
    let error = convert_text(block, "mol_v3000", "smiles").unwrap_err();
    assert!(error.contains("unsupported V3000 query semantics"));
    let roundtrip = convert_text(block, "mol_v3000", "mol_v3000").unwrap();
    assert!(roundtrip.contains("SUBST=3"));
}

#[test]
fn rejects_unknown_formats() {
    assert!(
        convert_text("CCO", "smiles", "nope")
            .unwrap_err()
            .contains("unsupported output format")
    );
}

#[test]
fn descriptors_are_json_and_include_core_fields() {
    let json: serde_json::Value = serde_json::from_str(&descriptors_json("CCO").unwrap()).unwrap();
    assert_eq!(json["heavy_atoms"], 3);
    assert!(json["molecular_weight"].as_f64().unwrap() > 40.0);
    assert!(json["rdkit_molecular_weight"].as_f64().unwrap() > 40.0);
    assert!(json["molar_refractivity"].as_f64().unwrap() > 10.0);
    assert_eq!(json["fsp3"], 1.0);
    assert_eq!(json["aromatic_ring_count"], 0);
    assert_eq!(json["formula"], "C2H6O");
}

#[test]
fn parse_reports_structure_without_expensive_descriptors() {
    let json: serde_json::Value = serde_json::from_str(&parse_json("C[NH3+]").unwrap()).unwrap();
    assert_eq!(json["input_smiles"], "C[NH3+]");
    assert_eq!(json["canonical_smiles"], "C[NH3+]");
    assert_eq!(json["atoms"], 2);
    assert_eq!(json["bonds"], 1);
    assert_eq!(json["formal_charge"], 1);
}

#[test]
fn batch_canonicalize_retains_input_order_and_rejections_without_descriptors() {
    let result: serde_json::Value = serde_json::from_str(
        &batch_canonicalize_json("CCO\nC1CC\nCCN\n", &BatchLimits::default()).unwrap(),
    )
    .unwrap();
    assert_eq!(result["operation"], "canonicalize");
    assert_eq!(result["valid_count"], 2);
    assert_eq!(result["error_count"], 1);
    assert_eq!(result["input_count"], 3);
    assert_eq!(result["accepted_count"], 2);
    assert_eq!(result["rejected_count"], 1);
    assert_eq!(result["refused_count"], 0);
    assert_eq!(result["skipped_count"], 0);
    assert_eq!(result["all_succeeded"], false);
    assert_eq!(result["records"][0]["input_index"], 0);
    assert_eq!(result["records"][0]["status"], "accepted");
    assert_eq!(result["records"][0]["canonical_smiles"], "C(C)O");
    assert!(result["records"][1]["canonical_smiles"].is_null());
    assert_eq!(result["records"][1]["status"], "rejected");
    assert_eq!(result["records"][1]["error_stage"], "parse");
    assert!(result["records"][1]["error"].as_str().is_some());
    assert_eq!(result["records"][2]["input_index"], 2);
    assert_eq!(result["records"][2]["canonical_smiles"], "C(C)N");
}

#[test]
fn parse_rejects_invalid_smiles() {
    assert!(parse_json("C1CC").is_err());
}

#[test]
fn fingerprint_json_is_stable_and_reports_set_bits() {
    let json: serde_json::Value =
        serde_json::from_str(&fingerprint_json("CCO", "ecfp4").unwrap()).unwrap();
    assert_eq!(json["algorithm"], "ecfp4");
    assert_eq!(json["n_bits"], 2048);
    assert!(json["popcount"].as_u64().unwrap() > 0);
    assert_eq!(
        json["set_bits"].as_array().unwrap().len(),
        json["popcount"].as_u64().unwrap() as usize
    );
}

#[test]
fn fingerprint_rejects_unknown_algorithms() {
    assert!(
        fingerprint_json("CCO", "unknown")
            .unwrap_err()
            .contains("unsupported fingerprint algorithm")
    );
}

#[test]
fn fingerprint_supports_rdkit_rdk_algorithm() {
    let json: serde_json::Value =
        serde_json::from_str(&fingerprint_json("CCO", "rdkit_rdk").unwrap()).unwrap();
    assert_eq!(json["algorithm"], "rdkit_rdk");
    assert_eq!(json["n_bits"], 2048);
    assert!(json["popcount"].as_u64().unwrap() > 0);
}

#[test]
fn fingerprint_supports_rdkit_torsion_algorithm() {
    let json: serde_json::Value =
        serde_json::from_str(&fingerprint_json("CCO", "rdkit_torsion").unwrap()).unwrap();
    assert_eq!(json["algorithm"], "rdkit_torsion");
    assert_eq!(json["n_bits"], 2048);
}

#[test]
fn similarity_is_one_for_identical_molecules() {
    let json: serde_json::Value =
        serde_json::from_str(&similarity_json("CCO", "CCO", "ecfp4").unwrap()).unwrap();
    assert_eq!(json["similarity"], 1.0);
    assert_eq!(json["algorithm"], "ecfp4");
}

#[test]
fn similarity_rejects_unknown_algorithms() {
    assert!(
        similarity_json("CCO", "CCN", "unknown")
            .unwrap_err()
            .contains("unsupported fingerprint algorithm")
    );
}

#[test]
fn substructure_reports_sorted_atom_mappings() {
    let json: serde_json::Value =
        serde_json::from_str(&substructure_json("c1ccccc1", "c", false).unwrap()).unwrap();
    assert_eq!(json["match_count"], 6);
    assert_eq!(
        json["matches"].as_array().unwrap()[0]
            .as_array()
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn substructure_rejects_invalid_smarts() {
    assert!(
        substructure_json("CCO", "[", false)
            .unwrap_err()
            .contains("SMARTS")
    );
}

#[test]
fn standardize_reports_pipeline_and_canonical_output() {
    let json: serde_json::Value =
        serde_json::from_str(&standardize_json("C[NH3+]", false).unwrap()).unwrap();
    assert_eq!(json["input_smiles"], "C[NH3+]");
    assert_eq!(json["status"], "modified");
    assert!(json["changed"].as_bool().unwrap());
    assert_eq!(json["steps"].as_array().unwrap().len(), 6);
    assert!(
        json["steps"]
            .as_array()
            .unwrap()
            .iter()
            .any(|step| step["step"] == "neutralize_charges" && step["changed"] == true)
    );
    assert!(json["warnings"].is_array());
}

#[test]
fn standardize_fragment_profile_is_explicit_and_default_preserves_fragments() {
    let default: serde_json::Value =
        serde_json::from_str(&standardize_json("CC.CCC", false).unwrap()).unwrap();
    let profile: serde_json::Value =
        serde_json::from_str(&standardize_json("CC.CCC", true).unwrap()).unwrap();
    assert!(default["output_smiles"].as_str().unwrap().contains('.'));
    assert!(!profile["output_smiles"].as_str().unwrap().contains('.'));
    assert_ne!(default["output"]["atoms"], profile["output"]["atoms"]);
}

#[test]
fn report_emits_complete_machine_readable_analysis() {
    let json: serde_json::Value =
        serde_json::from_str(&report_json("CC(=O)Oc1ccccc1C(=O)O").unwrap()).unwrap();
    assert_eq!(json["input_smiles"], "CC(=O)Oc1ccccc1C(=O)O");
    assert!(json["canonical_smiles"].as_str().is_some());
    assert!(json["descriptors"]["molecular_weight"].as_f64().unwrap() > 100.0);
    assert!(json["filters"]["pains_alerts"].is_array());
    assert!(json["functional_groups"].is_array());
}

#[test]
fn report_rejects_invalid_smiles() {
    assert!(report_json("C1CC").is_err());
}

#[test]
fn reaction_normalizes_components_and_reports_mapped_center() {
    let json: serde_json::Value =
        serde_json::from_str(&reaction_json("[CH3:1][OH:2]>>[CH3:1][O-:2]").unwrap()).unwrap();
    assert_eq!(json["reaction_smiles"], "[OH:2][CH3:1]>>[CH3:1][O-:2]");
    assert_eq!(json["reactants"], 1);
    assert_eq!(json["products"], 1);
    assert_eq!(json["mapped"], true);
    assert_eq!(
        json["reaction_center"]["changed_atoms"],
        serde_json::json!([1])
    );
}

#[test]
fn reaction_rejects_missing_arrow() {
    assert!(
        reaction_json("CCO")
            .unwrap_err()
            .contains("reaction SMILES")
    );
}

#[test]
fn reaction_match_reports_boolean_without_inference() {
    let json: serde_json::Value =
        serde_json::from_str(&reaction_match_json("CCO>>CCO", "[#6]>>[#6]").unwrap()).unwrap();
    assert_eq!(json["matched"], true);
    assert_eq!(json["query"], "[#6]>>[#6]");
}

#[test]
fn reaction_match_rejects_invalid_query() {
    assert!(
        reaction_match_json("CCO>>CCO", "[#6]")
            .unwrap_err()
            .contains("invalid reaction SMARTS")
    );
}

#[test]
fn reaction_balance_reports_counts_and_differences() {
    let balanced: serde_json::Value =
        serde_json::from_str(&reaction_balance_json("CO.CO>>COC.O").unwrap()).unwrap();
    assert_eq!(balanced["balanced"], true);
    assert!(balanced["differences"].as_array().unwrap().is_empty());

    let unbalanced: serde_json::Value =
        serde_json::from_str(&reaction_balance_json("C>>CC").unwrap()).unwrap();
    assert_eq!(unbalanced["balanced"], false);
    assert!(
        unbalanced["differences"]
            .as_array()
            .unwrap()
            .iter()
            .any(|difference| difference.as_str().unwrap().contains("C"))
    );
}

#[test]
fn reaction_fingerprint_reports_transformation_bits() {
    let json: serde_json::Value =
        serde_json::from_str(&reaction_fingerprint_json("CCO>>CC=O", "xor").unwrap()).unwrap();
    assert_eq!(json["mode"], "xor");
    assert_eq!(json["n_bits"], 2048);
    assert!(json["popcount"].as_u64().unwrap() > 0);
    assert_eq!(
        json["set_bits"].as_array().unwrap().len(),
        json["popcount"].as_u64().unwrap() as usize
    );
}

#[test]
fn reaction_fingerprint_rejects_unknown_mode() {
    assert!(
        reaction_fingerprint_json("CCO>>CCO", "bad")
            .unwrap_err()
            .contains("unsupported reaction fingerprint mode")
    );
}

#[test]
fn reaction_similarity_is_one_for_identical_reactions() {
    let json: serde_json::Value =
        serde_json::from_str(&reaction_similarity_json("CCO>>CC=O", "CCO>>CC=O").unwrap()).unwrap();
    assert_eq!(json["similarity"], 1.0);
    assert_eq!(json["fingerprint"], "reaction_ecfp4_xor");
}

#[test]
fn reaction_similarity_rejects_invalid_reaction() {
    assert!(reaction_similarity_json("CCO", "CCO>>CCO").is_err());
}

#[test]
fn batch_report_retains_partial_errors_in_input_order() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_report_json("CCO\nC1CC\n# comment\nCCN\n", &default_batch_limits()).unwrap(),
    )
    .unwrap();
    assert_eq!(json["schema_version"], 1);
    assert_eq!(json["operation"], "report");
    assert_eq!(json["status"], "complete");
    assert_eq!(json["record_count"], 3);
    let records = json["records"].as_array().unwrap();
    assert_eq!(records.len(), 3);
    assert_eq!(records[0]["input_index"], 0);
    assert!(records[0]["report"].is_object());
    assert!(records[1]["error"].as_str().is_some());
    assert_eq!(records[2]["input_index"], 2);
    assert!(records[2]["report"].is_object());
}

#[test]
fn batch_descriptors_returns_lightweight_partial_manifest() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_descriptors_json("CCO\nC1CC\n# comment\nCCN\n", &default_batch_limits()).unwrap(),
    )
    .unwrap();
    assert_eq!(json["valid_count"], 2);
    assert_eq!(json["error_count"], 1);
    assert_eq!(json["schema_version"], 1);
    assert_eq!(json["operation"], "descriptors");
    assert_eq!(json["status"], "complete");
    assert_eq!(json["record_count"], 3);
    assert_eq!(
        json["limits"]["max_records"],
        default_batch_limits().max_records
    );
    let records = json["records"].as_array().unwrap();
    assert!(records[0]["descriptors"].is_object());
    assert!(records[1]["error"].as_str().is_some());
    assert!(records[2]["descriptors"].is_object());
}

#[test]
fn batch_fingerprints_returns_algorithm_and_partial_manifest() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_fingerprints_json("CCO\nC1CC\nCCN\n", "ecfp4", &default_batch_limits()).unwrap(),
    )
    .unwrap();
    assert_eq!(json["algorithm"], "ecfp4");
    assert_eq!(json["valid_count"], 2);
    assert_eq!(json["error_count"], 1);
    assert!(json["records"][0]["fingerprint"]["set_bits"].is_array());
    assert!(json["records"][1]["error"].as_str().is_some());
}

#[test]
fn batch_fingerprints_rejects_unknown_algorithm_per_record() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_fingerprints_json("CCO\nCCN\n", "bad", &default_batch_limits()).unwrap(),
    )
    .unwrap();
    assert_eq!(json["valid_count"], 0);
    assert_eq!(json["error_count"], 2);
}

#[test]
fn batch_standardize_retains_audit_reports_and_errors() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_standardize_json("C[NH3+]\nC1CC\nCCO\n", false, &default_batch_limits()).unwrap(),
    )
    .unwrap();
    assert_eq!(json["valid_count"], 2);
    assert_eq!(json["error_count"], 1);
    assert!(json["records"][0]["standardization"]["steps"].is_array());
    assert!(json["records"][1]["error"].as_str().is_some());
    assert_eq!(json["records"][2]["standardization"]["status"], "unchanged");
}

#[test]
fn batch_similarity_processes_tsv_pairs_and_retains_errors() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_similarity_json(
            "CCO\tCCO\nCCO\tC1CC\nmissing\n",
            "ecfp4",
            &default_batch_limits(),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(json["algorithm"], "ecfp4");
    assert_eq!(json["valid_count"], 1);
    assert_eq!(json["error_count"], 2);
    assert_eq!(json["records"][0]["similarity"]["similarity"], 1.0);
    assert!(json["records"][1]["error"].as_str().is_some());
    assert!(json["records"][2]["error"].as_str().is_some());
}

#[test]
fn batch_substructure_processes_tsv_queries_and_retains_errors() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_substructure_json("c1ccccc1\tc\nCCO\t[\n", &default_batch_limits(), false).unwrap(),
    )
    .unwrap();
    assert_eq!(json["valid_count"], 1);
    assert_eq!(json["error_count"], 1);
    assert_eq!(json["records"][0]["substructure"]["match_count"], 6);
    assert!(json["records"][1]["error"].as_str().is_some());
}

#[test]
fn batch_reactions_returns_normalized_records_and_errors() {
    let json: serde_json::Value = serde_json::from_str(
        &batch_reactions_json("CCO>>CCO\nCCO\n", &default_batch_limits()).unwrap(),
    )
    .unwrap();
    assert_eq!(json["valid_count"], 1);
    assert_eq!(json["error_count"], 1);
    assert_eq!(json["records"][0]["reaction"]["reactants"], 1);
    assert!(json["records"][1]["error"].as_str().is_some());
}

#[test]
fn batch_limits_reject_oversized_input_and_records() {
    let mut limits = default_batch_limits();
    limits.max_input_bytes = 3;
    assert!(
        batch_reactions_json("CCO\n", &limits)
            .unwrap_err()
            .contains("max-input-bytes")
    );

    limits = default_batch_limits();
    limits.max_records = 1;
    assert!(
        batch_reactions_json("CCO>>CCO\nCCN>>CCN\n", &limits)
            .unwrap_err()
            .contains("max-records")
    );

    limits = default_batch_limits();
    limits.max_line_bytes = 3;
    assert!(
        batch_reactions_json("CCO>>CCO\n", &limits)
            .unwrap_err()
            .contains("max-line-bytes")
    );
}

#[test]
fn single_smiles_contract_rejects_oversized_bytes_before_parse() {
    let input = "C".repeat(MAX_SMILES_INPUT_BYTES + 1);
    let error = parse_cli_smiles(&input)
        .err()
        .expect("oversized input must fail");
    assert!(error.contains("maximum input size"));
}

#[test]
fn single_smiles_contract_rejects_oversized_molecules() {
    let input = "C".repeat(MAX_SMILES_ATOMS + 1);
    let error = parse_cli_smiles(&input)
        .err()
        .expect("oversized molecule must fail");
    assert!(error.contains("maximum atom count"));
}

#[test]
fn output_limit_is_explicit_and_bounded() {
    let oversized = "x".repeat(MAX_OUTPUT_BYTES + 1);
    assert!(
        super::write_output(Some(&std::path::PathBuf::from("/dev/null")), &oversized)
            .unwrap_err()
            .contains("maximum size")
    );
}

#[test]
fn gjf_common_route_retains_source_and_refuses_cross_format_projection() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gjf-cartesian-2026-10-22/source-cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        for (input, output) in [("gjf", "gjf"), (".COM", ".GJF"), ("gjf", "com")] {
            let result = convert_text(source, input, output);
            if case["expected"] == "retain" {
                assert_eq!(result.unwrap().as_bytes(), source.as_bytes());
            } else {
                assert!(result.is_err(), "{}", case["id"]);
            }
        }
    }
    assert!(
        convert_text("C", "smiles", "gjf")
            .unwrap_err()
            .contains("same-format")
    );
    assert!(
        convert_text("irrelevant", "gjf", "smiles")
            .unwrap_err()
            .contains("same-format")
    );
    let source = format!("# HF\n\nlimit\n\n0 1\n{}\n", "H 0 0 0\n".repeat(10_001));
    assert!(convert_text(&source, "gjf", "com").is_err());
}
#[test]
fn gjf_document_checks_every_job_without_changing_single_job_conversion() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gjf-document-2026-10-22/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        let result = crate::conversion_format::gjf_document(source, false, 128);
        if case["expected"] == "retain" {
            assert_eq!(result.unwrap(), source);
            let inspected: serde_json::Value = serde_json::from_str(
                &crate::conversion_format::gjf_document(source, true, 128).unwrap(),
            )
            .unwrap();
            assert_eq!(
                inspected["jobs"].as_array().unwrap().len(),
                case["jobs"].as_array().unwrap().len()
            );
            if case["jobs"].as_array().unwrap().len() > 1 {
                assert!(convert_text(source, "gjf", "gjf").is_err());
            }
        } else {
            assert!(result.is_err(), "{}", case["id"]);
            assert!(crate::conversion_format::gjf_document(source, true, 128).is_err());
        }
    }
    let source = include_str!("../../../validation/gjf-document-2026-10-22/two-jobs.gjf");
    for cap in [0, 1, 129] {
        assert!(crate::conversion_format::gjf_document(source, false, cap).is_err());
    }
}

#[test]
fn moljson_hydrogen_authority_survives_cli_conversion() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/moljson_hydrogens/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let input = case["input"].as_str().unwrap();
        let expected = chematic_smiles::canonical_smiles(
            &chematic_smiles::parse(case["expected_smiles"].as_str().unwrap()).unwrap(),
        );
        assert_eq!(
            convert_text(input, "moljson", "smiles").unwrap(),
            expected,
            "{}",
            case["id"]
        );
        let output = convert_text(input, "moljson", "moljson").unwrap();
        let value: serde_json::Value = serde_json::from_str(&output).unwrap();
        let counts: Vec<_> = value["atoms"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| a["hydrogens"].clone())
            .collect();
        assert_eq!(serde_json::json!(counts), case["expected_hydrogens"]);
    }
    for hydrogens in packet["invalid_hydrogens"].as_array().unwrap() {
        let input = serde_json::json!({"atoms":[{"id":"a1","element":"C","hydrogens":hydrogens}],"bonds":[]}).to_string();
        assert!(
            convert_text(&input, "moljson", "smiles")
                .unwrap_err()
                .to_string()
                .contains("hydrogens")
        );
    }
}
