//! Exercise the executable boundary, including stdin, stdout, files and errors.
use serde_json::{Value, json};
use std::io::{ErrorKind, Write};
use std::process::{Command, Output, Stdio};

fn invoke(args: &[&str], stdin: &str) -> Output {
    invoke_bytes(args, stdin.as_bytes())
}

fn invoke_bytes(args: &[&str], stdin: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let write_result = child.stdin.take().unwrap().write_all(stdin);
    let output = child.wait_with_output().unwrap();
    match write_result {
        Ok(()) => {}
        // Invalid arguments can make the CLI exit before it reads stdin.
        Err(error) if error.kind() == ErrorKind::BrokenPipe && !output.status.success() => {}
        Err(error) => panic!("failed to send CLI input: {error}"),
    }
    output
}

#[test]
fn cli_global_input_guards_reject_invalid_bytes_and_oversized_molecules() {
    let excessive_limit = invoke(
        &["batch-descriptors", "--max-input-bytes", "67108865"],
        "CCO\n",
    );
    assert!(!excessive_limit.status.success());
    assert!(String::from_utf8_lossy(&excessive_limit.stderr).contains("CLI maximum"));

    let invalid_utf8 = invoke_bytes(&["batch-descriptors"], &[0xff, b'\n']);
    assert!(!invalid_utf8.status.success());
    assert!(String::from_utf8_lossy(&invalid_utf8.stderr).contains("not UTF-8"));

    let too_many_atoms = "C".repeat(10_001);
    let oversized_molecule = invoke(
        &[
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "moljson",
        ],
        &too_many_atoms,
    );
    assert!(!oversized_molecule.status.success());
    assert!(String::from_utf8_lossy(&oversized_molecule.stderr).contains("maximum atom count"));
}
fn json_output(args: &[&str], stdin: &str) -> Value {
    let out = invoke(args, stdin);
    assert!(
        out.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}

#[test]
fn rdkit_parity_search_perceives_aromaticity_in_single_and_batch_modes() {
    // RDKit 2025.09.3 sanitizes aromatic cyclobutadiene to a non-aromatic
    // ring, while a Kekule benzene is perceived as aromatic.
    let native = json_output(&["substructure", "c1ccc1", "c"], "");
    assert_eq!(native["match_count"], 4);
    let parity = json_output(&["substructure", "c1ccc1", "c", "--rdkit-parity"], "");
    assert_eq!(parity["match_count"], 0);
    assert_eq!(parity["rdkit_parity"], true);
    let arsenic_query = json_output(&["substructure", "s1cccc1", "[as]", "--rdkit-parity"], "");
    assert_eq!(arsenic_query["match_count"], 0);
    let opposite_alkene = json_output(
        &[
            "substructure",
            "F/C=C/F",
            "F/C=C\\F",
            "--rdkit-parity",
            "--use-chirality",
        ],
        "",
    );
    assert_eq!(opposite_alkene["match_count"], 0);
    let explicit_h = json_output(
        &["substructure", "[H][C@](F)(Cl)Br", "[D3]", "--rdkit-parity"],
        "",
    );
    assert_eq!(explicit_h["matches"], json!([[0]]));
    let ammonia = json_output(
        &["substructure", "[H]N([H])[H]", "[D3]", "--rdkit-parity"],
        "",
    );
    assert_eq!(ammonia["match_count"], 0);
    for parity in [false, true] {
        let mut args = vec!["substructure", "[*H]", "[H1]"];
        if parity {
            args.push("--rdkit-parity");
        }
        let wildcard_h = json_output(&args, "");
        assert_eq!(wildcard_h["match_count"], 1);
    }

    let batch = json_output(
        &["batch-substructure", "--rdkit-parity"],
        "c1ccc1\tc\nC1=CC=CC=C1\tc\n",
    );
    assert_eq!(batch["valid_count"], 2);
    assert_eq!(batch["records"][0]["substructure"]["match_count"], 0);
    assert_eq!(batch["records"][1]["substructure"]["match_count"], 6);
}

#[test]
fn rdkit_parity_search_uses_2026_09_relevant_cycle_ring_counts() {
    // RDKit 2026.09.1 has no atom in two relevant cycles in this bridged
    // hydrocarbon. The earlier symmetrized-SSSR profile reports an [R2] hit.
    let smiles = "c1cc2ccc1CCc1ccc(cc1)CC2";
    let single = json_output(&["substructure", smiles, "[R2]", "--rdkit-parity"], "");
    assert_eq!(single["match_count"], 0);
    let batch = json_output(
        &["batch-substructure", "--rdkit-parity"],
        &format!("{smiles}\t[R2]\n"),
    );
    assert_eq!(batch["records"][0]["substructure"]["match_count"], 0);
}

#[test]
fn conversion_round_trips_all_advertised_molecular_formats() {
    for format in ["v3000", "cml", "cjson", "moljson", "cdxml"] {
        let out = invoke(
            &[
                "convert",
                "--input-format",
                "smiles",
                "--output-format",
                format,
            ],
            "CCO",
        );
        assert!(
            out.status.success(),
            "{format}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let encoded = String::from_utf8(out.stdout).unwrap();
        let out = invoke(
            &[
                "convert",
                "--input-format",
                format,
                "--output-format",
                "smiles",
            ],
            &encoded,
        );
        assert!(
            out.status.success(),
            "{format}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let decoded = String::from_utf8(out.stdout).unwrap();
        let report = json_output(&["parse", decoded.trim()], "");
        assert_eq!(report["formula"], "C2H6O", "{format}");
        assert_eq!(report["bonds"], 2, "{format}");
    }
}

#[test]
fn converts_kekule_ctabs_to_aromatic_smiles() {
    for smiles in [
        "c1ccccc1",
        "c1cc[nH]c1",
        "[cH-]1cccc1",
        "[cH+]1cccccc1",
        "[cH+]1cc1",
        "C1=CC=CC1",
    ] {
        let source = chematic_smiles::parse(smiles).unwrap();
        let expected = chematic_smiles::canonical_smiles(&source);
        for format in ["mol", "mol_v3000"] {
            let block = invoke(
                &[
                    "convert",
                    "--input-format",
                    "smiles",
                    "--output-format",
                    format,
                ],
                smiles,
            );
            assert!(block.status.success(), "{format}: {smiles}");
            let converted = invoke(
                &[
                    "convert",
                    "--input-format",
                    format,
                    "--output-format",
                    "smiles",
                ],
                &String::from_utf8(block.stdout).unwrap(),
            );
            assert!(converted.status.success(), "{format}: {smiles}");
            assert_eq!(
                String::from_utf8(converted.stdout).unwrap().trim(),
                expected,
                "{format}: {smiles}"
            );
        }
    }
}

#[test]
fn preserves_cumulated_ring_and_alkene_stereo_from_ctabs() {
    let expected = "C1(Oc2ccc(O/C=C/C)cc2)=CC=C=C=C1";
    for (format, block) in [
        (
            "mol",
            include_str!("fixtures/cumulated_ring_stereo_v2000.mol"),
        ),
        (
            "mol_v3000",
            include_str!("fixtures/cumulated_ring_stereo_v3000.mol"),
        ),
    ] {
        let converted = invoke(
            &[
                "convert",
                "--input-format",
                format,
                "--output-format",
                "smiles",
            ],
            block,
        );
        assert!(converted.status.success(), "{format}");
        assert_eq!(
            String::from_utf8(converted.stdout).unwrap().trim(),
            expected
        );
    }
}

#[test]
fn molecule_and_reaction_commands_emit_inspectable_json() {
    for (args, path, expected) in [
        (vec!["parse", "CCO"], "/atoms", json!(3)),
        (vec!["descriptors", "CCO"], "/formula", json!("C2H6O")),
        (vec!["cxsmiles", "CCO"], "/atom_count", json!(3)),
        (
            vec!["fingerprint", "CCO", "--algorithm", "ecfp6"],
            "/algorithm",
            json!("ecfp6"),
        ),
        (vec!["similarity", "CCO", "CCO"], "/similarity", json!(1.)),
        (vec!["substructure", "CCO", "[O]"], "/matches", json!([[2]])),
        (
            vec!["standardize", "CCO.[Na+]", "--largest-fragment-only"],
            "/output/atoms",
            json!(3),
        ),
        (vec!["reaction", "CCO>>CC=O"], "/products", json!(1)),
        (
            vec!["reaction-match", "CCO>>CC=O", "CO>>C=O"],
            "/matched",
            json!(true),
        ),
        (
            vec!["reaction-balance", "CCO>>CCO"],
            "/balanced",
            json!(true),
        ),
        (
            vec!["reaction-fingerprint", "CCO>>CC=O", "--mode", "or"],
            "/mode",
            json!("or"),
        ),
        (
            vec!["reaction-similarity", "CCO>>CC=O", "CCO>>CC=O"],
            "/similarity",
            json!(1.),
        ),
    ] {
        let output = json_output(&args, "");
        assert_eq!(
            output.pointer(path).unwrap(),
            &expected,
            "{args:?}: {output}"
        );
    }
    let report = json_output(&["report", "CCO"], "");
    assert!(report.is_object());
    assert!(report.to_string().contains("C2H6O"));
}

#[test]
fn batch_commands_keep_valid_records_and_input_order() {
    for (command, input) in [
        ("batch-descriptors", "CCO\ninvalid?\n"),
        ("batch-canonicalize", "CCO\ninvalid?\n"),
        ("batch-fingerprints", "CCO\ninvalid?\n"),
        ("batch-standardize", "CCO\ninvalid?\n"),
        ("batch-similarity", "CCO\tCCO\ninvalid?\tCCO\n"),
        ("batch-substructure", "CCO\t[O]\ninvalid?\t[O]\n"),
        ("batch-reactions", "CCO>>CC=O\ninvalid?\n"),
    ] {
        let output = json_output(&[command], input);
        assert_eq!(output["valid_count"], 1, "{command}: {output}");
        assert_eq!(output["error_count"], 1, "{command}: {output}");
        assert_eq!(output["records"][0]["input_index"], 0);
        assert_eq!(output["records"][1]["input_index"], 1);
        assert!(output["records"][1]["error"].is_string());
    }
    let output = json_output(&["batch-report"], "CCO\n");
    assert!(output.is_object());
    assert!(output.to_string().contains("C2H6O"));
}

#[test]
fn conversion_supports_files_and_errors_use_stderr() {
    let dir = std::env::temp_dir().join(format!("chematic-cli-contracts-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.smi");
    let output = dir.join("output.mol");
    std::fs::write(&input, "CCO").unwrap();
    let result = invoke(
        &[
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "mol",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "",
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(result.stdout.is_empty());
    let mol = std::fs::read_to_string(&output).unwrap();
    let back = invoke(
        &[
            "convert",
            "--input-format",
            "mol",
            "--output-format",
            "smiles",
        ],
        &mol,
    );
    assert!(back.status.success());
    let smiles = String::from_utf8(back.stdout).unwrap();
    let parsed = json_output(&["parse", smiles.trim()], "");
    assert_eq!(parsed["formula"], "C2H6O");
    assert_eq!(parsed["atoms"], 3);
    assert_eq!(parsed["bonds"], 2);
    for args in [
        vec!["parse", "invalid?"],
        vec!["fingerprint", "CCO", "--algorithm", "unknown"],
        vec!["reaction-fingerprint", "CCO>>CC=O", "--mode", "unknown"],
        vec!["batch-descriptors", "--max-records", "0"],
        vec![
            "convert",
            "--input-format",
            "unknown",
            "--output-format",
            "smiles",
        ],
        vec![
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "mol",
            "--input",
            dir.to_str().unwrap(),
        ],
    ] {
        let result = invoke(&args, "CCO\n");
        assert!(!result.status.success(), "{args:?}");
        assert!(result.stdout.is_empty());
        assert!(!result.stderr.is_empty());
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn every_cli_algorithm_and_resource_limit_has_an_executable_contract() {
    for algorithm in [
        "ecfp4",
        "ecfp6",
        "maccs",
        "rdkit_morgan",
        "rdkit_pattern",
        "rdkit_rdk",
        "rdkit_torsion",
    ] {
        let fingerprint = json_output(&["fingerprint", "c1ccccc1O", "--algorithm", algorithm], "");
        assert_eq!(fingerprint["algorithm"], algorithm);
        assert!(fingerprint["popcount"].as_u64().unwrap() > 0);
        let similarity = json_output(
            &[
                "similarity",
                "c1ccccc1O",
                "c1ccccc1O",
                "--algorithm",
                algorithm,
            ],
            "",
        );
        assert_eq!(similarity["similarity"], 1.0);
    }

    for (args, stdin) in [
        (vec!["batch-descriptors", "--max-input-bytes", "2"], "CCO\n"),
        (vec!["batch-descriptors", "--max-line-bytes", "2"], "CCO\n"),
        (vec!["batch-descriptors", "--max-records", "1"], "C\nO\n"),
        (
            vec!["batch-fingerprints", "--algorithm", "unknown"],
            "CCO\n",
        ),
        (
            vec!["batch-similarity", "--algorithm", "unknown"],
            "CCO\tCCO\n",
        ),
        (vec!["batch-similarity"], "missing-tab\n"),
        (vec!["batch-substructure"], "missing-tab\n"),
    ] {
        let output = invoke(&args, stdin);
        if args[0].starts_with("batch-") && !matches!(args[0], "batch-descriptors") {
            // Per-record failures are retained in a successful manifest.
            assert!(output.status.success(), "{args:?}");
            let value: Value = serde_json::from_slice(&output.stdout).unwrap();
            assert_eq!(value["error_count"], 1, "{args:?}: {value}");
            assert_eq!(value["all_succeeded"], false);
        } else {
            assert!(!output.status.success(), "{args:?}");
        }
    }
}

#[test]
fn rdkit_compatible_fingerprints_match_independent_ethanol_oracle() {
    // RDKit 2025.09.3: GetMorganGenerator(radius=2, fpSize=2048) and
    // Chem.PatternFingerprint(mol, fpSize=2048), both on MolFromSmiles("CCO").
    let morgan = json_output(&["fingerprint", "CCO", "--algorithm", "rdkit_morgan"], "");
    assert_eq!(morgan["set_bits"], json!([80, 222, 294, 807, 1057, 1410]));
    let pattern = json_output(&["fingerprint", "CCO", "--algorithm", "rdkit_pattern"], "");
    assert_eq!(
        pattern["set_bits"],
        json!([
            54, 173, 217, 429, 622, 759, 778, 874, 946, 967, 1022, 1033, 1061, 1236, 1289, 1295
        ])
    );
    let batch = json_output(
        &["batch-fingerprints", "--algorithm", "rdkit_morgan"],
        "CCO\n",
    );
    assert_eq!(
        batch["records"][0]["fingerprint"]["set_bits"],
        morgan["set_bits"]
    );
}

#[test]
fn stereo_search_option_matches_rdkit_enantiomer_behavior() {
    // RDKit GetSubstructMatches(query, useChirality=True) distinguishes these
    // alanine enantiomers; its default useChirality=False accepts both.
    let target = "N[C@@H](C)C(=O)O";
    let opposite = "N[C@H](C)C(=O)O";
    let default = json_output(&["substructure", target, opposite], "");
    assert_eq!(default["match_count"], 1);
    assert_eq!(default["use_chirality"], false);
    let strict = json_output(&["substructure", target, opposite, "--use-chirality"], "");
    assert_eq!(strict["match_count"], 0);
    assert_eq!(strict["use_chirality"], true);
    let same = json_output(&["substructure", target, target, "--use-chirality"], "");
    assert_eq!(same["matches"], json!([[0, 1, 2, 3, 4, 5]]));

    let batch = json_output(
        &["batch-substructure", "--use-chirality"],
        &format!("{target}\t{opposite}\n{target}\t{target}\n"),
    );
    assert_eq!(batch["records"][0]["substructure"]["match_count"], 0);
    assert_eq!(batch["records"][1]["substructure"]["match_count"], 1);

    // The same option must also enforce E/Z and reject an @ tag that RDKit
    // removes from a symmetric cyclobutane carbon on MolFromSmiles.
    let ez_opposite = json_output(
        &["substructure", "F/C=C/F", "F/C=C\\F", "--use-chirality"],
        "",
    );
    assert_eq!(ez_opposite["match_count"], 0);
    let symmetric_ring = json_output(
        &[
            "substructure",
            "C1[C@H](F)CC1",
            "C1[C@H](F)CC1",
            "--use-chirality",
        ],
        "",
    );
    assert_eq!(symmetric_ring["match_count"], 0);
    let nonstereogenic_alkene = json_output(
        &["substructure", "C/C=C(/C)C", "C/C=C/C", "--use-chirality"],
        "",
    );
    assert_eq!(nonstereogenic_alkene["match_count"], 0);
}

#[test]
fn checked_conversion_and_mol2_round_trip_are_exposed_by_the_binary() {
    let mol2 = invoke(
        &[
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "mol2",
        ],
        "CCO",
    );
    assert!(mol2.status.success());
    let round_trip = invoke(
        &[
            "convert",
            "--input-format",
            "mol2",
            "--output-format",
            "smiles",
        ],
        &String::from_utf8(mol2.stdout).unwrap(),
    );
    assert!(round_trip.status.success());

    let checked = json_output(
        &[
            "convert-with-report",
            "--input-format",
            "mol2",
            "--output-format",
            "mol2",
            "--policy",
            "preserve_source",
        ],
        "@<TRIPOS>MOLECULE\nwater\n1 0 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 O 0 0 0 O.3 1 HOH 0\n",
    );
    assert_ne!(checked["status"], "refused");
}
