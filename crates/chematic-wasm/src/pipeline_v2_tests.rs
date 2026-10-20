#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_depict::parse_smiles;

    fn safe_config_json(
        force_field: &str,
        stereo_policy: &str,
        ring_torsion_policy: &str,
    ) -> String {
        format!(
            r#"{{
                "embedSeed": 7,
                "maxAttempts": 8,
                "embedTimeoutMs": null,
                "useExpTorsions": false,
                "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false,
                "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false,
                "stereoPolicy": "{stereo_policy}",
                "failOnUnevaluableStereo": false,
                "forceFieldPolicy": "{force_field}",
                "forceFieldMaxIterations": 200,
                "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false,
                "ringTorsionPolicy": "{ring_torsion_policy}",
                "totalTimeoutMs": null
            }}"#
        )
    }

    #[test]
    fn success_path_has_expected_envelope_shape() {
        let mol = parse_smiles("CCCCCCCCCC").expect("decane"); // decane
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], true);
        let result = &value["result"];
        for key in [
            "coords",
            "embedStats",
            "boundAdjustmentReport",
            "torsionKnowledgeReport",
            "ringTorsionEvidence",
            "torsionOptimizationReport",
            "stereoBefore",
            "stereoRepair",
            "stereoAfterRepair",
            "forceField",
            "finalStereo",
            "finalValidation",
            "elapsedMsByStage",
        ] {
            assert!(result.get(key).is_some(), "missing result.{key}");
        }
        assert_eq!(result["coords"].as_array().unwrap().len(), 10);
        assert_eq!(result["boundAdjustmentReport"], serde_json::Value::Null);
        assert_eq!(result["stereoRepair"], serde_json::Value::Null);
        assert_eq!(result["torsionOptimizationReport"], serde_json::Value::Null);
        assert_eq!(result["forceField"]["requestedForceField"], "none");
        assert_eq!(result["forceField"]["actualForceFieldUsed"], "none");
        assert_eq!(
            result["forceField"]["mmff94Termination"],
            serde_json::Value::Null
        );
        assert_eq!(
            result["forceField"]["fallbackReason"],
            serde_json::Value::Null
        );
    }

    #[test]
    fn mmff94_termination_is_exposed_in_json() {
        let mol = parse_smiles("CCCC").expect("butane");
        let config = safe_config_json("mmff94_bond_angle_strict", "ignore", "diagnostic_only");
        let value: serde_json::Value =
            serde_json::from_str(&embed_pipeline_v2_json(&mol, &config)).expect("valid JSON");
        assert_eq!(value["ok"], true, "{value}");
        let ff = &value["result"]["forceField"];
        let reason = ff["mmff94Termination"].as_str().expect("MMFF94 reason");
        assert!(matches!(
            reason,
            "gradient_converged" | "iteration_limit" | "constraint_rejected_fallback"
        ));
        assert_eq!(ff["converged"] == true, reason == "gradient_converged");
    }

    #[test]
    fn mmff94_failure_json_preserves_optimizer_stop_reason() {
        use chematic_3d::minimize::{
            ForceFieldBridgeError, MinimizationFailureDetail, MinimizationFailureReason,
        };
        use chematic_ff::Mmff94TerminationReason;

        let error =
            ForceFieldBridgeError::MinimizationFailed(Box::new(MinimizationFailureDetail {
                policy: ForceFieldPolicy::Mmff94BondAngleStrict,
                reason: MinimizationFailureReason::ExcessiveResidualForce,
                converged: false,
                iterations: 20,
                mmff94_termination: Some(Mmff94TerminationReason::ConstraintRejectedFallback),
                max_residual_force: 557.0,
                worst_bond_length: 2.1,
                distance_geometry_v2_retry_attempted: false,
            }));
        let value = serde_json::to_value(force_field_bridge_error_json(&error))
            .expect("serializable typed refusal");
        assert_eq!(value["kind"], "minimization_failed");
        assert_eq!(value["reason"], "excessive_residual_force");
        assert_eq!(value["mmff94Termination"], "constraint_rejected_fallback");
    }

    #[test]
    fn same_seed_is_reproducible() {
        // Wall-clock timing (`elapsedMsByStage`) is explicitly excluded from this
        // determinism requirement (it varies with system load, not with the
        // computation itself) -- strip it before comparing, matching the spec's
        // "wall-clock timing numeric equality is not required" instruction.
        fn without_timings(json: &str) -> serde_json::Value {
            let mut value: serde_json::Value = serde_json::from_str(json).unwrap();
            value["result"]["elapsedMsByStage"].take();
            value
        }
        let mol = parse_smiles("CCCCCCCCCC").expect("decane");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json1 = embed_pipeline_v2_json(&mol, &config);
        let json2 = embed_pipeline_v2_json(&mol, &config);
        assert_eq!(without_timings(&json1), without_timings(&json2));
    }

    #[test]
    fn failure_path_has_expected_envelope_shape() {
        // A saturated small ring (cyclohexane) fused to an acyclic tail: requesting
        // small-ring torsions under the fail-closed policy is a reliably-reachable
        // typed failure (see the Python binding's own test suite for the same
        // fixture/rationale).
        let mol = parse_smiles("C1CCCCC1CCCCCCCCCCCC").expect("cyclohexane+chain");
        let config = r#"{
            "embedSeed": 7,
            "maxAttempts": 8,
            "embedTimeoutMs": null,
            "useExpTorsions": false,
            "useSmallRingTorsions": true,
            "useMacrocycleTorsions": false,
            "useMacrocycle14Bounds": false,
            "includeLegacyTorsionHeuristic": false,
            "stereoPolicy": "ignore",
            "failOnUnevaluableStereo": false,
            "forceFieldPolicy": "dreiding",
            "forceFieldMaxIterations": 200,
            "gateMmff94TorsionOop": false,
            "gateMmff94StretchBend": false,
            "ringTorsionPolicy": "fail_closed",
            "totalTimeoutMs": null
        }"#;
        let json = embed_pipeline_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], false);
        let error = &value["error"];
        assert_eq!(
            error["cause"]["kind"],
            "ring_torsion_application_unsupported"
        );
        assert_eq!(error["stage"], "torsion_optimization");
        assert_eq!(error["coordsAreDiagnosticOnly"], true);
        assert!(error["lastKnownCoords"].as_array().is_some());
        assert!(!value["error"].as_object().unwrap().contains_key("coords"));
        assert!(error["diagnostics"]["embedStats"].is_object());
    }

    #[test]
    fn unknown_config_field_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["notARealField"] = serde_json::json!(true);
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["stage"], "wasm_input_validation");
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn missing_required_field_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config.as_object_mut().unwrap().remove("maxAttempts");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn pre_priority2_config_json_without_gate_stretch_bend_still_parses() {
        // Backward-compat regression (Priority 2, issue #227): a caller's
        // pre-existing 15-field config JSON (no `gateMmff94StretchBend` at
        // all) must keep working exactly as before, defaulting to `false`
        // -- not become a "missing required field" error just because a new
        // gate dimension was added. Contrast with
        // `missing_required_field_is_rejected` above: that one IS still a
        // real required field.
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config
            .as_object_mut()
            .unwrap()
            .remove("gateMmff94StretchBend");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(
            value["ok"], true,
            "old caller config without gateMmff94StretchBend must still succeed: {value:?}"
        );
    }

    #[test]
    fn missing_nullable_timeout_field_is_still_rejected() {
        // embedTimeoutMs/totalTimeoutMs may be `null`, but the key itself must be
        // present -- a caller omitting it entirely must not be silently treated
        // as null (the double-Option serde trick this test pins).
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config.as_object_mut().unwrap().remove("embedTimeoutMs");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn null_timeout_value_is_accepted() {
        let mol = parse_smiles("CC").expect("ethane");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
    }

    #[test]
    fn present_timeout_value_is_accepted() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["totalTimeoutMs"] = serde_json::json!(60_000);
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
    }

    #[test]
    fn unknown_enum_value_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["stereoPolicy"] = serde_json::json!("not_a_real_policy");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn out_of_range_integer_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["maxAttempts"] = serde_json::json!(-1);
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn malformed_json_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let json = embed_pipeline_v2_json(&mol, "{not valid json");
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    /// Builds a `MolHandle` directly from `chematic_core`'s own builder, bypassing
    /// `parse_smiles`'s `enforce_wasm_molecule_size` check entirely -- that check
    /// itself constructs a `JsValue` on the error path, which aborts the process
    /// when run as a native (non-wasm32) test outside a real JS host. This lets
    /// `embed_pipeline_v2_json`'s *own* atom-limit check (which returns a plain
    /// `String`, never a `JsValue`) be exercised in isolation.
    fn oversized_mol_handle(atom_count: usize) -> MolHandle {
        let mut builder = chematic_core::MoleculeBuilder::new();
        for _ in 0..atom_count {
            builder.add_atom(chematic_core::Atom::new(chematic_core::Element::C));
        }
        MolHandle {
            inner: std::rc::Rc::new(builder.build()),
        }
    }

    #[test]
    fn atom_limit_exceeded_is_fail_closed() {
        let mol = oversized_mol_handle(WASM_MAX_ATOMS + 1);
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "atom_limit_exceeded");
        assert_eq!(value["error"]["cause"]["limit"], WASM_MAX_ATOMS);
        assert_eq!(value["error"]["stage"], "wasm_input_validation");
    }

    #[test]
    fn oversized_config_json_is_fail_closed() {
        let mol = parse_smiles("CC").expect("ethane");
        let padding = " ".repeat(WASM_MAX_INPUT_BYTES + 1);
        let json = embed_pipeline_v2_json(&mol, &padding);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "input_too_large");
    }

    #[test]
    fn finite_f64_serializes_as_number_non_finite_as_null() {
        assert_eq!(serde_json::to_string(&FiniteF64(1.5)).unwrap(), "1.5");
        assert_eq!(serde_json::to_string(&FiniteF64(f64::NAN)).unwrap(), "null");
        assert_eq!(
            serde_json::to_string(&FiniteF64(f64::INFINITY)).unwrap(),
            "null"
        );
        assert_eq!(
            serde_json::to_string(&FiniteF64(f64::NEG_INFINITY)).unwrap(),
            "null"
        );
        // Full round-trippable precision -- not rounded to a fixed number of
        // decimal places, unlike this crate's older 3D JSON helpers.
        let precise = 1.0 / 3.0;
        let rendered = serde_json::to_string(&FiniteF64(precise)).unwrap();
        let parsed: f64 = rendered.parse().unwrap();
        assert_eq!(parsed, precise);
    }

    fn coord_distance(coords: &serde_json::Value, i: usize, j: usize) -> f64 {
        let p = |idx: usize| {
            let c = &coords[idx];
            (
                c[0].as_f64().unwrap(),
                c[1].as_f64().unwrap(),
                c[2].as_f64().unwrap(),
            )
        };
        let (x1, y1, z1) = p(i);
        let (x2, y2, z2) = p(j);
        ((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2)).sqrt()
    }

    /// Same property `test_conformer_ensemble.py`'s `_assert_consistent_indexing`
    /// checks on the Python side (issue #172): every real bond must have a
    /// plausible length in the returned coords, and no two atoms may coincide --
    /// exactly what breaks if a canonicalize-then-reparse bug (like
    /// `conformer_ensemble_json`'s still-open one) silently desyncs the atom index
    /// space between `mol` and the returned coordinates.
    fn assert_consistent_indexing(mol: &MolHandle, coords: &serde_json::Value) {
        let coords_arr = coords.as_array().expect("coords must be an array");
        assert_eq!(coords_arr.len(), mol.inner.atom_count());
        for (_bond_idx, bond) in mol.inner.bonds() {
            let d = coord_distance(coords, bond.atom1.0 as usize, bond.atom2.0 as usize);
            assert!(
                (0.6..=2.6).contains(&d),
                "bond {}-{} has length {d:.3} -- outside sane range [0.6, 2.6]; \
                 likely atom-index mismatch",
                bond.atom1.0,
                bond.atom2.0
            );
        }
        let n = coords_arr.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = coord_distance(coords, i, j);
                assert!(
                    d >= 0.4,
                    "atoms {i},{j} are {d:.3} apart -- degenerate/clashing geometry"
                );
            }
        }
    }

    #[test]
    fn atom_order_is_preserved_naphthalene() {
        // Naphthalene reorders under canonicalization (same fixture the Python
        // binding's atom-order test suite uses) -- coords must stay indexed to
        // THIS Mol's own atom order, not a reparsed/canonicalized copy.
        let mol = parse_smiles("c1ccc2ccccc2c1").expect("naphthalene");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
        assert_consistent_indexing(&mol, &value["result"]["coords"]);
    }

    #[test]
    fn atom_order_is_preserved_branched_and_aspirin() {
        for smiles in ["CCC(C)C", "CC(=O)Oc1ccccc1C(=O)O"] {
            let mol = parse_smiles(smiles).unwrap_or_else(|_| panic!("{smiles}"));
            let config = safe_config_json("none", "ignore", "fail_closed");
            let json = embed_pipeline_v2_json(&mol, &config);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            assert_eq!(value["ok"], true, "{smiles}");
            assert_consistent_indexing(&mol, &value["result"]["coords"]);
        }
    }

    #[test]
    fn negative_control_reproduces_old_reparse_bug_shape() {
        // Hand-reconstructs the exact anti-pattern `conformer_ensemble_json` still
        // has (canonicalize -> reparse -> compute -> return coords indexed to the
        // REPARSED molecule) and confirms cross-indexing it against the ORIGINAL
        // mol is caught -- proving `assert_consistent_indexing` can actually detect
        // a real atom-order mismatch, not just that it happens to pass.
        let mol = parse_smiles("c1ccc2ccccc2c1").expect("naphthalene");
        let canonical = mol.canonical_smiles();
        assert_ne!(
            canonical, "c1ccc2ccccc2c1",
            "expected canonicalization to reorder naphthalene's atoms"
        );
        let reparsed = parse_smiles(&canonical).expect("reparse canonical form");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&reparsed, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);

        let result = std::panic::catch_unwind(|| {
            assert_consistent_indexing(&mol, &value["result"]["coords"]);
        });
        assert!(
            result.is_err(),
            "cross-indexing the reparsed molecule's coords against the original \
             mol should have failed a bond-length sanity check"
        );
    }

    // -----------------------------------------------------------------------
    // enforce_chirality (v0.14.0, issue #285's E/Z bound fix) -- WASM parity
    // -----------------------------------------------------------------------

    fn enforce_chirality_config_json(stereo_policy: &str, enforce_chirality: bool) -> String {
        format!(
            r#"{{
                "embedSeed": 0,
                "maxAttempts": 1,
                "embedTimeoutMs": null,
                "useExpTorsions": false,
                "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false,
                "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false,
                "stereoPolicy": "{stereo_policy}",
                "failOnUnevaluableStereo": false,
                "forceFieldPolicy": "none",
                "forceFieldMaxIterations": 200,
                "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false,
                "ringTorsionPolicy": "fail_closed",
                "totalTimeoutMs": null,
                "enforceChirality": {enforce_chirality}
            }}"#
        )
    }

    #[test]
    fn enforce_chirality_true_fixes_but2ene_z_raw_embed() {
        // Direct WASM-level confirmation that `enforceChirality` in the JSON
        // config reaches distance_geometry_v2.rs's apply_declared_ez_bounds
        // (issue #285): but2ene_Z is the exact molecule that fix targets --
        // raw embedding (no force field) must satisfy declared E/Z once
        // enforceChirality is set, across multiple seeds, matching the
        // Rust-level corpus measurement and the Python binding's parity test.
        let mol = parse_smiles(r"C/C=C\C").expect("but2ene_Z");
        for seed in 0..5u64 {
            let config = format!(
                r#"{{
                    "embedSeed": {seed},
                    "maxAttempts": 1,
                    "embedTimeoutMs": null,
                    "useExpTorsions": false,
                    "useSmallRingTorsions": false,
                    "useMacrocycleTorsions": false,
                    "useMacrocycle14Bounds": false,
                    "includeLegacyTorsionHeuristic": false,
                    "stereoPolicy": "ignore",
                    "failOnUnevaluableStereo": false,
                    "forceFieldPolicy": "none",
                    "forceFieldMaxIterations": 200,
                    "gateMmff94TorsionOop": false,
                    "gateMmff94StretchBend": false,
                    "ringTorsionPolicy": "fail_closed",
                    "totalTimeoutMs": null,
                    "enforceChirality": true
                }}"#
            );
            let json = embed_pipeline_v2_json(&mol, &config);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            assert_eq!(value["ok"], true, "seed {seed}: {json}");
            assert_eq!(
                value["result"]["finalStereo"]["isFullySatisfied"], true,
                "seed {seed}: raw embed must already satisfy declared E/Z"
            );
        }
    }

    #[test]
    fn enforce_chirality_defaults_false_missing_field_still_parses() {
        // #[serde(default)] precedent (matches gateMmff94StretchBend): configs
        // written before this field existed must keep parsing, as `false`.
        let mol = parse_smiles("CCCCCCCCCC").expect("decane");
        let config = safe_config_json("none", "ignore", "fail_closed"); // no enforceChirality key
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "{json}");
    }

    #[test]
    fn enforce_chirality_with_repair_and_verify_is_allowed() {
        // Revised 2026-08-24 (issue #291 Step A): this combination was
        // previously rejected as invalid_configuration -- now validated (see
        // `chematic_3d::pipeline_v2`'s revised Stage 1 doc entry and
        // `crates/chematic-3d/examples/issue291_repair_policy_measurement.rs`).
        let mol = parse_smiles(r"C/C=C\C").expect("but2ene_Z");
        let config = enforce_chirality_config_json("repair_and_verify", true);
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "{json}");
        assert_eq!(
            value["result"]["finalStereo"]["isFullySatisfied"], true,
            "{json}"
        );
    }

    // -----------------------------------------------------------------------
    // expand_implicit_h_through_pipeline / pipeline_v2_stereo_safe_config_json
    // (issue #291/#383) -- WASM parity
    // -----------------------------------------------------------------------

    const TESTOSTERONE: &str = "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O";

    #[test]
    fn expand_implicit_h_through_pipeline_requires_enforce_chirality() {
        let mol = parse_smiles(TESTOSTERONE).expect("testosterone");
        let config = r#"{
                "embedSeed": 0,
                "maxAttempts": 8,
                "embedTimeoutMs": null,
                "useExpTorsions": false,
                "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false,
                "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false,
                "stereoPolicy": "repair_and_verify",
                "failOnUnevaluableStereo": false,
                "forceFieldPolicy": "mmff94_with_uff_fallback",
                "forceFieldMaxIterations": 200,
                "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false,
                "ringTorsionPolicy": "diagnostic_only",
                "totalTimeoutMs": null,
                "enforceChirality": false,
                "expandImplicitHThroughPipeline": true
            }"#;
        let json = embed_pipeline_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false, "{json}");
        assert_eq!(
            value["error"]["cause"]["kind"], "invalid_configuration",
            "{json}"
        );
        assert_eq!(value["error"]["stage"], "validate_config", "{json}");
    }

    #[test]
    fn stereo_safe_config_json_has_expected_shape() {
        let json = pipeline_v2_stereo_safe_config_json("mmff94_with_uff_fallback", "fail_closed");
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], true);
        let config = &value["config"];
        assert_eq!(config["stereoPolicy"], "repair_and_verify");
        assert_eq!(config["enforceChirality"], true);
        assert_eq!(config["expandImplicitHThroughPipeline"], true);
        assert_eq!(config["forceFieldPolicy"], "mmff94_with_uff_fallback");
        assert_eq!(config["ringTorsionPolicy"], "fail_closed");
    }

    #[test]
    fn stereo_safe_config_json_rejects_unknown_force_field() {
        let json = pipeline_v2_stereo_safe_config_json("not_a_real_force_field", "fail_closed");
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false, "{json}");
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config", "{json}");
    }

    #[cfg(target_arch = "wasm32")]
    #[test]
    fn stereo_safe_config_json_fixes_testosterone_via_wasm_binding() {
        // Same seed/configuration already Rust-level tested and cross-checked
        // against an independent oracle in pipeline_v2.rs's own
        // stereo_safe_matches_the_hand_built_configuration_above test.
        let mol = parse_smiles(TESTOSTERONE).expect("testosterone");
        let generated =
            pipeline_v2_stereo_safe_config_json("mmff94_with_uff_fallback", "diagnostic_only");
        let mut config: serde_json::Value = serde_json::from_str(&generated).unwrap();
        config["config"]["embedSeed"] = serde_json::json!(0);
        let config_json = config["config"].to_string();

        let json = embed_pipeline_v2_json(&mol, &config_json);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "{json}");
        assert_eq!(
            value["result"]["finalStereo"]["isFullySatisfied"], true,
            "{json}"
        );
        assert_eq!(value["result"]["coords"].as_array().unwrap().len(), 20);
    }

    #[test]
    fn expand_implicit_h_through_pipeline_is_noop_without_declared_stereo() {
        let mol = parse_smiles("CCCCCCCCCC").expect("decane"); // no declared stereo
        let base = embed_pipeline_v2_json(
            &mol,
            &enforce_chirality_config_json("repair_and_verify", true),
        );
        let expanded_config = {
            let mut c: serde_json::Value =
                serde_json::from_str(&enforce_chirality_config_json("repair_and_verify", true))
                    .unwrap();
            c["expandImplicitHThroughPipeline"] = serde_json::json!(true);
            c.to_string()
        };
        let expanded = embed_pipeline_v2_json(&mol, &expanded_config);
        let base_value: serde_json::Value = serde_json::from_str(&base).unwrap();
        let expanded_value: serde_json::Value = serde_json::from_str(&expanded).unwrap();
        assert_eq!(
            base_value["result"]["coords"],
            expanded_value["result"]["coords"]
        );
    }

    // -----------------------------------------------------------------------
    // Cross-binding parity: validation/pipeline_v2_wasm_parity_fixtures.json
    //
    // Generated once by scripts/generators/gen_pipeline_v2_wasm_parity_fixtures.py via the
    // Python binding (crates/chematic-py/src/pipeline_v2.rs). Both the WASM
    // binding under test here and Python ultimately call the exact same
    // `chematic_3d::pipeline_v2::embed_pipeline_v2` -- so this is less "does the
    // algorithm agree" (trivially true) and more "does each binding's own JSON
    // conversion faithfully represent the same underlying evidence, for a fixed
    // molecule/config/seed." Node's pipeline_v2_parity.test.mjs reads the same
    // file and checks the WASM JS entry point the same way.
    // -----------------------------------------------------------------------

    fn parity_fixtures() -> serde_json::Value {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../validation/pipeline_v2_wasm_parity_fixtures.json");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
        serde_json::from_str(&text).expect("valid JSON")
    }

    fn config_json_from_fixture(fixture: &serde_json::Value) -> String {
        fixture["config"].to_string()
    }

    #[test]
    fn wasm_binding_matches_python_reference_fixtures() {
        let fixtures = parity_fixtures();
        for fixture in fixtures["fixtures"].as_array().unwrap() {
            let name = fixture["name"].as_str().unwrap();
            let smiles = fixture["smiles"].as_str().unwrap();
            let config_json = config_json_from_fixture(fixture);
            let mol = parse_smiles(smiles).unwrap_or_else(|_| panic!("{name}: parse {smiles}"));

            assert_eq!(
                mol.inner.atom_count() as i64,
                fixture["atomCount"].as_i64().unwrap(),
                "{name}: atom count"
            );

            let json = embed_pipeline_v2_json(&mol, &config_json);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            let expected_ok = fixture["ok"].as_bool().unwrap();
            assert_eq!(value["ok"], expected_ok, "{name}: ok");

            if expected_ok {
                let result = &value["result"];
                assert_eq!(
                    result["coords"].as_array().unwrap().len() as i64,
                    fixture["coordsLength"].as_i64().unwrap(),
                    "{name}: coords length"
                );
                assert_consistent_indexing(&mol, &result["coords"]);
                assert_eq!(
                    result["stereoBefore"]["nDeclared"], fixture["stereoBeforeDeclared"],
                    "{name}: stereoBefore.nDeclared"
                );
                assert_eq!(
                    result["finalStereo"]["nDeclared"], fixture["stereoAfterDeclared"],
                    "{name}: finalStereo.nDeclared"
                );
                assert_eq!(
                    result["finalStereo"]["nViolations"], fixture["stereoAfterViolations"],
                    "{name}: finalStereo.nViolations"
                );
                assert_eq!(
                    result["forceField"]["requestedForceField"], fixture["forceFieldRequested"],
                    "{name}: forceField.requestedForceField"
                );
                assert_eq!(
                    result["forceField"]["actualForceFieldUsed"], fixture["forceFieldActual"],
                    "{name}: forceField.actualForceFieldUsed"
                );
                assert_eq!(
                    !result["forceField"]["fallbackReason"].is_null(),
                    fixture["hasFallback"].as_bool().unwrap(),
                    "{name}: forceField fallback presence"
                );
                assert_eq!(
                    result["finalValidation"]["sound"], fixture["sound"],
                    "{name}: finalValidation.sound"
                );
            } else {
                assert_eq!(
                    value["error"]["stage"], fixture["stage"],
                    "{name}: error.stage"
                );
                assert_eq!(
                    value["error"]["cause"]["kind"], fixture["causeKind"],
                    "{name}: error.cause.kind"
                );
            }
        }
    }

    /// Same fixtures, but computed via the raw `chematic_3d::pipeline_v2` API
    /// directly (no JSON at all) -- confirms the WASM JSON conversion (tested
    /// above) is not itself hiding a divergence from what the underlying Rust
    /// call actually returns.
    #[test]
    fn raw_rust_api_matches_python_reference_fixtures() {
        let fixtures = parity_fixtures();
        for fixture in fixtures["fixtures"].as_array().unwrap() {
            let name = fixture["name"].as_str().unwrap();
            let smiles = fixture["smiles"].as_str().unwrap();
            let cfg = &fixture["config"];
            let mol = chematic_smiles::parse(smiles).unwrap_or_else(|_| panic!("{name}: parse"));

            let stereo_policy = match cfg["stereoPolicy"].as_str().unwrap() {
                "ignore" => pv2::StereoPolicy::Ignore,
                "verify_only" => pv2::StereoPolicy::VerifyOnly,
                "repair_and_verify" => pv2::StereoPolicy::RepairAndVerify,
                other => panic!("unknown stereoPolicy {other}"),
            };
            let ring_torsion_policy = match cfg["ringTorsionPolicy"].as_str().unwrap() {
                "fail_closed" => pv2::RingTorsionApplicationPolicy::FailClosed,
                "diagnostic_only" => pv2::RingTorsionApplicationPolicy::DiagnosticOnly,
                other => panic!("unknown ringTorsionPolicy {other}"),
            };
            let force_field_policy = match cfg["forceFieldPolicy"].as_str().unwrap() {
                "mmff94_bond_angle_strict" => ForceFieldPolicy::Mmff94BondAngleStrict,
                "mmff94_with_uff_fallback" => ForceFieldPolicy::Mmff94WithUffFallback,
                "uff_only" => ForceFieldPolicy::UffOnly,
                "dreiding" => ForceFieldPolicy::Dreiding,
                "none" => ForceFieldPolicy::None,
                other => panic!("unknown forceFieldPolicy {other}"),
            };
            let config = pv2::PipelineV2Config {
                embed: EmbedParameters {
                    random_seed: cfg["embedSeed"].as_u64().unwrap(),
                    max_attempts: cfg["maxAttempts"].as_u64().unwrap() as usize,
                    timeout_ms: cfg["embedTimeoutMs"].as_u64(),
                    use_exp_torsions: cfg["useExpTorsions"].as_bool().unwrap(),
                    use_small_ring_torsions: cfg["useSmallRingTorsions"].as_bool().unwrap(),
                    use_macrocycle_torsions: cfg["useMacrocycleTorsions"].as_bool().unwrap(),
                    use_macrocycle_14_bounds: cfg["useMacrocycle14Bounds"].as_bool().unwrap(),
                    enforce_chirality: cfg["enforceChirality"].as_bool().unwrap_or(false),
                    ..EmbedParameters::default()
                },
                torsion_optimization: TorsionOptimizationConfig::default(),
                include_legacy_torsion_heuristic: cfg["includeLegacyTorsionHeuristic"]
                    .as_bool()
                    .unwrap(),
                stereo_policy,
                fail_on_unevaluable_stereo: cfg["failOnUnevaluableStereo"].as_bool().unwrap(),
                force_field_policy,
                force_field_max_iterations: cfg["forceFieldMaxIterations"].as_u64().unwrap()
                    as usize,
                gate_mmff94_torsion_oop: cfg["gateMmff94TorsionOop"].as_bool().unwrap(),
                gate_mmff94_stretch_bend: cfg["gateMmff94StretchBend"].as_bool().unwrap_or(false),
                ring_torsion_policy,
                total_timeout_ms: cfg["totalTimeoutMs"].as_u64(),
                expand_implicit_h_through_pipeline: false,
            };

            let expected_ok = fixture["ok"].as_bool().unwrap();
            match pv2::embed_pipeline_v2(&mol, &config) {
                Ok(result) => {
                    assert!(expected_ok, "{name}: expected failure, got success");
                    assert_eq!(
                        result.coords.atom_count() as i64,
                        fixture["coordsLength"].as_i64().unwrap(),
                        "{name}: coords length"
                    );
                    assert_eq!(
                        result.stereo_before.n_declared() as i64,
                        fixture["stereoBeforeDeclared"].as_i64().unwrap(),
                        "{name}: stereo_before.n_declared"
                    );
                    assert_eq!(
                        force_field_policy_str(result.force_field.actual_force_field_used),
                        fixture["forceFieldActual"].as_str().unwrap(),
                        "{name}: force_field.actual_force_field_used"
                    );
                    assert_eq!(
                        result.force_field.fallback_reason.is_some(),
                        fixture["hasFallback"].as_bool().unwrap(),
                        "{name}: force_field fallback presence"
                    );
                    assert_eq!(
                        result.final_validation.sound,
                        fixture["sound"].as_bool().unwrap(),
                        "{name}: final_validation.sound"
                    );
                }
                Err(failure) => {
                    assert!(!expected_ok, "{name}: expected success, got failure");
                    assert_eq!(
                        snake_case_debug(&failure.stage),
                        fixture["stage"].as_str().unwrap(),
                        "{name}: stage"
                    );
                    assert_eq!(
                        failure_cause_kind_str(&failure.cause),
                        fixture["causeKind"].as_str().unwrap(),
                        "{name}: cause kind"
                    );
                }
            }
        }
    }

    /// `FailureCauseJson`'s serde tag value for a cause, without needing to
    /// serialize the whole struct first -- used only by the raw-API parity test.
    fn failure_cause_kind_str(cause: &pv2::PipelineV2FailureCause) -> String {
        let json = serde_json::to_value(failure_cause_json(cause)).unwrap();
        json["kind"].as_str().unwrap().to_string()
    }
}
