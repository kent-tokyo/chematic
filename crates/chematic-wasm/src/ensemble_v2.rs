//! WASM binding for `chematic_3d::embed_ensemble_v2` (A2.1).
//!
//! Mirrors the Python binding (`crates/chematic-py/src/ensemble_v2.rs`, already
//! on `main`) as closely as JS/JSON allow, and reuses the WASM `pipeline_v2.rs`
//! binding's own JSON conventions (envelope shape, camelCase keys, snake_case
//! policy/force-field strings, `FiniteF64`, per-attempt failure shape) rather
//! than inventing a second scheme.
//!
//! # Config JSON shape
//!
//! `{"perConformer": <the same object embed_pipeline_v2_json's config_json
//! parses>, "count": <u32>, "baseSeed": <u64>, "rmsdThreshold": <f64>,
//! "useSymmetricRmsdPruning": <bool>, "ensembleTimeoutMs": <u64 | null>}`.
//! Every field is required (deny_unknown_fields, no silent defaults) --
//! matches `PipelineV2ConfigJson`'s own closed-config convention in this
//! crate rather than the Python binding's optional-keyword-argument defaults
//! (`rmsd_threshold = 0.5`, etc.): this is a brand-new WASM API with no
//! existing callers to preserve compatibility for, so it follows this
//! crate's dominant "explicit, nothing silently defaulted" JSON convention
//! instead. `ensembleTimeoutMs` uses the same "key must be present, value may
//! be null" convention as `embedTimeoutMs`/`totalTimeoutMs`.
//!
//! # Error asymmetry with `embed_pipeline_v2_json`
//!
//! Same asymmetry as the Python binding: `embed_ensemble_v2` itself only
//! rejects a config that could never succeed regardless of the molecule
//! (currently: an invalid `rmsdThreshold`) -- an ensemble where every attempt
//! failed and zero conformers were kept is still a normal `{"ok": true,
//! "result": {...}}` response, with per-attempt detail in `result.attempts`.
//! A caller must not assume `ok: true` implies `result.conformers.length >
//! 0`.
//!
//! Every JSON output is a tagged union with `schemaVersion: 1`, same as
//! `embed_pipeline_v2_json`. Never throws.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use chematic_3d::{
    ConformerAttempt, ConformerDisposition, EnsembleV2Config, EnsembleV2Result, embed_ensemble_v2,
};

use crate::pipeline_v2::{
    ErrorEnvelopeJson, FailureCauseJson, FiniteF64, ForceFieldBridgeErrorJson,
    PipelineV2ConfigJson, coords_to_json, deserialize_present, error_envelope_json,
    force_field_bridge_error_json, force_field_policy_str, snake_case_debug, wasm_input_error_json,
};
use crate::{MolHandle, WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES};

const SCHEMA_VERSION: u32 = 1;

// ---------------------------------------------------------------------------
// Config: input JSON shape
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
struct EnsembleV2ConfigJson {
    per_conformer: PipelineV2ConfigJson,
    count: usize,
    base_seed: u64,
    rmsd_threshold: f64,
    use_symmetric_rmsd_pruning: bool,
    #[serde(deserialize_with = "deserialize_present")]
    ensemble_timeout_ms: Option<Option<u64>>,
}

impl EnsembleV2ConfigJson {
    fn into_ensemble_config(self) -> Result<EnsembleV2Config, String> {
        let ensemble_timeout_ms = self.ensemble_timeout_ms.ok_or_else(|| {
            "missing field `ensembleTimeoutMs` (must be present; value may be null)".to_string()
        })?;
        Ok(EnsembleV2Config {
            per_conformer: self.per_conformer.into_pipeline_config()?,
            count: self.count,
            base_seed: self.base_seed,
            rmsd_threshold: self.rmsd_threshold,
            use_symmetric_rmsd_pruning: self.use_symmetric_rmsd_pruning,
            ensemble_timeout_ms,
        })
    }
}

fn parse_ensemble_config(config_json: &str) -> Result<EnsembleV2Config, String> {
    let parsed: EnsembleV2ConfigJson =
        serde_json::from_str(config_json).map_err(|e| e.to_string())?;
    parsed.into_ensemble_config()
}

// ---------------------------------------------------------------------------
// Result: JSON shape
// ---------------------------------------------------------------------------

#[derive(Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum ConformerDispositionJson {
    Kept {
        #[serde(rename = "conformerIndex")]
        conformer_index: usize,
    },
    PrunedAsDuplicate {
        #[serde(rename = "representativeAttemptIndex")]
        representative_attempt_index: usize,
        rmsd: FiniteF64,
        symmetric: bool,
    },
}

fn conformer_disposition_json(d: &ConformerDisposition) -> ConformerDispositionJson {
    match d {
        ConformerDisposition::Kept { conformer_index } => ConformerDispositionJson::Kept {
            conformer_index: *conformer_index,
        },
        ConformerDisposition::PrunedAsDuplicate {
            representative_attempt_index,
            rmsd,
            symmetric,
        } => ConformerDispositionJson::PrunedAsDuplicate {
            representative_attempt_index: *representative_attempt_index,
            rmsd: (*rmsd).into(),
            symmetric: *symmetric,
        },
        // `ConformerDisposition` is `#[non_exhaustive]`: mirrors the Python
        // binding's own refusal to paper over an unhandled future variant --
        // there is no sensible JSON to emit for it, so this stays a hard
        // panic. Not caught anywhere on this path -- same convention
        // `embed_pipeline_v2_json` itself already uses for its own internal
        // invariant violations (e.g. `result_to_json`'s
        // `expect("index is within conformer_count()")` below), not a new
        // failure mode introduced by this binding.
        other => panic!("unhandled ConformerDisposition variant in WASM binding: {other:?}"),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ConformerSuccessJson {
    energy: Option<FiniteF64>,
    actual_force_field_used: String,
    fallback_reason: Option<ForceFieldBridgeErrorJson>,
    disposition: ConformerDispositionJson,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ConformerAttemptJson {
    attempt_index: usize,
    seed: u64,
    outcome: &'static str,
    success: Option<ConformerSuccessJson>,
    failure: Option<ErrorEnvelopeJson>,
}

fn conformer_attempt_json(a: &ConformerAttempt) -> ConformerAttemptJson {
    match &a.outcome {
        Ok(success) => ConformerAttemptJson {
            attempt_index: a.attempt_index,
            seed: a.seed,
            outcome: "success",
            success: Some(ConformerSuccessJson {
                energy: success.energy.map(Into::into),
                actual_force_field_used: force_field_policy_str(success.actual_force_field_used)
                    .to_string(),
                fallback_reason: success
                    .fallback_reason
                    .as_ref()
                    .map(force_field_bridge_error_json),
                disposition: conformer_disposition_json(&success.disposition),
            }),
            failure: None,
        },
        Err(failure) => ConformerAttemptJson {
            attempt_index: a.attempt_index,
            seed: a.seed,
            outcome: "failure",
            success: None,
            failure: Some(error_envelope_json(failure)),
        },
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ConformerProvenanceJson {
    attempt_index: usize,
    seed: u64,
    energy: Option<FiniteF64>,
    actual_force_field_used: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct EnsembleV2ResultJson {
    conformers: Vec<Vec<[FiniteF64; 3]>>,
    conformer_provenance: Vec<ConformerProvenanceJson>,
    attempts: Vec<ConformerAttemptJson>,
    mixed_force_field: bool,
    termination: String,
    requested_count: usize,
}

#[derive(Serialize)]
struct SuccessEnvelopeJson {
    #[serde(rename = "schemaVersion")]
    schema_version: u32,
    ok: bool,
    result: EnsembleV2ResultJson,
}

/// Mirrors `run_embed_ensemble_v2::ensemble_v2_result_dict`'s
/// `conformer_provenance` construction: reverse-maps each `Kept{conformer_index}`
/// disposition back to the attempt that produced it, so a caller never has to
/// scan `attempts` to find out what a given `result.conformers[i]` actually is.
fn result_to_json(r: &EnsembleV2Result) -> String {
    let conformer_count = r.ensemble.conformer_count();
    let conformers: Vec<Vec<[FiniteF64; 3]>> = (0..conformer_count)
        .map(|i| {
            coords_to_json(
                r.ensemble
                    .get_conformer(i)
                    .expect("index is within conformer_count()"),
            )
        })
        .collect();

    let mut provenance_by_index: Vec<Option<ConformerProvenanceJson>> =
        (0..conformer_count).map(|_| None).collect();
    for attempt in &r.attempts {
        if let Ok(success) = &attempt.outcome
            && let ConformerDisposition::Kept { conformer_index } = &success.disposition
        {
            provenance_by_index[*conformer_index] = Some(ConformerProvenanceJson {
                attempt_index: attempt.attempt_index,
                seed: attempt.seed,
                energy: success.energy.map(Into::into),
                actual_force_field_used: force_field_policy_str(success.actual_force_field_used)
                    .to_string(),
            });
        }
    }
    let conformer_provenance: Vec<ConformerProvenanceJson> = provenance_by_index
        .into_iter()
        .enumerate()
        .map(|(i, entry)| {
            entry.unwrap_or_else(|| {
                panic!(
                    "conformer {i} of {conformer_count} has no Kept disposition in `attempts` \
                     -- embed_ensemble_v2's own invariant is broken, not a case to paper over"
                )
            })
        })
        .collect();

    let envelope = SuccessEnvelopeJson {
        schema_version: SCHEMA_VERSION,
        ok: true,
        result: EnsembleV2ResultJson {
            conformers,
            conformer_provenance,
            attempts: r.attempts.iter().map(conformer_attempt_json).collect(),
            mixed_force_field: r.mixed_force_field,
            termination: snake_case_debug(&r.termination),
            requested_count: r.requested_count,
        },
    };
    serde_json::to_string(&envelope).expect("success envelope must always serialize")
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Run `embed_ensemble_v2` on `mol`'s own atom order (never canonicalizes/
/// reparses, same convention as `embed_pipeline_v2_json`). See the module doc
/// for the config JSON shape and the success/failure envelope shape.
///
/// Never throws. Always returns a JSON string tagged `schemaVersion: 1` and
/// `ok: true`/`false`. `ok: false` covers both a WASM-level rejection
/// (oversized input, too many atoms, malformed/incomplete config JSON) and a
/// config `embed_ensemble_v2` itself rejects (currently only an invalid
/// `rmsdThreshold`) -- distinguished by `error.stage` /
/// `error.cause.kind`, same as `embed_pipeline_v2_json`.
#[wasm_bindgen]
pub fn embed_ensemble_v2_json(mol: &MolHandle, config_json: &str) -> String {
    if config_json.len() > WASM_MAX_INPUT_BYTES {
        return wasm_input_error_json(FailureCauseJson::InputTooLarge {
            limit_bytes: WASM_MAX_INPUT_BYTES,
            actual_bytes: config_json.len(),
        });
    }
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return wasm_input_error_json(FailureCauseJson::AtomLimitExceeded {
            limit: WASM_MAX_ATOMS,
            actual: mol.inner.atom_count(),
        });
    }

    let config = match parse_ensemble_config(config_json) {
        Ok(c) => c,
        Err(message) => {
            return wasm_input_error_json(FailureCauseJson::InvalidConfig { message });
        }
    };

    match embed_ensemble_v2(&mol.inner, &config) {
        Ok(result) => result_to_json(&result),
        Err(e) => wasm_input_error_json(FailureCauseJson::InvalidConfig {
            message: e.to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_depict::parse_smiles;

    fn safe_config_json(count: usize, base_seed: u64, ensemble_timeout_ms: &str) -> String {
        format!(
            r#"{{
                "perConformer": {{
                    "embedSeed": 7,
                    "maxAttempts": 8,
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
                    "totalTimeoutMs": null
                }},
                "count": {count},
                "baseSeed": {base_seed},
                "rmsdThreshold": 0.5,
                "useSymmetricRmsdPruning": true,
                "ensembleTimeoutMs": {ensemble_timeout_ms}
            }}"#
        )
    }

    #[test]
    fn success_path_has_expected_envelope_shape() {
        let mol = parse_smiles("CCCCCCCCCC").expect("decane");
        let config = safe_config_json(3, 20260828, "null");
        let json = embed_ensemble_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], true);
        let result = &value["result"];
        for key in [
            "conformers",
            "conformerProvenance",
            "attempts",
            "mixedForceField",
            "termination",
            "requestedCount",
        ] {
            assert!(result.get(key).is_some(), "missing result.{key}");
        }
        assert_eq!(result["requestedCount"], 3);
        assert_eq!(result["attempts"].as_array().unwrap().len(), 3);
        assert_eq!(result["termination"], "completed");
        let conformers = result["conformers"].as_array().unwrap();
        assert!(!conformers.is_empty(), "decane should embed at least once");
        assert_eq!(
            conformers.len(),
            result["conformerProvenance"].as_array().unwrap().len()
        );
        for conformer in conformers {
            assert_eq!(conformer.as_array().unwrap().len(), 10, "decane atom count");
        }
    }

    #[test]
    fn same_base_seed_is_reproducible() {
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").expect("aspirin");
        let config = safe_config_json(4, 42, "null");
        let a = embed_ensemble_v2_json(&mol, &config);
        let b = embed_ensemble_v2_json(&mol, &config);
        let av: serde_json::Value = serde_json::from_str(&a).unwrap();
        let bv: serde_json::Value = serde_json::from_str(&b).unwrap();
        assert_eq!(
            av["result"]["conformers"], bv["result"]["conformers"],
            "same base_seed must reproduce identical conformers"
        );
    }

    #[test]
    fn zero_count_returns_ok_with_empty_ensemble() {
        let mol = parse_smiles("CC").expect("ethane");
        let config = safe_config_json(0, 1, "null");
        let json = embed_ensemble_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
        assert_eq!(value["result"]["conformers"].as_array().unwrap().len(), 0);
        assert_eq!(value["result"]["attempts"].as_array().unwrap().len(), 0);
        assert_eq!(value["result"]["termination"], "completed");
    }

    #[test]
    fn invalid_rmsd_threshold_rejected_as_typed_error_not_panic() {
        let mol = parse_smiles("CC").expect("ethane");
        let config = r#"{
                "perConformer": {
                    "embedSeed": 7, "maxAttempts": 8, "embedTimeoutMs": null,
                    "useExpTorsions": false, "useSmallRingTorsions": false,
                    "useMacrocycleTorsions": false, "useMacrocycle14Bounds": false,
                    "includeLegacyTorsionHeuristic": false, "stereoPolicy": "ignore",
                    "failOnUnevaluableStereo": false, "forceFieldPolicy": "none",
                    "forceFieldMaxIterations": 200, "gateMmff94TorsionOop": false,
                    "gateMmff94StretchBend": false, "ringTorsionPolicy": "fail_closed",
                    "totalTimeoutMs": null
                },
                "count": 2, "baseSeed": 1, "rmsdThreshold": -1.0,
                "useSymmetricRmsdPruning": true, "ensembleTimeoutMs": null
            }"#;
        let json = embed_ensemble_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["stage"], "wasm_input_validation");
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn malformed_and_incomplete_config_json_rejected_never_throws() {
        let mol = parse_smiles("CC").expect("ethane");
        let cases = [
            "{not valid json",
            r#"{"perConformer": {}, "count": 1, "baseSeed": 1, "rmsdThreshold": 0.5, "useSymmetricRmsdPruning": true, "ensembleTimeoutMs": null, "notARealField": true}"#,
        ];
        for config in cases {
            let json = embed_ensemble_v2_json(&mol, config);
            let value: serde_json::Value =
                serde_json::from_str(&json).expect("must be valid JSON, never throw");
            assert_eq!(value["ok"], false);
            assert_eq!(value["error"]["stage"], "wasm_input_validation");
            assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
        }
    }

    #[test]
    fn missing_present_but_nullable_ensemble_timeout_field_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let config = r#"{
            "perConformer": {
                "embedSeed": 7, "maxAttempts": 8, "embedTimeoutMs": null,
                "useExpTorsions": false, "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false, "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false, "stereoPolicy": "ignore",
                "failOnUnevaluableStereo": false, "forceFieldPolicy": "none",
                "forceFieldMaxIterations": 200, "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false, "ringTorsionPolicy": "fail_closed",
                "totalTimeoutMs": null
            },
            "count": 1, "baseSeed": 1, "rmsdThreshold": 0.5,
            "useSymmetricRmsdPruning": true
        }"#;
        let json = embed_ensemble_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn every_attempt_failing_is_still_ok_true_with_zero_kept_conformers() {
        // Mirrors the module doc's documented error asymmetry: a config that
        // makes every attempt fail (here: a near-zero per-attempt timeout) is
        // still a normal `ok: true` ensemble result with zero kept conformers,
        // never surfaced as a top-level failure.
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").expect("aspirin");
        let config = r#"{
                "perConformer": {
                    "embedSeed": 7, "maxAttempts": 8, "embedTimeoutMs": null,
                    "useExpTorsions": false, "useSmallRingTorsions": false,
                    "useMacrocycleTorsions": false, "useMacrocycle14Bounds": false,
                    "includeLegacyTorsionHeuristic": false, "stereoPolicy": "ignore",
                    "failOnUnevaluableStereo": false, "forceFieldPolicy": "none",
                    "forceFieldMaxIterations": 200, "gateMmff94TorsionOop": false,
                    "gateMmff94StretchBend": false, "ringTorsionPolicy": "fail_closed",
                    "totalTimeoutMs": 0
                },
                "count": 2, "baseSeed": 1, "rmsdThreshold": 0.5,
                "useSymmetricRmsdPruning": true, "ensembleTimeoutMs": null
            }"#;
        let json = embed_ensemble_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "an all-failed ensemble is still ok:true");
        assert_eq!(value["result"]["conformers"].as_array().unwrap().len(), 0);
        let attempts = value["result"]["attempts"].as_array().unwrap();
        assert_eq!(attempts.len(), 2);
        for attempt in attempts {
            assert_eq!(attempt["outcome"], "failure");
            assert!(attempt["failure"].is_object());
            assert_eq!(attempt["success"], serde_json::Value::Null);
        }
    }

    #[test]
    fn oversized_config_json_rejected_before_parsing() {
        let mol = parse_smiles("CC").expect("ethane");
        let padded = " ".repeat(WASM_MAX_INPUT_BYTES + 1);
        let json = embed_ensemble_v2_json(&mol, &padded);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "input_too_large");
    }
}
