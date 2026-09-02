//! Cross-language shared-corpus test for the promoted RDKit-exact Morgan/ECFP API.
//!
//! Reads `validation/ecfp4_rdkit_stable_api_fixtures.json` -- generated once from a
//! live RDKit oracle by `scripts/gen_ecfp4_rdkit_stable_api_fixtures.py` (RDKit
//! version/commit recorded in the file itself) -- and checks
//! `rdkit_morgan_ecfp4_experimental`/`rdkit_morgan_fingerprint`'s output against it
//! exactly. The Python (`crates/chematic-py/tests/test_rdkit_ecfp4_stable_api.py`) and
//! WASM (`crates/chematic-wasm/tests/rdkit_ecfp4_stable_api.test.mjs`) test suites read
//! this *same* file, so all three surfaces are checked against one shared source of
//! truth rather than three independently-maintained expectation lists.

use chematic_fp::{RdkitMorganConfig, RdkitMorganError, RdkitMorganFpSize, RdkitMorganRadius};
use chematic_smiles::parse;
use serde_json::Value;
use std::collections::BTreeMap;

const CORPUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/ecfp4_rdkit_stable_api_fixtures.json"
));

fn radius_from_int(r: u64) -> RdkitMorganRadius {
    match r {
        0 => RdkitMorganRadius::R0,
        1 => RdkitMorganRadius::R1,
        2 => RdkitMorganRadius::R2,
        3 => RdkitMorganRadius::R3,
        other => panic!("unsupported radius in fixture corpus: {other}"),
    }
}

fn fp_size_from_int(n: u64) -> RdkitMorganFpSize {
    match n {
        128 => RdkitMorganFpSize::B128,
        256 => RdkitMorganFpSize::B256,
        512 => RdkitMorganFpSize::B512,
        1024 => RdkitMorganFpSize::B1024,
        2048 => RdkitMorganFpSize::B2048,
        other => panic!("unsupported fp_size in fixture corpus: {other}"),
    }
}

fn json_u32_map(v: &Value) -> BTreeMap<u32, u32> {
    v.as_object()
        .unwrap()
        .iter()
        .map(|(k, val)| (k.parse::<u32>().unwrap(), val.as_u64().unwrap() as u32))
        .collect()
}

fn json_bit_info_map(v: &Value) -> BTreeMap<u32, Vec<(u32, u32)>> {
    v.as_object()
        .unwrap()
        .iter()
        .map(|(k, val)| {
            let mut pairs: Vec<(u32, u32)> = val
                .as_array()
                .unwrap()
                .iter()
                .map(|pair| {
                    let a = pair[0].as_u64().unwrap() as u32;
                    let r = pair[1].as_u64().unwrap() as u32;
                    (a, r)
                })
                .collect();
            pairs.sort_unstable();
            (k.parse::<u32>().unwrap(), pairs)
        })
        .collect()
}

fn json_bit_list(v: &Value) -> Vec<usize> {
    let mut bits: Vec<usize> = v
        .as_array()
        .unwrap()
        .iter()
        .map(|b| b.as_u64().unwrap() as usize)
        .collect();
    bits.sort_unstable();
    bits
}

#[test]
fn rdkit_morgan_ecfp4_experimental_matches_shared_rdkit_oracle_corpus() {
    let doc: Value = serde_json::from_str(CORPUS).expect("corpus JSON must parse");
    let rdkit_version = doc["rdkit_version"].as_str().unwrap();
    assert!(
        !rdkit_version.is_empty(),
        "corpus must record the RDKit version used to generate it"
    );

    let fixtures = doc["fixtures"].as_array().unwrap();
    assert!(
        fixtures.len() >= 30,
        "expected a real fixture corpus, not a stub"
    );

    let mut checked_ok = 0usize;
    let mut checked_error = 0usize;

    for fx in fixtures {
        let id = fx["id"].as_str().unwrap();
        let smiles = fx["smiles"].as_str().unwrap();
        let expect = fx["expect"].as_str().unwrap();
        let mol = parse(smiles).unwrap_or_else(|e| panic!("fixture {id} ({smiles}): {e}"));

        match expect {
            "ok" => {
                let result = chematic_fp::rdkit_morgan_ecfp4_experimental(&mol)
                    .unwrap_or_else(|e| panic!("fixture {id} ({smiles}) expected ok, got {e}"));

                let got_bits: Vec<usize> = (0..2048usize)
                    .filter(|&b| result.fingerprint.get(b))
                    .collect();
                assert_eq!(
                    got_bits,
                    json_bit_list(&fx["folded_bits"]),
                    "folded_bits mismatch for fixture {id}"
                );

                let got_counts: BTreeMap<u32, u32> =
                    result.sparse_counts.iter().map(|(&k, &v)| (k, v)).collect();
                assert_eq!(
                    got_counts,
                    json_u32_map(&fx["sparse_counts"]),
                    "sparse_counts mismatch for fixture {id}"
                );

                let got_raw_bit_info: BTreeMap<u32, Vec<(u32, u32)>> = result
                    .raw_bit_info
                    .iter()
                    .map(|(&k, v)| {
                        let mut pairs = v.clone();
                        pairs.sort_unstable();
                        (k, pairs)
                    })
                    .collect();
                assert_eq!(
                    got_raw_bit_info,
                    json_bit_info_map(&fx["raw_bit_info"]),
                    "raw_bit_info mismatch for fixture {id}"
                );

                let got_folded_bit_info: BTreeMap<u32, Vec<(u32, u32)>> = result
                    .folded_bit_info
                    .iter()
                    .map(|(&k, v)| {
                        let mut pairs = v.clone();
                        pairs.sort_unstable();
                        (k as u32, pairs)
                    })
                    .collect();
                assert_eq!(
                    got_folded_bit_info,
                    json_bit_info_map(&fx["folded_bit_info"]),
                    "folded_bit_info mismatch for fixture {id}"
                );

                checked_ok += 1;
            }
            "error" => {
                let error_kind = fx["error_kind"].as_str().unwrap();
                let result = chematic_fp::rdkit_morgan_ecfp4_experimental(&mol);
                match (error_kind, &result) {
                    ("Aromaticity", Err(RdkitMorganError::Aromaticity(_))) => {}
                    _ => panic!(
                        "fixture {id} ({smiles}) expected error kind {error_kind}, got {result:?}"
                    ),
                }
                checked_error += 1;
            }
            other => panic!("fixture {id}: unknown expect '{other}'"),
        }
    }

    assert!(
        checked_ok >= 30,
        "expected real success coverage, got {checked_ok}"
    );
    assert!(
        checked_error >= 1,
        "expected real error-path coverage, got {checked_error}"
    );
}

#[test]
fn rdkit_morgan_fingerprint_matches_shared_oracle_corpus_across_radius_and_fp_size_axes() {
    let doc: Value = serde_json::from_str(CORPUS).expect("corpus JSON must parse");
    let fixtures = doc["fixtures"].as_array().unwrap();

    let mut radius_cells_checked = 0usize;
    let mut fp_size_cells_checked = 0usize;

    for fx in fixtures {
        if fx["expect"].as_str().unwrap() != "ok" {
            continue;
        }
        let id = fx["id"].as_str().unwrap();
        let smiles = fx["smiles"].as_str().unwrap();
        let mol = parse(smiles).unwrap();

        for cell in fx["radius_axis"].as_array().unwrap() {
            let radius_int = cell["radius"].as_u64().unwrap();
            let config = RdkitMorganConfig {
                radius: radius_from_int(radius_int),
                fp_size: RdkitMorganFpSize::B2048,
                include_chirality: false,
            };
            let result = chematic_fp::rdkit_morgan_fingerprint(&mol, &config)
                .unwrap_or_else(|e| panic!("fixture {id} radius={radius_int}: {e}"));
            let got_bits: Vec<usize> = (0..2048usize)
                .filter(|&b| result.fingerprint.get(b))
                .collect();
            assert_eq!(
                got_bits,
                json_bit_list(&cell["folded_bits"]),
                "radius axis mismatch for fixture {id} at radius={radius_int}"
            );
            radius_cells_checked += 1;
        }

        for cell in fx["fp_size_axis"].as_array().unwrap() {
            let fp_size_int = cell["fp_size"].as_u64().unwrap();
            let fp_size = fp_size_from_int(fp_size_int);
            let config = RdkitMorganConfig {
                radius: RdkitMorganRadius::R2,
                fp_size,
                include_chirality: false,
            };
            let result = chematic_fp::rdkit_morgan_fingerprint(&mol, &config)
                .unwrap_or_else(|e| panic!("fixture {id} fp_size={fp_size_int}: {e}"));
            let got_bits: Vec<usize> = (0..fp_size.bits())
                .filter(|&b| result.fingerprint.get(b))
                .collect();
            assert_eq!(
                got_bits,
                json_bit_list(&cell["folded_bits"]),
                "fp_size axis mismatch for fixture {id} at fp_size={fp_size_int}"
            );
            fp_size_cells_checked += 1;
        }
    }

    assert!(
        radius_cells_checked >= 30 * 4,
        "expected full radius-axis coverage"
    );
    assert!(
        fp_size_cells_checked >= 30 * 5,
        "expected full fp_size-axis coverage"
    );
}
