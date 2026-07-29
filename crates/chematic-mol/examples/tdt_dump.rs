//! IO-2 acceptance-gate dump: runs `TdtRecordReader` over every scenario
//! file named in `scripts/gen_tdt_fixtures.py`'s manifest and dumps each
//! row's extracted `(status, name, properties, coordinates)` as JSONL for
//! comparison against a real RDKit oracle
//! (`scripts/gen_rdkit_tdt_oracle.py`).
//!
//! Usage:
//! ```text
//! cargo run -p chematic-mol --release --example tdt_dump -- \
//!     <manifest.json> <fixtures_dir> <out.jsonl>
//! ```

use chematic_mol::{TdtError, TdtReaderOptions, TdtRecordReader};
use serde_json::{Value, json};
use std::fs;
use std::io::{BufReader, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let manifest_path = args
        .get(1)
        .expect("usage: tdt_dump <manifest.json> <fixtures_dir> <out.jsonl>");
    let fixtures_dir = args
        .get(2)
        .expect("usage: tdt_dump <manifest.json> <fixtures_dir> <out.jsonl>");
    let out_path = args
        .get(3)
        .expect("usage: tdt_dump <manifest.json> <fixtures_dir> <out.jsonl>");

    let manifest_text =
        fs::read_to_string(manifest_path).unwrap_or_else(|e| panic!("read manifest: {e}"));
    let manifest: Value =
        serde_json::from_str(&manifest_text).unwrap_or_else(|e| panic!("parse manifest: {e}"));

    let mut out = fs::File::create(out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    let mut total_rows = 0usize;

    let scenarios = manifest["scenarios"].as_object().expect("scenarios object");
    let mut scenario_names: Vec<&String> = scenarios.keys().collect();
    scenario_names.sort();

    for name in scenario_names {
        let scenario = &scenarios[name];
        let file_name = scenario["file"].as_str().unwrap();
        let path = format!("{fixtures_dir}/{file_name}");
        let file = fs::File::open(&path).unwrap_or_else(|e| panic!("open {path}: {e}"));

        let known_rows = scenario["rows"].as_array().cloned().unwrap_or_default();

        let read_coords = name == "coordinates";
        let reader_options = TdtReaderOptions {
            name_tag: Some("NAME".to_string()),
            read_2d: read_coords,
            read_3d: read_coords,
            strict_parsing: false,
            ..Default::default()
        };
        let reader = TdtRecordReader::new(BufReader::new(file), reader_options);

        for (row_index, result) in reader.enumerate() {
            let row = match result {
                Ok(rec) => {
                    let mut props: Vec<(String, String)> = rec.properties;
                    props.sort_unstable();

                    let self_consistent = known_rows
                        .get(row_index)
                        .and_then(|r| r["smiles"].as_str())
                        .map(|known_smi| match chematic_smiles::parse(known_smi) {
                            Ok(known_mol) => {
                                chematic_smiles::write(&known_mol)
                                    == chematic_smiles::write(&rec.mol)
                            }
                            Err(_) => false,
                        });

                    json!({
                        "scenario": name,
                        "row_index": row_index,
                        "status": "success",
                        "name": rec.name,
                        "properties": props,
                        "coordinates_2d": rec.coordinates_2d,
                        "coordinates_3d": rec.coordinates_3d,
                        "self_consistent_with_known_smiles": self_consistent,
                    })
                }
                Err(TdtError::Io(msg)) => panic!("unexpected IO error in {name}: {msg}"),
                Err(e) => json!({
                    "scenario": name,
                    "row_index": row_index,
                    "status": "error",
                    "error": e.to_string(),
                }),
            };
            writeln!(out, "{row}").unwrap();
            total_rows += 1;
        }
    }

    eprintln!("total_rows={total_rows} out={out_path}");
}
