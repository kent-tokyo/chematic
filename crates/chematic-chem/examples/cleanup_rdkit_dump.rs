//! Emit Cleanup results using the RDKit-compatible SMILES writer.

use std::io::{self, BufRead, Write};

use chematic_chem::{cleanup_unchecked, normalize_groups, reionize};
use serde_json::json;

fn main() {
    let legacy = std::env::args().skip(1).any(|arg| arg == "--legacy");
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let line = line.expect("read SMILES row");
        let Some(smiles) = line.split_whitespace().next() else {
            continue;
        };
        let record = match chematic_smiles::parse(smiles) {
            Ok(mol) => {
                let cleaned = if legacy {
                    Ok(cleanup_unchecked(&mol))
                } else {
                    chematic_chem::cleanup(&mol)
                        .map_err(|error| format!("cleanup_metal_hydrogen: {error}"))
                };
                let cleanup_smiles = cleaned.and_then(|cleaned| {
                    chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned)
                        .map_err(|error| error.to_string())
                });
                json!({
                "index": index,
                "input": smiles,
                "rdkit_input": chematic_smiles::rdkit_canonical_smiles(&mol).map_err(|error| error.to_string()),
                "normalize": chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&mol)).map_err(|error| error.to_string()),
                "cleanup": cleanup_smiles,
                "reionize": chematic_smiles::rdkit_canonical_smiles(&reionize(&mol)).map_err(|error| error.to_string()),
                })
            }
            Err(error) => json!({"index": index, "input": smiles, "error": error.to_string()}),
        };
        serde_json::to_writer(&mut out, &record).expect("serialize result");
        writeln!(out).expect("write result");
    }
}
