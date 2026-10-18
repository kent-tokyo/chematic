//! Emit Cleanup results using the RDKit-compatible SMILES writer.

use std::io::{self, BufRead, Write};

use chematic_chem::{cleanup, normalize_groups, reionize};
use serde_json::json;

fn main() {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let line = line.expect("read SMILES row");
        let Some(smiles) = line.split_whitespace().next() else {
            continue;
        };
        let record = match chematic_smiles::parse(smiles) {
            Ok(mol) => json!({
                "index": index,
                "input": smiles,
                "rdkit_input": chematic_smiles::rdkit_canonical_smiles(&mol).map_err(|error| error.to_string()),
                "normalize": chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&mol)).map_err(|error| error.to_string()),
                "cleanup": chematic_smiles::rdkit_cleanup_canonical_smiles(&cleanup(&mol)).map_err(|error| error.to_string()),
                "reionize": chematic_smiles::rdkit_canonical_smiles(&reionize(&mol)).map_err(|error| error.to_string()),
            }),
            Err(error) => json!({"index": index, "input": smiles, "error": error.to_string()}),
        };
        serde_json::to_writer(&mut out, &record).expect("serialize result");
        writeln!(out).expect("write result");
    }
}
