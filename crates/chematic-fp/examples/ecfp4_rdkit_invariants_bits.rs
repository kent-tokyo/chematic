//! ECFP RDKit-invariant-mode PR: reference (non-gating) fingerprint-level
//! measurement input. Dumps `ecfp4_rdkit_invariants(mol)`'s set-bit list for
//! every corpus molecule, for a downstream Tanimoto-correlation-vs-RDKit
//! comparison (`scripts/ecfp_rdkit_invariant_parity_fingerprint_ref.py`).
//! This is reference data only -- not part of the atom-invariant-partition
//! merge gate.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example ecfp4_rdkit_invariants_bits \
//!     -- <SMILES.csv> <out.tsv>
//! ```

use chematic_fp::ecfp4_rdkit_invariants;
use chematic_smiles::parse;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: ecfp4_rdkit_invariants_bits <SMILES.csv> <out.tsv>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "snapshot.tsv".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));

    let mut lines: Vec<String> = Vec::new();
    let mut parse_fail = 0usize;
    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(_) => {
                parse_fail += 1;
                continue;
            }
        };
        let fp = ecfp4_rdkit_invariants(&mol);
        let bits: Vec<String> = (0..2048)
            .filter(|&i| fp.get(i))
            .map(|i| i.to_string())
            .collect();
        lines.push(format!("{smi}\t{}", bits.join(",")));
    }

    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    for l in &lines {
        writeln!(f, "{l}").unwrap();
    }

    eprintln!(
        "input_lines={} parse_fail={parse_fail} rows={} out={out_path}",
        content.lines().filter(|l| !l.trim().is_empty()).count(),
        lines.len(),
    );
}
