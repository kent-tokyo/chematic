//! ECFP RDKit-invariant-mode PR: per-atom `EcfpInvariantMode::RdkitMorgan`
//! invariant dump, for a partition-agreement comparison against RDKit's
//! `GetConnectivityInvariants` (see `scripts/ecfp_rdkit_invariant_parity.py`).
//!
//! One row per atom: `smiles\tatom_idx\tinvariant`. Raw invariant values are
//! NOT meant to be compared cross-implementation directly (chematic uses
//! FNV-1a, RDKit uses its own hash) -- only whether two atoms get the SAME
//! value as each other (the equivalence partition) is a meaningful
//! comparison.
//!
//! A SMILES chematic fails to parse gets a marker row instead of being
//! silently dropped: `smiles\tPARSE_FAIL\t<error>` -- so a downstream
//! comparator can explicitly count and gate on chematic parse completeness,
//! not just eyeball this tool's own stderr summary.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example rdkit_invariant_snapshot \
//!     -- <SMILES.csv> <out.tsv>
//! ```

use chematic_fp::{EcfpInvariantMode, atom_invariants};
use chematic_smiles::parse;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: rdkit_invariant_snapshot <SMILES.csv> <out.tsv>"));
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
            Err(e) => {
                parse_fail += 1;
                lines.push(format!("{smi}\tPARSE_FAIL\t{e}"));
                continue;
            }
        };
        for (idx, inv) in atom_invariants(&mol, EcfpInvariantMode::RdkitMorgan)
            .into_iter()
            .enumerate()
        {
            lines.push(format!("{smi}\t{idx}\t{inv}"));
        }
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
