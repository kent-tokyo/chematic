//! Canonical-Stereo-D0: snapshot each molecule's FIRST-PARSE legacy CIP
//! element set, keyed by input line number. Used to diff old-vs-new parser
//! behavior on a real corpus and confirm the `resolve_aromatic_direction_stash`
//! fix has zero effect on first-parse output -- it only changes behavior on
//! re-parsing chematic's own canonical output (ring-closure and
//! branch-attachment bond paths, previously missing the guard the plain
//! chain-edge path already had).
//!
//! Verification recipe used for the D0 PR: run this once on `main` (before
//! the fix), once on the fix branch, diff the two outputs -- 5,000/5,000
//! lines were byte-identical.
//!
//! Run:
//! ```text
//! cargo run -p chematic-chem --release --example canonical_stereo_d0_corpus_snapshot \
//!     -- ~/Downloads/SMILES.csv > /tmp/d0_snapshot.tsv
//! ```

use chematic_chem::assign_cip;
use chematic_smiles::parse;
use std::fs;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| panic!("usage: canonical_stereo_d0_corpus_snapshot <smiles_file>"));
    let content = fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));

    for (i, line) in content.lines().enumerate() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                println!("{i}\tPARSE_FAIL\t{e}");
                continue;
            }
        };
        let mut codes: Vec<String> = assign_cip(&mol)
            .assignments
            .iter()
            .map(|(_, c)| format!("{c:?}"))
            .collect();
        codes.sort();
        println!("{i}\t{}\t{}", codes.join(","), smi);
    }
}
