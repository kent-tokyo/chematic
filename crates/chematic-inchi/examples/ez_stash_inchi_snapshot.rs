//! EZ-S1 sibling gap: snapshot `standard_inchi()` (native-inchi feature)
//! output across a SMILES corpus, for a before/after diff of
//! `find_stereo_sub`'s fix (`crates/chematic-inchi/src/native/convert.rs`) --
//! the same stashed-aromatic-direction gap fixed in `chematic-chem`'s
//! `substituent_is_up` (see `crates/chematic-chem/src/cip.rs`), applied to
//! InChI `Stereo0D` double-bond descriptors.
//!
//! One row per line: `smiles\tinchi_or_error`.
//!
//! Run once on `main` (before the fix) and once on the fix branch, diff the
//! two TSVs:
//! ```text
//! cargo run -p chematic-inchi --release --features native-inchi \
//!     --example ez_stash_inchi_snapshot -- ~/Downloads/SMILES.csv baseline.tsv
//! ```

use chematic_inchi::standard_inchi;
use chematic_smiles::parse;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: ez_stash_inchi_snapshot <SMILES.csv> <out.tsv>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "snapshot.tsv".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));

    let mut lines: Vec<String> = Vec::new();
    let mut smiles_parse_fail = 0usize;
    let mut inchi_fail = 0usize;

    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(_) => {
                smiles_parse_fail += 1;
                continue;
            }
        };
        let out = match standard_inchi(&mol) {
            Ok(s) => s,
            Err(e) => {
                inchi_fail += 1;
                format!("ERR:{e:?}")
            }
        };
        lines.push(format!("{smi}\t{out}"));
    }

    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    for l in &lines {
        writeln!(f, "{l}").unwrap();
    }

    eprintln!(
        "input_lines={} smiles_parse_fail={smiles_parse_fail} inchi_fail={inchi_fail} rows={} out={out_path}",
        content.lines().filter(|l| !l.trim().is_empty()).count(),
        lines.len(),
    );
}
