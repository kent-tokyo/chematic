//! EZ-S1: snapshot `assign_cip`'s (legacy engine) full R/S/E/Z output across a
//! SMILES corpus, for a before/after diff of the `substituent_is_up` fix that
//! makes `assign_ez`/`assign_allene` read `Molecule::bond_direction` (the
//! stashed `/`/`\` side channel used when a marker lands on a bond between
//! two aromatic-flagged atoms, e.g. a ring bond flanking an exocyclic C=N).
//!
//! For E/Z rows, `partner_idx` (the double bond's other atom) is also
//! emitted so a downstream diff can key against RDKit's bond-level
//! `_CIPCode` oracle; it's empty for R/S/r/s rows.
//!
//! Run once on `main` (before the fix) and once on the fix branch, diff the
//! two TSVs:
//! ```text
//! cargo run -p chematic-chem --release --example ez_stash_gap_snapshot \
//!     -- ~/Downloads/SMILES.csv baseline.tsv
//! ```

use chematic_chem::assign_cip;
use chematic_core::{AtomIdx, BondOrder, CipCode};
use chematic_smiles::parse;
use std::fs;
use std::io::Write;

fn code_str(c: CipCode) -> &'static str {
    match c {
        CipCode::R => "R",
        CipCode::S => "S",
        CipCode::E => "E",
        CipCode::Z => "Z",
        CipCode::LowerR => "r",
        CipCode::LowerS => "s",
    }
}

/// The other atom of the double bond at `idx`, if any (used to key E/Z rows
/// against RDKit's bond-level oracle downstream).
fn double_bond_partner(mol: &chematic_core::Molecule, idx: AtomIdx) -> Option<AtomIdx> {
    mol.neighbors(idx)
        .find(|(_, bidx)| mol.bond(*bidx).order == BondOrder::Double)
        .map(|(nb, _)| nb)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: ez_stash_gap_snapshot <SMILES.csv> <out.tsv>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "snapshot.tsv".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));

    let mut lines: Vec<String> = Vec::new();
    let mut parse_fail = 0usize;
    let mut rows_written = 0usize;

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

        for (idx, code) in &assign_cip(&mol).assignments {
            let partner = match code {
                CipCode::E | CipCode::Z => double_bond_partner(&mol, *idx)
                    .map(|p| p.0.to_string())
                    .unwrap_or_default(),
                _ => String::new(),
            };
            lines.push(format!("{smi}\t{}\t{partner}\t{}", idx.0, code_str(*code)));
            rows_written += 1;
        }
    }

    lines.sort();
    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    for l in &lines {
        writeln!(f, "{l}").unwrap();
    }

    eprintln!(
        "input_lines={} parse_fail={parse_fail} rows_written={rows_written} out={out_path}",
        content.lines().filter(|l| !l.trim().is_empty()).count(),
    );
}
