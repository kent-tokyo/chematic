//! EZ-S1: snapshot `assign_cip`'s (legacy engine) full R/S/E/Z output across a
//! SMILES corpus, for a before/after diff of the `substituent_is_up` fix that
//! makes `assign_ez`/`assign_allene` read `Molecule::bond_direction` (the
//! stashed `/`/`\` side channel used when a marker lands on a bond between
//! two aromatic-flagged atoms, e.g. a ring bond flanking an exocyclic C=N),
//! and the `highest_stereo_sub` fix (true-highest-priority substituent, not
//! the highest-priority-among-marked one).
//!
//! Each row is `smiles\tkind\tatom_idx\tpartner_idx\tcode`:
//! - `kind` is `tetra` (R/S/r/s), `ez` (a plain double bond's E/Z), `allene`
//!   (axial chirality's E/Z), or `parse_fail` (chematic couldn't parse this
//!   SMILES at all -- `code` carries the error message, `atom_idx`/
//!   `partner_idx` are empty).
//! - `ez` vs `allene` is a structural distinction, not the enum discriminant
//!   (`CipCode::E`/`Z` is shared by both): an E/Z atom's double-bond partner
//!   is checked for the same "allene central atom" shape `cip.rs`'s
//!   `is_allene_central` uses (exactly 2 double bonds, exactly 2 heavy-atom
//!   neighbors) -- if it matches, this row's atom is an allene terminal
//!   reporting through its bond to the central atom, so `kind=allene`;
//!   otherwise it's a plain stereogenic double bond, `kind=ez`. Reimplemented
//!   here rather than exposed from `chematic-chem`'s public API.
//! - `partner_idx` is the double bond's other atom, non-empty only for `ez`/
//!   `allene` rows -- keys a downstream diff against RDKit's bond-level
//!   `_CIPCode` oracle via `(smiles, frozenset({atom_idx, partner_idx}))`.
//!
//! Run once on `main` (before the fix) and once on the fix branch, diff the
//! two TSVs:
//! ```text
//! cargo run -p chematic-chem --release --example ez_stash_gap_snapshot \
//!     -- ~/Downloads/SMILES.csv baseline.tsv
//! ```

use chematic_chem::assign_cip;
use chematic_core::{AtomIdx, BondOrder, CipCode, Molecule};
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

/// The other atom of the double bond at `idx`, if any.
fn double_bond_partner(mol: &Molecule, idx: AtomIdx) -> Option<AtomIdx> {
    mol.neighbors(idx)
        .find(|(_, bidx)| mol.bond(*bidx).order == BondOrder::Double)
        .map(|(nb, _)| nb)
}

/// Mirrors `chematic_chem::cip`'s private `is_allene_central`: exactly 2
/// double bonds and exactly 2 heavy-atom neighbors.
fn is_allene_central(mol: &Molecule, idx: AtomIdx) -> bool {
    let dbl_count = mol
        .neighbors(idx)
        .filter(|(_, bidx)| mol.bond(*bidx).order == BondOrder::Double)
        .count();
    dbl_count == 2 && mol.neighbors(idx).count() == 2
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
            Err(e) => {
                parse_fail += 1;
                lines.push(format!("{smi}\tparse_fail\t\t\t{e}"));
                continue;
            }
        };

        for (idx, code) in &assign_cip(&mol).assignments {
            let (kind, partner) = match code {
                CipCode::E | CipCode::Z => match double_bond_partner(&mol, *idx) {
                    Some(p) if is_allene_central(&mol, p) => ("allene", p.0.to_string()),
                    Some(p) => ("ez", p.0.to_string()),
                    None => ("ez", String::new()),
                },
                _ => ("tetra", String::new()),
            };
            lines.push(format!(
                "{smi}\t{kind}\t{}\t{partner}\t{}",
                idx.0,
                code_str(*code)
            ));
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
