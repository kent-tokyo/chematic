//! Emit `smiles\tstatus\tname` for every line of a SMILES corpus, where
//! `status` is `PARSE_FAIL`, `OK`, or `NOT_SUPPORTED` (`name` is empty for
//! the latter two). Used to diff `chematic_iupac::name()` behavior
//! before/after a change against a fixed corpus -- see
//! `scripts/iupac_snapshot_diff.py`.
//!
//! Usage: cargo run -p chematic-iupac --release --example iupac_name_snapshot -- <SMILES.csv> <out.tsv>

use std::env;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = env::args().collect();
    let in_path = args.get(1).expect("usage: <SMILES.csv> <out.tsv>");
    let out_path = args.get(2).expect("usage: <SMILES.csv> <out.tsv>");

    let input = fs::read_to_string(in_path).expect("read input csv");
    let mut out = fs::File::create(out_path).expect("create output tsv");

    for line in input.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let row = match chematic_smiles::parse(smi) {
            Err(_) => format!("{smi}\tPARSE_FAIL\t"),
            Ok(mol) => match chematic_iupac::name(&mol) {
                Ok(n) => format!("{smi}\tOK\t{n}"),
                Err(_) => format!("{smi}\tNOT_SUPPORTED\t"),
            },
        };
        writeln!(out, "{row}").expect("write row");
    }
}
