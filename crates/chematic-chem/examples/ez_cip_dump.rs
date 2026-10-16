//! Dump source-built E/Z CIP labels for comparison with an external oracle.
//! Usage: cargo run -p chematic-chem --example ez_cip_dump -- corpus.smi

use std::io::BufRead;

use chematic_chem::{CipMode, assign_ez_bonds_with_mode};

fn main() {
    let path = std::env::args().nth(1).expect("pass a SMILES corpus");
    let file = std::fs::File::open(path).expect("open corpus");
    for (index, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line.expect("read corpus");
        let Some(smiles) = line.split_whitespace().next() else {
            continue;
        };
        let Ok(mol) = chematic_smiles::parse(smiles) else {
            println!("{index}\tERROR");
            continue;
        };
        for mode in [CipMode::LegacyFast, CipMode::Accurate] {
            for (bond_idx, code) in assign_ez_bonds_with_mode(&mol, mode) {
                let bond = mol.bond(bond_idx);
                println!(
                    "{index}\t{mode:?}\t{}\t{}\t{code:?}",
                    bond.atom1.0, bond.atom2.0
                );
            }
        }
    }
}
