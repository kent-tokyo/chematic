//! One atom-indexed fingerprint result per SMILES on stdin for the RDKit ring audit.

use std::io::{self, BufRead};

fn main() {
    for line in io::stdin().lock().lines() {
        let smiles = line.expect("read SMILES");
        let mol = chematic_smiles::parse(&smiles).expect("parse SMILES");
        let sparse = chematic_fp::rdkit_torsion_sparse_counts(&mol);
        let bitvec = chematic_fp::rdkit_torsion_fp(&mol);
        let bits: Vec<_> = (0..2048).filter(|&bit| bitvec.get(bit)).collect();
        println!("{smiles}\t{sparse:?}\t{bits:?}");
    }
}
