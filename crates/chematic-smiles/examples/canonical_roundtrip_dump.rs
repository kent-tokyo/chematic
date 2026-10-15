//! Print the native canonical SMILES for every line of a SMILES corpus.
//!
//! Usage: `cargo run -p chematic-smiles --example canonical_roundtrip_dump -- corpus.smi`

use std::io::BufRead;

fn main() {
    let path = std::env::args().nth(1).expect("pass a SMILES corpus");
    let file = std::fs::File::open(path).expect("open corpus");
    for (index, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line.expect("read corpus");
        let Some(smiles) = line.split_whitespace().next() else {
            continue;
        };
        match chematic_smiles::parse(smiles) {
            Ok(mol) => println!("{index}\t{}", chematic_smiles::canonical_smiles(&mol)),
            Err(error) => println!("{index}\tERROR:{error}"),
        }
    }
}
