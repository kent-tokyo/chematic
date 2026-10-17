//! Print randomized SMILES for each row of a SMILES corpus.
//!
//! Usage: `random_roundtrip_dump corpus.smi [seed-count]`.

use std::io::BufRead;

fn main() {
    let path = std::env::args().nth(1).expect("pass a SMILES corpus");
    let seeds: u64 = std::env::args()
        .nth(2)
        .map(|value| value.parse().expect("seed count must be an integer"))
        .unwrap_or(16);
    let file = std::fs::File::open(path).expect("open corpus");
    for (index, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line.expect("read corpus");
        let Some(smiles) = line.split_whitespace().next() else {
            continue;
        };
        match chematic_smiles::parse(smiles) {
            Ok(mol) => {
                for seed in 1..=seeds {
                    println!(
                        "{index}\t{seed}\t{}",
                        chematic_smiles::random_smiles(&mol, seed)
                    );
                }
            }
            Err(error) => println!("{index}\t0\tERROR:{error}"),
        }
    }
}
