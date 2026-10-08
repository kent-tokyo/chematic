//! Compare `rdkit_canonical_smiles` with RDKit's output for a TSV of
//! `input<TAB>expected` lines (expected = `Chem.MolToSmiles(Chem.MolFromSmiles(input))`).
//!
//! `cargo run --release -p chematic-smiles --example rdkit_smiles_check -- file.tsv [max_shown]`

use std::io::BufRead;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: rdkit_smiles_check FILE.tsv [max_shown]");
    let max_shown: usize = args.next().map_or(20, |s| s.parse().expect("number"));
    let file = std::fs::File::open(&path).expect("open tsv");
    let (mut ok, mut bad, mut err) = (0usize, 0usize, 0usize);
    let mut shown = 0;
    for line in std::io::BufReader::new(file).lines() {
        let line = line.expect("read");
        let mut parts = line.split('\t');
        let (Some(input), Some(expected)) = (parts.next(), parts.next()) else {
            continue;
        };
        let got = chematic_smiles::parse(input)
            .map_err(|e| e.to_string())
            .and_then(|m| chematic_smiles::rdkit_canonical_smiles(&m).map_err(|e| e.to_string()));
        if expected == "<NONE>" {
            // RDKit rejects the input: we must not return a string.
            match got {
                Err(_) => ok += 1,
                Ok(s) => {
                    bad += 1;
                    if shown < max_shown {
                        shown += 1;
                        println!("RDKIT-REJECTS {input}\n  ours  {s}");
                    }
                }
            }
            continue;
        }
        match got {
            Ok(s) if s == expected => ok += 1,
            Ok(s) => {
                bad += 1;
                if shown < max_shown {
                    shown += 1;
                    println!("MISMATCH {input}\n  rdkit {expected}\n  ours  {s}");
                }
            }
            Err(e) => {
                err += 1;
                if shown < max_shown {
                    shown += 1;
                    println!("ERROR {input}\n  rdkit {expected}\n  error {e}");
                }
            }
        }
    }
    println!("match {ok} mismatch {bad} error {err}");
}
