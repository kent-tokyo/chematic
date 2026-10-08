//! Compare `rdkit_canonical_smiles(rdkit_mol_from_inchi(inchi))` with RDKit's
//! `MolToSmiles(MolFromInchi(inchi))` from a TSV `inchi<TAB>expected`
//! (`<NONE>`: RDKit returns None).
//!
//! `cargo run --release -p chematic-inchi --features native-inchi --example rdkit_inchi_read_check -- FILE.tsv`

use std::io::BufRead;

fn main() {
    let path = std::env::args().nth(1).expect("FILE.tsv");
    let max_shown: usize = std::env::args().nth(2).map_or(20, |s| s.parse().unwrap());
    let (mut ok, mut bad, mut shown) = (0, 0, 0);
    for line in std::io::BufReader::new(std::fs::File::open(path).unwrap()).lines() {
        let line = line.unwrap();
        let Some((inchi, expected)) = line.split_once('\t') else {
            continue;
        };
        let got = chematic_inchi::rdkit_mol_from_inchi(inchi)
            .map_err(|e| e.to_string())
            .and_then(|m| chematic_smiles::rdkit_canonical_smiles(&m).map_err(|e| e.to_string()));
        let good = match (&got, expected) {
            (Ok(s), e) => s == e,
            (Err(_), "<NONE>") => true,
            _ => false,
        };
        if good {
            ok += 1;
        } else {
            bad += 1;
            if shown < max_shown {
                shown += 1;
                println!("MISMATCH {inchi}\n  rdkit {expected}\n  ours  {got:?}");
            }
        }
    }
    println!("match {ok} mismatch {bad}");
}
