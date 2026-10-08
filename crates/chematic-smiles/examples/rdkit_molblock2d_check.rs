//! Compare `rdkit_mol_block_2d` with `Chem.MolToMolBlock` after
//! `rdDepictor.Compute2DCoords` for a TSV of `smiles<TAB>escaped block`
//! lines (newlines written as `\n`, backslashes as `\\`).
//!
//! `cargo run --release -p chematic-smiles --example rdkit_molblock2d_check -- file.tsv [max_shown]`

use std::io::BufRead;

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut it = s.chars();
    while let Some(c) = it.next() {
        if c == '\\' {
            match it.next() {
                Some('n') => out.push('\n'),
                Some('\\') => out.push('\\'),
                Some(o) => {
                    out.push('\\');
                    out.push(o);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: rdkit_molblock2d_check FILE [max_shown]");
    let max_shown: usize = args.next().map_or(10, |s| s.parse().expect("number"));
    let file = std::fs::File::open(&path).expect("open");
    let (mut ok, mut bad, mut err, mut skipped) = (0usize, 0usize, 0usize, 0usize);
    let mut shown = 0;
    for line in std::io::BufReader::new(file).lines() {
        let line = line.expect("read");
        let Some((smi, block)) = line.split_once('\t') else {
            continue;
        };
        if block.starts_with('<') {
            skipped += 1;
            continue;
        }
        let expected = unescape(block);
        let got = chematic_smiles::parse(smi)
            .map_err(|e| e.to_string())
            .and_then(|m| chematic_smiles::rdkit_mol_block_2d(&m).map_err(|e| e.to_string()));
        match got {
            Ok(g) if g == expected => ok += 1,
            Ok(g) => {
                bad += 1;
                if shown < max_shown {
                    shown += 1;
                    println!("MISMATCH {smi}");
                    for (a, b) in g.lines().zip(expected.lines()) {
                        if a != b {
                            println!("  ours  {a}\n  rdkit {b}");
                        }
                    }
                    if g.lines().count() != expected.lines().count() {
                        println!(
                            "  line counts {} vs {}",
                            g.lines().count(),
                            expected.lines().count()
                        );
                    }
                }
            }
            Err(e) => {
                err += 1;
                if shown < max_shown {
                    shown += 1;
                    println!("ERROR {smi}: {e}");
                }
            }
        }
    }
    println!("match {ok} mismatch {bad} error {err} skipped {skipped}");
}
