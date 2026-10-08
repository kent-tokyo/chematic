//! Compare `rdkit_smiles` with RDKit's `MolToSmiles` option profiles from a
//! TSV of `input<TAB>bits:root<TAB>expected` lines, where `bits` is
//! isomeric, kekule, canonical, allBondsExplicit, allHsExplicit as 0/1 and
//! `root` is `rootedAtAtom` (-1: none).
//!
//! `cargo run --release -p chematic-smiles --example rdkit_smiles_profiles_check -- file.tsv [max_shown]`

use std::collections::BTreeMap;
use std::io::BufRead;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: FILE.tsv [max_shown]");
    let max_shown: usize = args.next().map_or(20, |s| s.parse().expect("number"));
    let file = std::fs::File::open(&path).expect("open tsv");
    let mut per_profile: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let (mut ok, mut bad) = (0usize, 0usize);
    let mut shown = 0;
    for line in std::io::BufReader::new(file).lines() {
        let line = line.expect("read");
        let parts: Vec<&str> = line.split('\t').collect();
        let [input, profile, expected] = parts[..] else {
            continue;
        };
        let (bits, root) = profile.split_once(':').expect("bits:root");
        let b: Vec<bool> = bits.bytes().map(|c| c == b'1').collect();
        let root: i64 = root.parse().expect("root");
        let params = chematic_smiles::RdkitSmilesParams {
            isomeric: b[0],
            kekule: b[1],
            canonical: b[2],
            all_bonds_explicit: b[3],
            all_hs_explicit: b[4],
            rooted_at_atom: (root >= 0).then_some(root as usize),
        };
        let got = chematic_smiles::parse(input)
            .map_err(|e| e.to_string())
            .and_then(|m| chematic_smiles::rdkit_smiles(&m, &params).map_err(|e| e.to_string()));
        let good = match (&got, expected) {
            (Ok(s), e) => s == e,
            (Err(_), "<ERROR>") => true,
            _ => false,
        };
        let key = format!(
            "{bits}:{}",
            if root < 0 {
                "none"
            } else if root == 0 {
                "first"
            } else {
                "last"
            }
        );
        let entry = per_profile.entry(key).or_default();
        if good {
            ok += 1;
            entry.0 += 1;
        } else {
            bad += 1;
            entry.1 += 1;
            if shown < max_shown {
                shown += 1;
                println!("MISMATCH [{profile}] {input}\n  rdkit {expected}\n  ours  {got:?}");
            }
        }
    }
    for (k, (g, b)) in &per_profile {
        if *b > 0 {
            println!("profile {k}: {g} ok, {b} bad");
        }
    }
    println!("profiles {} match {ok} mismatch {bad}", per_profile.len());
}
