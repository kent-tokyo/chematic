//! Stage timings of the RDKit-model pipeline on a corpus (run with
//! `RDKIT_PROFILE_CORPUS=path cargo test --release -p chematic-smiles
//! rdkit_profile -- --ignored --nocapture`).

use std::time::{Duration, Instant};

use super::*;

#[test]
#[ignore]
fn rdkit_profile() {
    let Ok(path) = std::env::var("RDKIT_PROFILE_CORPUS") else {
        return;
    };
    let text = std::fs::read_to_string(path).unwrap();
    let mols: Vec<_> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter_map(|s| crate::parse(s).ok())
        .collect();
    let mut t = [Duration::ZERO; 12];
    let names = [
        "from_chematic",
        "remove_hs+sanitize",
        "  find_rings",
        "  kekulize",
        "legacy_stereo",
        "rank_mol_atoms",
        "mol_to_smiles",
        "mol_to_smarts",
        "murcko_decompose",
        "murcko perception+smiles",
        "rdkit_smiles total",
        "total",
    ];
    let only = std::env::var("RDKIT_PROFILE_ONLY").ok();
    let t_all = Instant::now();
    if let Some(only) = only {
        for mol in &mols {
            match only.as_str() {
                "smiles" => drop(rdkit_canonical_smiles(mol)),
                "smarts" => drop(rdkit_smarts(mol, true, None)),
                "murcko" => drop(rdkit_murcko_scaffold(mol)),
                _ => {}
            }
        }
        eprintln!("{only} {:9.1} ms", t_all.elapsed().as_secs_f64() * 1e3);
        return;
    }
    for mol in &mols {
        let t0 = Instant::now();
        let Ok(mut m) = parse::from_chematic(mol) else {
            continue;
        };
        t[0] += t0.elapsed();
        let t0 = Instant::now();
        if sanitize::remove_hs_and_sanitize(&mut m).is_err() {
            continue;
        }
        t[1] += t0.elapsed();
        let mut c = m.clone();
        let t0 = Instant::now();
        let _ = c.find_rings();
        t[2] += t0.elapsed();
        let t0 = Instant::now();
        let _ = kekulize::kekulize(&mut c);
        t[3] += t0.elapsed();
        let t0 = Instant::now();
        stereo::legacy_stereo_perception(&mut m, true, true);
        t[4] += t0.elapsed();
        let t0 = Instant::now();
        let _ = rank::rank_mol_atoms_with(&m, true);
        t[5] += t0.elapsed();
        let t0 = Instant::now();
        let _ = write::mol_to_smiles(&m, &RdkitSmilesParams::default());
        t[6] += t0.elapsed();
        let t0 = Instant::now();
        let _ = smarts_write::mol_to_smarts(&m, true, None, true);
        t[7] += t0.elapsed();
        let t0 = Instant::now();
        let s = murcko::murcko_decompose(&m);
        t[8] += t0.elapsed();
        if let Ok(mut s) = s {
            let t0 = Instant::now();
            stereo::legacy_stereo_perception(&mut s, true, false);
            let _ = write::mol_to_smiles(&s, &RdkitSmilesParams::default());
            t[9] += t0.elapsed();
        }
        let t0 = Instant::now();
        let _ = rdkit_canonical_smiles(mol);
        t[10] += t0.elapsed();
    }
    t[11] = t_all.elapsed();
    for (n, d) in names.iter().zip(t) {
        eprintln!("{n:28} {:9.1} ms", d.as_secs_f64() * 1e3);
    }
}
