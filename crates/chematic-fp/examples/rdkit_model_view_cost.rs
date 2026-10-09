//! Cost of chematic-smiles' RDKit sanitization model next to the RDKit
//! fingerprints that consult it: `cargo run --release --example
//! rdkit_model_view_cost -- corpus.smi`.

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("corpus path");
    let smis: Vec<String> = std::fs::read_to_string(path)
        .expect("read")
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_string))
        .collect();
    let fresh = || -> Vec<chematic_core::Molecule> {
        smis.iter()
            .filter_map(|s| chematic_smiles::parse(s).ok())
            .collect()
    };
    let mols = fresh();
    let n = mols.len() as f64;
    let t = Instant::now();
    for m in &mols {
        let _ = chematic_smiles::rdkit_sanitized_model(m);
    }
    println!(
        "rdkit_sanitized_model {:.1} us",
        t.elapsed().as_secs_f64() / n * 1e6
    );
    let mols = fresh();
    let t = Instant::now();
    for m in &mols {
        let _ = chematic_perception::apply_aromaticity_rdkit_parity_shared(m);
    }
    println!(
        "parity view          {:.1} us",
        t.elapsed().as_secs_f64() / n * 1e6
    );
    let mols = fresh();
    let t = Instant::now();
    for m in &mols {
        let _ = chematic_fp::rdkit_morgan_ecfp4_bitvec(m);
    }
    println!(
        "morgan ecfp4 (fresh) {:.1} us",
        t.elapsed().as_secs_f64() / n * 1e6
    );
    let mols = fresh();
    let t = Instant::now();
    for m in &mols {
        let _ = chematic_fp::rdkit_atom_pair_fp(m);
    }
    println!(
        "atom pair (fresh)    {:.1} us",
        t.elapsed().as_secs_f64() / n * 1e6
    );
}
