//! How often the RDKit-model hook runs and what it costs next to the
//! RDKit fingerprints: `cargo run --release --example rdkit_model_view_cost
//! -- corpus.smi [-v]`.

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("corpus path");
    let verbose = std::env::args().any(|a| a == "-v");
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
    let gated = mols
        .iter()
        .filter(|m| chematic_perception::rdkit_model_may_disagree(m))
        .count();
    let gate_us = t.elapsed().as_secs_f64() / n * 1e6;
    let mut corrected = 0;
    for (m, s) in mols.iter().zip(&smis) {
        if chematic_perception::rdkit_model_may_disagree(m) {
            if verbose {
                println!("gated\t{s}");
            }
            // corrected relative to the parity-only view
            if let Ok(v) = chematic_perception::aromaticity_rdkit_parity_only(m)
                && chematic_smiles::rdkit_model_correct_view(m, &v).is_some()
            {
                corrected += 1;
                if verbose {
                    println!("corrected\t{s}");
                }
            }
        }
    }
    println!(
        "rows {} gated {gated} corrected {corrected} gate {gate_us:.2} us",
        mols.len()
    );
    let mols = fresh();
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
        "parity view (hooked)  {:.1} us",
        t.elapsed().as_secs_f64() / n * 1e6
    );
    let mols = fresh();
    let t = Instant::now();
    for m in &mols {
        let _ = chematic_fp::rdkit_morgan_ecfp4_bitvec(m);
    }
    println!(
        "morgan ecfp4 (fresh)  {:.1} us",
        t.elapsed().as_secs_f64() / n * 1e6
    );
}
