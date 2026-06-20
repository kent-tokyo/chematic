//! Criterion benchmarks for molecular descriptor calculation.
//!
//! Run: `cargo bench -p chematic-chem`
//! HTML report: `target/criterion/descriptor_bench/report/index.html`

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use chematic_chem::{
    admet_profile, hba_count, hbd_count, logp_crippen, molecular_weight, pka_acid, pka_base,
    predict_pka, qed, tpsa,
};
use chematic_smiles::parse;

/// Same 10-molecule set as chematic-smiles and chematic-fp benchmarks.
const BENCH_SMILES: &[&str] = &[
    "c1ccccc1",                   // benzene
    "Cc1ccccc1",                  // toluene
    "CC(=O)Oc1ccccc1C(=O)O",      // aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C", // caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", // ibuprofen
    "c1ccncc1",                   // pyridine
    "C1CCNCC1",                   // piperidine
    "CC(=O)Nc1ccc(O)cc1",         // paracetamol
    "NCC(=O)O",                   // glycine (amphoteric)
    "CN(C)C(=N)NC(=N)N",          // metformin (multiple basic N)
];

// ── descriptor batch (5 descriptors × 10 molecules) ──────────────────────────

fn bench_descriptors(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("descriptors_5x10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(molecular_weight(black_box(mol)));
                let _ = black_box(logp_crippen(black_box(mol)));
                let _ = black_box(tpsa(black_box(mol)));
                let _ = black_box(hbd_count(black_box(mol)));
                let _ = black_box(hba_count(black_box(mol)));
            }
        })
    });
}

// ── QED (computationally heavier) ────────────────────────────────────────────

fn bench_qed(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("qed_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(qed(black_box(mol)));
            }
        })
    });
}

// ── pKa prediction (SMARTS-based) ────────────────────────────────────────────

fn bench_pka_predict(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("pka_predict_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(predict_pka(black_box(mol)));
            }
        })
    });
}

fn bench_pka_scalar(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("pka_acid_base_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(pka_acid(black_box(mol)));
                let _ = black_box(pka_base(black_box(mol)));
            }
        })
    });
}

// ── ADMET profile (full bundle) ───────────────────────────────────────────────

fn bench_admet_profile(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("admet_profile_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(admet_profile(black_box(mol)));
            }
        })
    });
}

criterion_group!(
    benches,
    bench_descriptors,
    bench_qed,
    bench_pka_predict,
    bench_pka_scalar,
    bench_admet_profile,
);
criterion_main!(benches);
