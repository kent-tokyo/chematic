//! Criterion benchmarks for SSSR ring perception and aromatic ring counting.
//!
//! Run: `cargo bench -p chematic-perception`
//! HTML report: `target/criterion/sssr_bench/report/index.html`
//!
//! This crate had no Criterion coverage at all before this file (see
//! project history: a 47x perf regression in `find_sssr`'s candidate sort
//! landed with Horton's algorithm and went unnoticed by
//! `bench-pr-gate.yml` for several rounds, since SSSR was never in its
//! benchmark list). `find_sssr_polycyclic_5mol` specifically targets that
//! regression class: Horton's algorithm generates O(V*E) candidate cycles,
//! so its candidate-sort cost scales with ring *count*, not just atom
//! count — a corpus of mostly single-ring/small molecules would not have
//! caught the original regression.

use chematic_perception::{count_aromatic_rings, find_sssr};
use chematic_smiles::parse;
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

/// Same 10-molecule set as chematic-smiles/chematic-fp/chematic-chem benchmarks.
const BENCH_SMILES: &[&str] = &[
    "c1ccccc1",                   // benzene
    "Cc1ccccc1",                  // toluene
    "CC(=O)Oc1ccccc1C(=O)O",      // aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C", // caffeine (fused bicyclic)
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", // ibuprofen
    "c1ccncc1",                   // pyridine
    "C1CCNCC1",                   // piperidine
    "CC(=O)Nc1ccc(O)cc1",         // paracetamol
    "NCC(=O)O",                   // glycine (amphoteric)
    "CN(C)C(=N)NC(=N)N",          // metformin (multiple basic N)
];

/// Polycyclic/fused/cage molecules chosen specifically to stress Horton's
/// O(V*E) candidate-cycle generation and sort — the code path that
/// regressed 47x. Deliberately NOT the shared 10-molecule set above, which
/// is mostly single-ring and would not exercise this path meaningfully.
const POLYCYCLIC_SMILES: &[&str] = &[
    "c1cc2ccc3cccc4ccc(c1)c2c34", // pyrene (4 fused aromatic rings)
    "C1C2CC3CC1CC(C2)C3",         // adamantane (cage, O(V*E) candidate blowup)
    "C12C3C4C1C5C4C3C25",         // cubane (cage)
    "c1ccc2ccc3ccccc3c2c1",       // anthracene (3 linearly fused rings)
    "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C", // cholesterol (4 fused rings)
];

fn bench_find_sssr(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("find_sssr_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(find_sssr(black_box(mol)));
            }
        })
    });
}

fn bench_find_sssr_polycyclic(c: &mut Criterion) {
    let mols: Vec<_> = POLYCYCLIC_SMILES
        .iter()
        .map(|s| parse(s).unwrap())
        .collect();
    c.bench_function("find_sssr_polycyclic_5mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(find_sssr(black_box(mol)));
            }
        })
    });
}

fn bench_count_aromatic_rings(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("count_aromatic_rings_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(count_aromatic_rings(black_box(mol)));
            }
        })
    });
}

criterion_group!(
    benches,
    bench_find_sssr,
    bench_find_sssr_polycyclic,
    bench_count_aromatic_rings
);
criterion_main!(benches);
