use chematic_3d::{generate_coords_etkdg, minimize_mmff94};
use chematic_smiles::parse;
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const BENCH_SMILES: &[&str] = &[
    "CC",
    "CCC",
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
];

fn etkdg_generation(c: &mut Criterion) {
    let molecules: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("etkdg_generation_5mol", |b| {
        b.iter(|| {
            for mol in &molecules {
                black_box(generate_coords_etkdg(black_box(mol)));
            }
        })
    });
}

fn mmff94_minimization(c: &mut Criterion) {
    let molecules: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let starts: Vec<_> = molecules.iter().map(generate_coords_etkdg).collect();
    c.bench_function("mmff94_minimization_5mol", |b| {
        b.iter(|| {
            for (mol, start) in molecules.iter().zip(&starts) {
                black_box(minimize_mmff94(black_box(mol), black_box(start.clone())));
            }
        })
    });
}

criterion_group!(benches, etkdg_generation, mmff94_minimization);
criterion_main!(benches);
