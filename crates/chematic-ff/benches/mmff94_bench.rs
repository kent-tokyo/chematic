use chematic_ff::{Mmff94EnergyModel, minimize_mmff94_lbfgs, mmff94_total_energy};
use chematic_smiles::parse;
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const BENCH_SMILES: &[&str] = &[
    "CC",
    "CCC",
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
];

fn prepared_energy(c: &mut Criterion) {
    let molecules: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let models: Vec<_> = molecules
        .iter()
        .map(|m| Mmff94EnergyModel::new(m).unwrap())
        .collect();
    let coordinates: Vec<Vec<[f64; 3]>> = molecules
        .iter()
        .map(|m| {
            (0..m.atom_count())
                .map(|i| [i as f64 * 1.2, 0.0, 0.0])
                .collect()
        })
        .collect();
    c.bench_function("mmff94_prepared_energy_6mol", |b| {
        b.iter(|| {
            for (model, coords) in models.iter().zip(&coordinates) {
                black_box(model.energy(black_box(coords)));
            }
        })
    });
}

fn one_shot_energy(c: &mut Criterion) {
    let molecules: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let coordinates: Vec<Vec<[f64; 3]>> = molecules
        .iter()
        .map(|m| {
            (0..m.atom_count())
                .map(|i| [i as f64 * 1.2, 0.0, 0.0])
                .collect()
        })
        .collect();
    c.bench_function("mmff94_one_shot_energy_6mol", |b| {
        b.iter(|| {
            for (mol, coords) in molecules.iter().zip(&coordinates) {
                black_box(mmff94_total_energy(mol, black_box(coords)).unwrap());
            }
        })
    });
}

fn lbfgs_minimize(c: &mut Criterion) {
    let molecules: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let coordinates: Vec<Vec<[f64; 3]>> = molecules
        .iter()
        .map(|m| {
            (0..m.atom_count())
                .map(|i| [i as f64 * 1.2, 0.1, 0.0])
                .collect()
        })
        .collect();
    c.bench_function("mmff94_lbfgs_6mol_8iter", |b| {
        b.iter(|| {
            for (mol, initial) in molecules.iter().zip(&coordinates) {
                let mut coords = initial.clone();
                black_box(minimize_mmff94_lbfgs(mol, &mut coords, 8).unwrap());
            }
        })
    });
}

criterion_group!(benches, prepared_energy, one_shot_energy, lbfgs_minimize);
criterion_main!(benches);
