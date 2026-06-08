use chematic_smiles::{canonical_smiles, parse};
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const BENCH_SMILES: &[&str] = &[
    "c1ccccc1",
    "Cc1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",                  // aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",             // caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",             // ibuprofen
    "c1ccncc1",                               // pyridine
    "c1ccoc1",                                // furan
    "C1CCCCC1",                               // cyclohexane
    "CC(=O)Nc1ccc(O)cc1",                     // paracetamol
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", // glucose
];

fn bench_parse(c: &mut Criterion) {
    c.bench_function("parse_smiles_10mol", |b| {
        b.iter(|| {
            for s in BENCH_SMILES {
                let _ = black_box(parse(black_box(s)));
            }
        })
    });
}

fn bench_canonical(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("canonical_smiles_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(canonical_smiles(black_box(mol)));
            }
        })
    });
}

criterion_group!(benches, bench_parse, bench_canonical);
criterion_main!(benches);
