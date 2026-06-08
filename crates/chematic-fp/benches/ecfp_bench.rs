use chematic_fp::{ecfp4, tanimoto_ecfp4};
use chematic_smiles::parse;
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

fn bench_ecfp4(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("ecfp4_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(ecfp4(black_box(mol)));
            }
        })
    });
}

fn bench_tanimoto(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("tanimoto_ecfp4_pairs", |b| {
        b.iter(|| {
            for i in 0..mols.len() {
                for j in i..mols.len() {
                    let _ = black_box(tanimoto_ecfp4(black_box(&mols[i]), black_box(&mols[j])));
                }
            }
        })
    });
}

criterion_group!(benches, bench_ecfp4, bench_tanimoto);
criterion_main!(benches);
