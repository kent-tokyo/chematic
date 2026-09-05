use chematic_fp::{
    FpType, Map4Config, PreparedFingerprintIndex, ecfp4, map4, nearest_neighbors, tanimoto_ecfp4,
    tanimoto_matrix, tanimoto_matrix_parallel,
};
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

fn bench_map4(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let config = Map4Config::default();
    c.bench_function("map4_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(map4(black_box(mol), black_box(&config)));
            }
        })
    });
}

fn bench_prepared_search(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let index = PreparedFingerprintIndex::new(&mols, FpType::Ecfp4);
    c.bench_function("prepared_search_10mol", |b| {
        b.iter(|| {
            for query in &mols {
                let _ = black_box(index.search(black_box(query), 5));
            }
        })
    });
    c.bench_function("rebuild_search_10mol", |b| {
        b.iter(|| {
            for query in &mols {
                let _ = black_box(nearest_neighbors(black_box(query), &mols, 5, FpType::Ecfp4));
            }
        })
    });
}

fn bench_tanimoto_matrix(c: &mut Criterion) {
    let mols: Vec<_> = (0..256)
        .map(|i| parse(BENCH_SMILES[i % BENCH_SMILES.len()]).unwrap())
        .collect();
    let fps: Vec<_> = mols.iter().map(ecfp4).collect();
    c.bench_function("tanimoto_matrix_serial_256x256", |b| {
        b.iter(|| black_box(tanimoto_matrix(black_box(&fps), black_box(&fps))))
    });
    c.bench_function("tanimoto_matrix_parallel_256x256", |b| {
        b.iter(|| black_box(tanimoto_matrix_parallel(black_box(&fps), black_box(&fps))))
    });
}

criterion_group!(
    benches,
    bench_ecfp4,
    bench_tanimoto,
    bench_map4,
    bench_prepared_search,
    bench_tanimoto_matrix
);
criterion_main!(benches);
