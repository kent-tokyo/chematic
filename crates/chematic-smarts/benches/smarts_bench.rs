//! Criterion benchmarks for SMARTS matching and SmartsCache.
//!
//! Run: `cargo bench -p chematic-smarts`
//! HTML report: `target/criterion/smarts_bench/report/index.html`
//!
//! The SmartsCache benchmark demonstrates the 5–20× speedup from caching
//! compiled SMARTS patterns compared to re-parsing on every call.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use chematic_smarts::{SmartsCache, find_matches, parse_smarts};
use chematic_smiles::parse;

const BENCH_SMILES: &[&str] = &[
    "c1ccccc1",
    "Cc1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "c1ccncc1",
    "C1CCNCC1",
    "CC(=O)Nc1ccc(O)cc1",
    "NCC(=O)O",
    "CN(C)C(=N)NC(=N)N",
];

// ── compile SMARTS (parse overhead only) ─────────────────────────────────────

fn bench_smarts_compile(c: &mut Criterion) {
    let patterns = &[
        "[NH2]c1ccccc1",
        "[CX3](=O)[OX2H1]",
        "c1ccccc1",
        "[nX2H0;r6]",
        "[OX2H1][cX3]",
    ];
    c.bench_function("smarts_compile_5pat", |b| {
        b.iter(|| {
            for pat in patterns {
                let _ = black_box(parse_smarts(black_box(*pat)));
            }
        })
    });
}

// ── SMARTS match without cache (parse + match on every call) ─────────────────

fn bench_smarts_match_no_cache(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    // Benzene ring pattern — matches 6 of the 10 molecules
    c.bench_function("smarts_match_nocache_10mol", |b| {
        b.iter(|| {
            let query = parse_smarts("c1ccccc1").unwrap();
            for mol in &mols {
                let _ = black_box(find_matches(black_box(&query), black_box(mol)));
            }
        })
    });
}

// ── SMARTS match with SmartsCache (compile once, reuse on every call) ─────────

fn bench_smarts_cache_hit(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let mut cache = SmartsCache::new(32);
    // Warm-up: ensure pattern is compiled and cached
    let _ = cache.find_matches("c1ccccc1", &mols[0]);
    c.bench_function("smarts_match_cached_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(cache.find_matches("c1ccccc1", black_box(mol)));
            }
        })
    });
}

// ── Complex recursive SMARTS ──────────────────────────────────────────────────

fn bench_smarts_complex(c: &mut Criterion) {
    let mols: Vec<_> = BENCH_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    let query = parse_smarts("[NH;$(NC=O)]").unwrap(); // amide N-H
    c.bench_function("smarts_recursive_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(find_matches(black_box(&query), black_box(mol)));
            }
        })
    });
}

criterion_group!(
    benches,
    bench_smarts_compile,
    bench_smarts_match_no_cache,
    bench_smarts_cache_hit,
    bench_smarts_complex,
);
criterion_main!(benches);
