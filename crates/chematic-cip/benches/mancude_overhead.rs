//! Criterion benchmark for the perf cost of Milestone 3B-1b's live MANCUDE wiring:
//! `assign_cip_accurate_experimental` (candidate -- builds a Kekulé-form clone +
//! `MancudeContext` per molecule) vs `assign_cip_accurate_experimental_without_mancude`
//! (baseline -- the pre-Milestone-3B-1b digraph shape, no Kekulé/MANCUDE work at all).
//!
//! Run: `cargo bench -p chematic-cip`
//! HTML report: `target/criterion/mancude_overhead/report/index.html`
//!
//! Criterion's own statistical regression/no-change verdict (printed on every run once a
//! `--save-baseline`/prior run exists to compare against) is the acceptance mechanism for
//! this benchmark -- see `docs/rfcs/cip_accurate_rfc.md`'s Milestone 3B closeout entry for why
//! a bespoke p95 harness isn't needed here, and for the full-corpus (~5,000 molecule)
//! measurement this in-repo 10-molecule set can't by itself substitute for.
//!
//! Molecules are the first 10 unique SMILES from `validation/cip_label_corpus.jsonl`'s
//! `aromatic_mancude` bucket (fused polycyclic tetracycline-like alkaloids and simple
//! monosubstituted-phenol stereocenters) -- deliberately not the shared 10-molecule
//! benchmark set used by other crates (mostly single-ring, per `sssr_bench.rs`'s own
//! precedent), chosen specifically to exercise the Kekulé-respelling + MANCUDE-typing
//! path this benchmark exists to measure.

use chematic_cip::{
    CipBudget, assign_cip_accurate_experimental, assign_cip_accurate_experimental_without_mancude,
};
use chematic_smiles::parse;
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const MANCUDE_SMILES: &[&str] = &[
    "C=CCCC[C@H](c1ccc(O)cc1)[C@@](C)(CC)c1ccc(O)cc1",
    "CC(=O)C1=C(O)[C@@]2(O)C(O)C3C(=O)c4c(O)cccc4[C@@](C)(O)C3CC2[C@H](N(C)C)C1=O",
    "CC(=O)O[C@@]12C(O)=C(C(N)=O)C(=O)[C@@H](N(C)C)C1[C@@H](O)C1C(C(=O)c3c(O)cccc3[C@@]1(C)O)C2O",
    "CC(=O)Oc1cccc2c1C(=O)C1C(O)[C@]3(OC(C)=O)C(O)=C(C(N)=O)C(=O)[C@@H](N(C)C)C3CC1[C@]2(C)O",
    "CC(=O)SC[C@H]1c2cccc(O)c2C(=O)C2C1[C@H](O)C1[C@H](N(C)C)C(=O)C(C(N)=O)=C(O)[C@@]1(O)C2O",
    "CC1c2ccccc2Oc2cccc([C@@](N)(C(=O)O)[C@H]3C[C@@H]3C(=O)O)c21",
    "CCCCCCCC[C@H](c1ccc(O)cc1)[C@@H](CC)c1ccc(O)cc1",
    "CCCCSC(=O)N[C@@H](c1ccco1)[C@@H](O)C(=O)O[C@H]1CC2(O)[C@@H](OC(=O)c3ccccc3)[C@H]3[C@](C)(C(=O)[C@H](OC(C)=O)C(=C1C)C2(C)C)[C@@H](O)CC1OC[C@]13OC(C)=O",
    "CCCC[C@H](c1ccc(O)cc1)[C@@H](CC)c1ccc(O)cc1",
    "CCCC[C@H](c1ccc(O)cc1)[C@@](C)(CC)c1ccc(O)cc1",
];

fn bench_candidate(c: &mut Criterion) {
    let mols: Vec<_> = MANCUDE_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function("assign_cip_accurate_experimental_10mol", |b| {
        b.iter(|| {
            for mol in &mols {
                let _ = black_box(assign_cip_accurate_experimental(
                    black_box(mol),
                    CipBudget::default_budget(),
                ));
            }
        })
    });
}

fn bench_baseline(c: &mut Criterion) {
    let mols: Vec<_> = MANCUDE_SMILES.iter().map(|s| parse(s).unwrap()).collect();
    c.bench_function(
        "assign_cip_accurate_experimental_without_mancude_10mol",
        |b| {
            b.iter(|| {
                for mol in &mols {
                    let _ = black_box(assign_cip_accurate_experimental_without_mancude(
                        black_box(mol),
                        CipBudget::default_budget(),
                    ));
                }
            })
        },
    );
}

criterion_group!(benches, bench_baseline, bench_candidate);
criterion_main!(benches);
