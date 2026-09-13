//! Reusable target-side reaction matching benchmark for issue #531.
//!
//! Run with:
//! `cargo bench -p chematic-rxn --bench reaction_match_context`
//!
//! The three groups intentionally separate one-time context construction from
//! repeated query work.  The ring-aware arm is the existing compatible
//! baseline; the context arm owns the target snapshot and reuses its rings.

use std::hint::black_box;

use chematic_perception::find_sssr;
use chematic_rxn::{PreparedReaction, ReactionMatchContext};
use chematic_smiles::parse;
use criterion::{Criterion, criterion_group, criterion_main};

const TARGET: &str = "CC(=O)Oc1ccccc1C(=O)O";
const TEMPLATES: &[&str] = &[
    "[C:1](=[O:2])[O:3][c:4]>>[C:1](=[O:2])O.[O:3][c:4]",
    "[c:1]1[c:2][c:3][c:4][c:5][c:6]1>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1",
    "[O:1]=[C:2][O:3]>>[O:1]=[C:2][O:3]",
    "[C:1][C:2]>>[C:1][C:2]",
    "[O:1][C:2]>>[O:1][C:2]",
    "[c:1][c:2]>>[c:1][c:2]",
    "[C:1](=[O:2])>>[C:1](=[O:2])",
    "[O:1]=[C:2][c:3]>>[O:1]=[C:2][c:3]",
    "[C:1][O:2][c:3]>>[C:1][O:2][c:3]",
    "[c:1]>>[c:1]",
];

fn prepared_templates() -> Vec<PreparedReaction> {
    TEMPLATES
        .iter()
        .map(|template| PreparedReaction::new(template).unwrap())
        .collect()
}

fn bench_context_build(c: &mut Criterion) {
    let target = parse(TARGET).unwrap();
    c.bench_function("reaction_match_context_build", |b| {
        b.iter(|| black_box(ReactionMatchContext::new(black_box(&target), None)))
    });
}

fn bench_prepared_ring_aware(c: &mut Criterion) {
    let target = parse(TARGET).unwrap();
    let rings = find_sssr(&target);
    let prepared = prepared_templates();
    c.bench_function("reaction_match_prepared_ring_aware_10_queries", |b| {
        b.iter(|| {
            for query in &prepared {
                black_box(
                    query
                        .find_matches_with_rings(&[&target], &[&rings])
                        .unwrap(),
                );
            }
        })
    });
}

fn bench_prepared_context(c: &mut Criterion) {
    let target = parse(TARGET).unwrap();
    let context = ReactionMatchContext::new(&target, None);
    let prepared = prepared_templates();
    c.bench_function("reaction_match_prepared_context_10_queries", |b| {
        b.iter(|| {
            for query in &prepared {
                black_box(
                    query
                        .find_matches_with_context(black_box(&context))
                        .unwrap(),
                );
            }
        })
    });
}

criterion_group!(
    benches,
    bench_context_build,
    bench_prepared_ring_aware,
    bench_prepared_context
);
criterion_main!(benches);
