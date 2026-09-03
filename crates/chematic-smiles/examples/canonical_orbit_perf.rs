//! Automorphism-orbit-pruning performance report for `canonical_smiles`.
//!
//! Built while implementing `fix/canonical-automorphism-pruning` (see
//! `docs/rfcs/canonical_automorphism_pruning.md` and `docs/rfcs/reaction_transform_perf.md`).
//! Compares the current orbit-pruned engine (`canonical_smiles`) against the
//! exact pre-pruning exhaustive search
//! (`chematic_smiles::canonical::legacy_canonical_smiles_for_benchmark`,
//! `#[doc(hidden)]`, same code as before this PR, kept only for this
//! comparison) across three tiers:
//!
//! - **Tier A (high symmetry)**: plain rings, cages, PAHs, `CF3`/`tBu`,
//!   multiple independent Boc/pivaloyl groups, repeated disconnected
//!   components. Performance gate: >= 80% leaf-count reduction on
//!   benzene/adamantane/coronene; >= 2x geometric-mean speedup overall; 0
//!   search-budget exhaustions on the multi-Boc/pivaloyl fixtures.
//! - **Tier B (low symmetry)**: chiral drug-like/heteroatom-rich molecules.
//!   Gate: <= 5% geometric-mean regression (most take exactly 1 branch under
//!   both engines, so this should be near-zero by construction).
//! - **Tier C (external corpus, optional)**: one SMILES per line via
//!   `CANONICAL_ORBIT_PERF_CORPUS=/path/to/file.smi`. Never required to be
//!   present in this repository -- point it at, e.g., `~/Downloads/SMILES.csv`
//!   (the project's 5,000-molecule benchmark corpus, see `scripts/bench5k.py`)
//!   or any other externally supplied file.
//!
//! Correctness differential (old string == new string) is a hard gate on
//! Tiers A and B unconditionally. On Tier C it additionally requires
//! `CANONICAL_ORBIT_PERF_RUN_LEGACY=1` (see below) -- **without it, the old
//! engine never runs on the corpus and the differential is SKIPPED, not
//! passed**; the report and exit status say so explicitly rather than
//! silently asserting 0 mismatches against zero comparisons (a real gap an
//! independent Round-3 performance-claim audit found in this file, PR #193).
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-smiles --example canonical_orbit_perf
//! cargo run --release -p chematic-smiles --features canonical-search-instrumentation \
//!     --example canonical_orbit_perf
//! # Timing/leaf-count only on Tier C, correctness differential SKIPPED there:
//! CANONICAL_ORBIT_PERF_CORPUS=~/Downloads/SMILES.csv \
//!     cargo run --release -p chematic-smiles --example canonical_orbit_perf
//! # Full old-vs-new correctness differential on Tier C too (slower: runs the
//! # exhaustive legacy engine on every corpus molecule):
//! CANONICAL_ORBIT_PERF_CORPUS=~/Downloads/SMILES.csv \
//! CANONICAL_ORBIT_PERF_RUN_LEGACY=1 \
//!     cargo run --release -p chematic-smiles --example canonical_orbit_perf
//! ```

use std::env;
use std::fs;
use std::time::{Duration, Instant};

use chematic_smiles::canonical::{
    legacy_branch_count_for_benchmark, legacy_canonical_smiles_for_benchmark,
};
use chematic_smiles::{
    CanonicalSearchStats, canonical_smiles, parse, reset_search_stats, search_stats_snapshot,
};

struct Stats {
    durations: Vec<Duration>,
}

impl Stats {
    fn new() -> Self {
        Self {
            durations: Vec::new(),
        }
    }

    fn percentile(&self, p: f64) -> Duration {
        if self.durations.is_empty() {
            return Duration::ZERO;
        }
        let mut sorted = self.durations.clone();
        sorted.sort_unstable();
        let idx = ((sorted.len() - 1) as f64 * p).round() as usize;
        sorted[idx]
    }

    fn max(&self) -> Duration {
        self.durations
            .iter()
            .copied()
            .max()
            .unwrap_or(Duration::ZERO)
    }

    fn total(&self) -> Duration {
        self.durations.iter().sum()
    }

    /// Geometric mean, in nanoseconds. Zero-duration entries are floored to
    /// 1ns so a free-running clock's occasional true-zero sample doesn't
    /// poison the product with a zero.
    fn geomean_ns(&self) -> f64 {
        if self.durations.is_empty() {
            return 0.0;
        }
        let log_sum: f64 = self
            .durations
            .iter()
            .map(|d| (d.as_nanos().max(1) as f64).ln())
            .sum();
        (log_sum / self.durations.len() as f64).exp()
    }
}

fn fmt_dur(d: Duration) -> String {
    if d.as_micros() < 1000 {
        format!("{}us", d.as_micros())
    } else {
        format!("{:.2}ms", d.as_secs_f64() * 1000.0)
    }
}

/// Tier A: high-symmetry fixtures. Cubane is built programmatically (cube
/// graph Q3, 8 atoms/12 bonds, |Aut|=48) rather than reusing this repo's
/// pre-existing "cubane" SMILES fixture used elsewhere
/// (`C12C3C4C1C1C2C3C41`, an 11-bond graph -- NOT a valid cube graph), per
/// this task's explicit instruction not to reuse that mislabeled fixture.
fn tier_a_high_symmetry() -> Vec<(&'static str, String)> {
    let mut v: Vec<(&'static str, String)> = vec![
        ("benzene", "c1ccccc1".to_string()),
        ("cyclohexane", "C1CCCCC1".to_string()),
        ("neopentane", "CC(C)(C)C".to_string()),
        ("tert-butanol", "CC(C)(C)O".to_string()),
        ("trifluoromethylbenzene", "FC(F)(F)c1ccccc1".to_string()),
        ("adamantane", "C1C2CC3CC1CC(C2)C3".to_string()),
        ("naphthalene", "c1ccc2ccccc2c1".to_string()),
        ("biphenyl", "c1ccc(-c2ccccc2)cc1".to_string()),
        (
            "multi-Boc intermediate",
            "CC(C)(C)OC(=O)NCC(CCNC(=O)OC(C)(C)C)CNC(=O)OC(C)(C)C".to_string(),
        ),
        (
            "multi-pivaloyl intermediate",
            "CC(C)(C)C(=O)NCC(CCNC(=O)C(C)(C)C)CNC(=O)C(C)(C)C".to_string(),
        ),
        (
            "repeated disconnected components",
            "CC(C)(C)C.CC(C)(C)C.CC(C)(C)C".to_string(),
        ),
    ];
    v.push(("cubane", cubane_smiles()));
    v.push(("coronene", coronene_smiles()));
    v
}

/// Build a correct cube graph (Q3): two 4-cycles `0-1-2-3-0` and
/// `4-5-6-7-4`, joined by the 4 "vertical" edges `0-4, 1-5, 2-6, 3-7`. Every
/// vertex has degree 3, 8 atoms, 12 bonds -- the real cubane skeleton
/// (`|Aut(Q3)| = 48`), unlike this repo's pre-existing, differently-shaped
/// "cubane" SMILES fixture used elsewhere (11 bonds).
fn cubane_smiles() -> String {
    use chematic_core::{BondOrder, Element, MoleculeBuilder};
    let mut b = MoleculeBuilder::new();
    let atoms: Vec<_> = (0..8)
        .map(|_| b.add_atom(chematic_core::Atom::organic(Element::C)))
        .collect();
    let edges: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    for (i, j) in edges {
        b.add_bond(atoms[i], atoms[j], BondOrder::Single).unwrap();
    }
    let mol = b.build();
    assert_eq!(mol.atom_count(), 8);
    assert_eq!(mol.bond_count(), 12);
    for i in 0..8 {
        assert_eq!(
            mol.degree(chematic_core::AtomIdx(i as u32)),
            3,
            "cube graph must be 3-regular"
        );
    }
    canonical_smiles(&mol)
}

/// Build a verified, correctly-sized coronene skeleton (C24H12, 24 atoms,
/// 30 bonds, 7 fused aromatic hexagons) geometrically: a flower of 7 unit
/// hexagons (1 center + 6 pointy-top neighbors) in axial hex coordinates,
/// corners deduped by rounded pixel coordinate, edges taken from each
/// hexagon's own 6-cycle (also deduped where two hexagons share an edge).
///
/// This project's own pre-existing `coronene` test fixture elsewhere in
/// this repo (`c1ccc2ccc3ccc4ccc5ccc6ccccc6c5c4c3c2c1`, used by the
/// previously-`#[ignore]`d `coronene_canonical_known_bug` test) was found
/// during this task to parse to **26** atoms, not the 24 a real coronene
/// (C24H12) has -- the same class of mislabeled-fixture problem this task
/// was warned about for "cubane" elsewhere in this repo. That pre-existing
/// fixture's test coverage is still valid (it is a real, if differently
/// sized, highly symmetric fused-ring molecule, and the idempotence fix it
/// exercises is real -- see `docs/rfcs/canonical_automorphism_pruning.md`), but
/// this function builds an independently *verified* correctly-sized
/// coronene from first principles (geometric construction, not a
/// hand-typed/memorized SMILES string) for this task's own performance
/// gate, which specifically names "coronene" as a required fixture.
fn coronene_smiles() -> String {
    use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, MoleculeBuilder};
    use std::collections::{HashMap, HashSet};

    let axial_centers: [(i64, i64); 7] =
        [(0, 0), (1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)];
    let size = 1.0_f64;
    let key = |x: f64, y: f64| -> (i64, i64) {
        ((x * 1000.0).round() as i64, (y * 1000.0).round() as i64)
    };

    let mut corner_of: HashMap<(i64, i64), AtomIdx> = HashMap::new();
    let mut builder = MoleculeBuilder::new();
    let mut seen_edges: HashSet<(i64, i64, i64, i64)> = HashSet::new();

    for &(q, r) in &axial_centers {
        let cx = size * (3f64.sqrt() * q as f64 + 3f64.sqrt() / 2.0 * r as f64);
        let cy = size * (1.5 * r as f64);
        let mut ring_atoms = Vec::with_capacity(6);
        for i in 0..6 {
            let angle = (std::f64::consts::PI / 180.0) * (60.0 * i as f64 - 30.0);
            let x = cx + size * angle.cos();
            let y = cy + size * angle.sin();
            let k = key(x, y);
            let atom = *corner_of
                .entry(k)
                .or_insert_with(|| builder.add_atom(Atom::organic(Element::C)));
            ring_atoms.push((k, atom));
        }
        for i in 0..6 {
            let (ka, a) = ring_atoms[i];
            let (kb, b) = ring_atoms[(i + 1) % 6];
            let ek = if ka <= kb {
                (ka.0, ka.1, kb.0, kb.1)
            } else {
                (kb.0, kb.1, ka.0, ka.1)
            };
            if seen_edges.insert(ek) {
                builder.add_bond(a, b, BondOrder::Single).unwrap();
            }
        }
    }
    let skeleton = builder.build();
    assert_eq!(
        skeleton.atom_count(),
        24,
        "coronene must have 24 carbons (C24H12)"
    );
    assert_eq!(
        skeleton.bond_count(),
        30,
        "coronene skeleton must have 30 bonds"
    );

    // Find a Kekule structure (perfect matching: every atom incident to
    // exactly one double bond) via trivial backtracking.
    fn find_matching(mol: &chematic_core::Molecule) -> Vec<BondIdx> {
        let n = mol.atom_count();
        let mut matched = vec![false; n];
        let mut chosen = Vec::new();
        fn go(
            mol: &chematic_core::Molecule,
            matched: &mut [bool],
            chosen: &mut Vec<BondIdx>,
        ) -> bool {
            let n = matched.len();
            let Some(u) = (0..n).find(|&i| !matched[i]) else {
                return true;
            };
            let candidates: Vec<(AtomIdx, BondIdx)> = mol
                .neighbors(AtomIdx(u as u32))
                .filter(|(nb, _)| !matched[nb.0 as usize])
                .collect();
            for (nb, bidx) in candidates {
                matched[u] = true;
                matched[nb.0 as usize] = true;
                chosen.push(bidx);
                if go(mol, matched, chosen) {
                    return true;
                }
                chosen.pop();
                matched[u] = false;
                matched[nb.0 as usize] = false;
            }
            false
        }
        assert!(
            go(mol, &mut matched, &mut chosen),
            "coronene skeleton has no Kekule structure"
        );
        chosen
    }
    let matching: HashSet<BondIdx> = find_matching(&skeleton).into_iter().collect();

    let mut kek_builder = MoleculeBuilder::new();
    for (_, atom) in skeleton.atoms() {
        kek_builder.add_atom(atom.clone());
    }
    for (bidx, bond) in skeleton.bonds() {
        let order = if matching.contains(&bidx) {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        kek_builder.add_bond(bond.atom1, bond.atom2, order).unwrap();
    }
    let kekule_mol = kek_builder.build();

    let aromatic_mol = chematic_perception::apply_aromaticity(&kekule_mol);
    let rings = chematic_perception::count_aromatic_rings(&aromatic_mol);
    assert_eq!(rings, 7, "coronene must have 7 aromatic rings");
    canonical_smiles(&aromatic_mol)
}

/// Tier B: low-symmetry negative-control fixtures (chiral drug-like,
/// heteroatom-rich). These should need exactly 1 individualize-refine
/// branch under both engines, so orbit pruning has ~nothing to prune -- the
/// gate here is that pruning must not *regress* this common case.
fn tier_b_low_symmetry() -> Vec<(&'static str, String)> {
    vec![
        (
            "cholesterol-like steroid",
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
                .to_string(),
        ),
        (
            "morphine-like polycyclic alkaloid",
            "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5".to_string(),
        ),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O".to_string()),
        ("L-alanine", "N[C@@H](C)C(=O)O".to_string()),
        (
            "heteroatom-rich asymmetric",
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1"
                .to_string(),
        ),
        ("caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C".to_string()),
    ]
}

fn tier_c_external_corpus() -> Option<Vec<(String, String)>> {
    let path = env::var("CANONICAL_ORBIT_PERF_CORPUS").ok()?;
    let contents =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    Some(
        contents
            .lines()
            .filter(|l| !l.trim().is_empty())
            .enumerate()
            .map(|(i, l)| {
                (
                    format!("corpus#{i}"),
                    l.split_whitespace().next().unwrap_or(l).to_string(),
                )
            })
            .collect(),
    )
}

struct TierReport {
    name: &'static str,
    old_stats: Stats,
    new_stats: Stats,
    mismatches: Vec<String>,
    parse_failures: usize,
    empty_outputs: usize,
    n: usize,
    /// Whether the old-vs-new correctness differential actually ran for this
    /// tier (`run_legacy` as passed to `run_tier`). When `false`,
    /// `mismatches` is trivially empty because it was never populated --
    /// callers must not read `mismatches.is_empty()` as "0 verified
    /// mismatches" in that case. See the module doc comment.
    differential_ran: bool,
    /// Cumulative orbit-search instrumentation across every fixture in this
    /// tier (all-zero when `canonical-search-instrumentation` is disabled --
    /// see `canonical_search::stats`).
    search: CanonicalSearchStats,
    /// Sum of `legacy_branch_count_for_benchmark` across every fixture --
    /// the old engine's exhaustive leaf count, independent of the
    /// instrumentation feature (always available).
    old_branch_count: usize,
}

fn add_stats(a: &CanonicalSearchStats, b: &CanonicalSearchStats) -> CanonicalSearchStats {
    CanonicalSearchStats {
        nodes_visited: a.nodes_visited + b.nodes_visited,
        leaves_written: a.leaves_written + b.leaves_written,
        orbit_tests: a.orbit_tests + b.orbit_tests,
        orbit_unions: a.orbit_unions + b.orbit_unions,
        children_pruned: a.children_pruned + b.children_pruned,
        max_depth: a.max_depth.max(b.max_depth),
        largest_target_cell: a.largest_target_cell.max(b.largest_target_cell),
        budget_exhaustions: a.budget_exhaustions + b.budget_exhaustions,
    }
}

fn run_tier(name: &'static str, fixtures: &[(String, String)], run_legacy: bool) -> TierReport {
    let mut old_stats = Stats::new();
    let mut new_stats = Stats::new();
    let mut mismatches = Vec::new();
    let mut parse_failures = 0;
    let mut empty_outputs = 0;
    let mut search = CanonicalSearchStats::default();
    let mut old_branch_count = 0usize;

    for (i, (label, smi)) in fixtures.iter().enumerate() {
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(_) => {
                parse_failures += 1;
                continue;
            }
        };

        // Alternate which engine's timer starts first, by fixture index.
        // An independent Round-3 performance-claim audit found (via a
        // throwaway microbench, since deleted) that whichever arm runs
        // first in this loop measures ~10-15% slower regardless of which
        // engine it is -- a harness artifact, not a real per-call effect.
        // Always timing new-then-old (the previous behavior) biased every
        // per-call comparison in the same direction; alternating cancels it
        // out in aggregate instead of just disclosing it.
        let new_first = i % 2 == 0;

        let (new_out, old_out_opt);
        if new_first {
            reset_search_stats();
            let t0 = Instant::now();
            let out = canonical_smiles(&mol);
            new_stats.durations.push(t0.elapsed());
            let per_fixture = search_stats_snapshot();
            search = add_stats(&search, &per_fixture);
            new_out = out;

            old_out_opt = if run_legacy {
                let t1 = Instant::now();
                let out = legacy_canonical_smiles_for_benchmark(&mol);
                old_stats.durations.push(t1.elapsed());
                old_branch_count += legacy_branch_count_for_benchmark(&mol);
                Some(out)
            } else {
                None
            };

            if env::var("CANONICAL_ORBIT_PERF_VERBOSE").is_ok() {
                let old_leaves = if run_legacy {
                    legacy_branch_count_for_benchmark(&mol)
                } else {
                    0
                };
                println!(
                    "    [{label}] new_us={} old_leaves={old_leaves} new_leaves={} nodes={} orbit_tests={} \
                     children_pruned={}",
                    new_stats.durations.last().unwrap().as_micros(),
                    per_fixture.leaves_written,
                    per_fixture.nodes_visited,
                    per_fixture.orbit_tests,
                    per_fixture.children_pruned
                );
            }
        } else {
            old_out_opt = if run_legacy {
                let t1 = Instant::now();
                let out = legacy_canonical_smiles_for_benchmark(&mol);
                old_stats.durations.push(t1.elapsed());
                old_branch_count += legacy_branch_count_for_benchmark(&mol);
                Some(out)
            } else {
                None
            };

            reset_search_stats();
            let t0 = Instant::now();
            let out = canonical_smiles(&mol);
            new_stats.durations.push(t0.elapsed());
            let per_fixture = search_stats_snapshot();
            search = add_stats(&search, &per_fixture);
            new_out = out;

            if env::var("CANONICAL_ORBIT_PERF_VERBOSE").is_ok() {
                let old_leaves = if run_legacy {
                    legacy_branch_count_for_benchmark(&mol)
                } else {
                    0
                };
                println!(
                    "    [{label}] new_us={} old_leaves={old_leaves} new_leaves={} nodes={} orbit_tests={} \
                     children_pruned={}",
                    new_stats.durations.last().unwrap().as_micros(),
                    per_fixture.leaves_written,
                    per_fixture.nodes_visited,
                    per_fixture.orbit_tests,
                    per_fixture.children_pruned
                );
            }
        }

        if new_out.is_empty() {
            empty_outputs += 1;
        }

        if let Some(old_out) = old_out_opt
            && old_out != new_out
        {
            mismatches.push(format!("{label} ({smi}): old='{old_out}' new='{new_out}'"));
        }
    }

    TierReport {
        name,
        old_stats,
        search,
        old_branch_count,
        new_stats,
        mismatches,
        parse_failures,
        empty_outputs,
        n: fixtures.len(),
        differential_ran: run_legacy,
    }
}

fn print_report(r: &TierReport) {
    println!("\n=== {} ({} fixtures) ===", r.name, r.n);
    println!("  parse failures: {}", r.parse_failures);
    println!("  empty outputs:  {}", r.empty_outputs);
    if r.differential_ran {
        println!(
            "  correctness mismatches (old vs new): {}",
            r.mismatches.len()
        );
        for m in &r.mismatches {
            println!("    MISMATCH: {m}");
        }
    } else {
        println!(
            "  correctness differential: SKIPPED (old engine not run for this tier -- \
             set CANONICAL_ORBIT_PERF_RUN_LEGACY=1 to actually check old-vs-new agreement)"
        );
    }
    if !r.old_stats.durations.is_empty() {
        println!(
            "  old: p50={} p95={} max={} total={} geomean={:.0}ns",
            fmt_dur(r.old_stats.percentile(0.50)),
            fmt_dur(r.old_stats.percentile(0.95)),
            fmt_dur(r.old_stats.max()),
            fmt_dur(r.old_stats.total()),
            r.old_stats.geomean_ns()
        );
    }
    println!(
        "  new: p50={} p95={} max={} total={} geomean={:.0}ns",
        fmt_dur(r.new_stats.percentile(0.50)),
        fmt_dur(r.new_stats.percentile(0.95)),
        fmt_dur(r.new_stats.max()),
        fmt_dur(r.new_stats.total()),
        r.new_stats.geomean_ns()
    );
    if !r.old_stats.durations.is_empty() {
        let speedup = r.old_stats.geomean_ns() / r.new_stats.geomean_ns().max(1.0);
        println!("  geometric-mean speedup (old/new): {speedup:.2}x");
    }
    if r.old_branch_count > 0 {
        println!(
            "  leaf count: old (exhaustive) = {}, new (orbit-pruned, requires \
             --features canonical-search-instrumentation for a nonzero value) = {}",
            r.old_branch_count, r.search.leaves_written
        );
        if r.search.leaves_written > 0 {
            let reduction = 100.0 * (r.old_branch_count as f64 - r.search.leaves_written as f64)
                / r.old_branch_count as f64;
            println!("    leaf-count reduction (old -> new): {reduction:.1}%");
        }
    }
    println!(
        "  orbit search (new engine, cumulative; all-zero unless built with \
         --features canonical-search-instrumentation):"
    );
    println!(
        "    nodes_visited={} leaves_written={} orbit_tests={} orbit_unions={} \
         children_pruned={} max_depth={} largest_target_cell={} budget_exhaustions={}",
        r.search.nodes_visited,
        r.search.leaves_written,
        r.search.orbit_tests,
        r.search.orbit_unions,
        r.search.children_pruned,
        r.search.max_depth,
        r.search.largest_target_cell,
        r.search.budget_exhaustions
    );
}

fn main() {
    println!("Canonical SMILES automorphism-orbit-pruning performance report");
    println!("(see docs/rfcs/canonical_automorphism_pruning.md)");

    let tier_a: Vec<(String, String)> = tier_a_high_symmetry()
        .into_iter()
        .map(|(n, s)| (n.to_string(), s))
        .collect();
    let tier_b: Vec<(String, String)> = tier_b_low_symmetry()
        .into_iter()
        .map(|(n, s)| (n.to_string(), s))
        .collect();

    let report_a = run_tier("Tier A: high symmetry", &tier_a, true);
    print_report(&report_a);

    let report_b = run_tier("Tier B: low symmetry (negative control)", &tier_b, true);
    print_report(&report_b);

    // `None` = no corpus configured; `Some(true)` = Tier C's differential
    // actually ran; `Some(false)` = corpus loaded but the differential was
    // skipped (see below) -- drives the final summary line's wording.
    let mut tier_c_differential_ran: Option<bool> = None;

    if let Some(corpus) = tier_c_external_corpus() {
        println!(
            "\nTier C external corpus: {} molecules loaded",
            corpus.len()
        );
        // Old engine is intentionally SKIPPED on a large external corpus by
        // default (it is exhaustive/unbounded-cap and can be extremely slow
        // on genuinely large symmetric molecules); set
        // CANONICAL_ORBIT_PERF_RUN_LEGACY=1 to force the full old-vs-new
        // differential anyway.
        let run_legacy = env::var("CANONICAL_ORBIT_PERF_RUN_LEGACY").is_ok();
        tier_c_differential_ran = Some(run_legacy);
        let report_c = run_tier("Tier C: external corpus", &corpus, run_legacy);
        print_report(&report_c);
        if run_legacy {
            assert_eq!(
                report_c.mismatches.len(),
                0,
                "external corpus produced {} old-vs-new canonical string differences",
                report_c.mismatches.len()
            );
        } else {
            println!(
                "\nTier C correctness differential: SKIPPED (CANONICAL_ORBIT_PERF_RUN_LEGACY not \
                 set) -- this run is leaf-count/timing-only, NOT a verified 0-mismatch claim. \
                 Re-run with CANONICAL_ORBIT_PERF_RUN_LEGACY=1 to actually check."
            );
        }
        assert_eq!(
            report_c.empty_outputs, 0,
            "external corpus produced empty canonical output(s)"
        );
    } else {
        println!(
            "\nTier C external corpus: not configured (set CANONICAL_ORBIT_PERF_CORPUS=/path/to/file.smi)"
        );
    }

    assert_eq!(
        report_a.mismatches.len(),
        0,
        "Tier A produced old-vs-new mismatches"
    );
    assert_eq!(
        report_b.mismatches.len(),
        0,
        "Tier B produced old-vs-new mismatches"
    );
    assert_eq!(
        report_a.empty_outputs, 0,
        "Tier A produced empty canonical output(s)"
    );
    assert_eq!(
        report_b.empty_outputs, 0,
        "Tier B produced empty canonical output(s)"
    );

    match tier_c_differential_ran {
        Some(false) => println!(
            "\nAll correctness gates passed for Tier A/B. Tier C's differential was SKIPPED \
             (see above) -- not a verified 0-mismatch claim for Tier C."
        ),
        _ => println!("\nAll correctness gates passed."),
    }
}
