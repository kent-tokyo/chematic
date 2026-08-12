//! Corpus-weighted `run_reactants` performance report.
//!
//! Built while diagnosing the RENKIN-reported `run_reactants`/`apply_retro`
//! performance regression between chematic 0.4.25 and 0.4.30 (see
//! `docs/rfcs/reaction_transform_perf.md`). A single-molecule/single-template
//! microbenchmark already demonstrably failed to catch that regression (it
//! is specific to *symmetric* intermediate molecules -- plain rings, cages,
//! `CF3`/`tBu`-style substituents -- which a handful of drug-like root
//! targets rarely exercise); this report instead runs a cross-product of
//! many template rules against many probe molecules, **and** feeds the
//! fragments produced by round 1 back in as round-2 probe molecules, to put
//! the same kind of symmetric intermediates a real multi-step retrosynthesis
//! search encounters into the measured population.
//!
//! # External corpora (optional)
//!
//! By default this runs against the small, hand-authored, MIT/BSD-clean
//! witness fixtures committed in `fixtures/` (`witness_templates.smirks`,
//! `witness_molecules.smi`). Point at a larger, external corpus instead via:
//!
//! - `RENKIN_TEMPLATES=/path/to/templates.smi` -- one SMIRKS per line
//!   (optionally `SMIRKS<TAB>count`, matching RENKIN's own
//!   `templates_extracted_5000.smi` format; the count column is ignored here).
//! - `RENKIN_PROBE=/path/to/probe.smi` -- one SMILES per line (optionally
//!   `SMILES<space>name`).
//!
//! Neither this repository nor this file ever contains or ships RENKIN's own
//! corpus data -- only the paths are read, from wherever the caller points.
//!
//! Other env vars: `RXN_PERF_MAX_TEMPLATES` (cap template count),
//! `RXN_PERF_DEPTH` (how many fragment-feedback rounds, default 2).
//!
//! Run for comparability with RENKIN's own methodology (`RAYON_NUM_THREADS=2`
//! pinned on both arms of a before/after comparison) via:
//! ```text
//! RAYON_NUM_THREADS=2 cargo run --release -p chematic-rxn \
//!     --features perf-instrumentation --example reaction_transform_perf_report
//! ```
//! ponytail: this harness is single-threaded (chematic-rxn itself spawns no
//! threads), so `RAYON_NUM_THREADS` has no effect on it directly -- it is
//! set purely so the *process environment* matches RENKIN's own measurement
//! conditions exactly, per the task's comparability requirement. Add rayon
//! parallelism here if a future user needs this harness to also reproduce
//! RENKIN's cross-core contention characteristics.

use std::env;
use std::fs;
use std::time::{Duration, Instant};

use chematic_chem::standardize::{StandardizeOptions, ZwitterionHandling, standardize};
use chematic_rxn::run_reactants;
use chematic_smiles::{canonical_smiles, parse};

const DEFAULT_TEMPLATES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/fixtures/witness_templates.smirks"
);
const DEFAULT_PROBE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/fixtures/witness_molecules.smi"
);

struct CallStats {
    durations: Vec<Duration>,
}

impl CallStats {
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
        self.durations.iter().copied().max().unwrap_or_default()
    }

    fn total(&self) -> Duration {
        self.durations.iter().sum()
    }
}

/// Mirror of RENKIN's `apply_retro` -> `split_fragments` pattern: apply one
/// SMIRKS template to one molecule, then split every product on '.',
/// re-parse, and standardize each fragment -- the exact call sequence that
/// (per the bisect in docs/rfcs/reaction_transform_perf.md) actually carries the
/// regression, not `run_reactants` in isolation.
///
/// Note: this does *not* reimplement RENKIN's own aromatic-atom-without-a-
/// ring-closure fragment filter (their `split_fragments`'s BFS-leakage
/// guard) -- `fragments_parsed_ok` below counts successfully-`parse()`d
/// fragments, not post-filter survivors, and is named accordingly rather
/// than promising a filter this harness doesn't apply.
#[allow(clippy::too_many_arguments)]
fn apply_retro_like(
    smirks: &str,
    mol: &chematic_core::Molecule,
    opts: &StandardizeOptions,
    canonical_smiles_calls: &mut u64,
    standardize_calls: &mut u64,
    fragments_before_filter: &mut u64,
    fragments_parsed_ok: &mut u64,
    errors_swallowed: &mut u64,
) -> (Duration, Vec<String>) {
    let t0 = Instant::now();
    let mut fragment_smiles = Vec::new();
    match run_reactants(smirks, &[mol]) {
        Ok(product_sets) => {
            for product_set in product_sets {
                for product_mol in product_set {
                    let smi = canonical_smiles(&product_mol);
                    *canonical_smiles_calls += 1;
                    for frag in smi.split('.') {
                        *fragments_before_filter += 1;
                        if let Ok(m) = parse(frag) {
                            let std_mol = standardize(&m, opts);
                            *standardize_calls += 1;
                            let std_smi = canonical_smiles(&std_mol);
                            *canonical_smiles_calls += 1;
                            *fragments_parsed_ok += 1;
                            fragment_smiles.push(std_smi);
                        }
                    }
                }
            }
        }
        Err(_) => {
            // Matches RENKIN's `.unwrap_or_default()` on run_reactants' Result.
            *errors_swallowed += 1;
        }
    }
    (t0.elapsed(), fragment_smiles)
}

fn load_lines(path: &str, take_first_field: bool) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read '{path}': {e}"))
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            if take_first_field {
                l.split(['\t', ' ']).next().unwrap_or(l).to_string()
            } else {
                l.to_string()
            }
        })
        .collect()
}

fn main() {
    let templates_path = env::var("RENKIN_TEMPLATES").unwrap_or_else(|_| DEFAULT_TEMPLATES.into());
    let probe_path = env::var("RENKIN_PROBE").unwrap_or_else(|_| DEFAULT_PROBE.into());
    let max_templates: usize = env::var("RXN_PERF_MAX_TEMPLATES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let depth: usize = env::var("RXN_PERF_DEPTH")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);

    let mut templates = load_lines(&templates_path, true);
    templates.truncate(max_templates);
    let probe_smiles = load_lines(&probe_path, true);

    eprintln!(
        "templates_path={templates_path} ({} rules) probe_path={probe_path} ({} molecules) depth={depth}",
        templates.len(),
        probe_smiles.len()
    );

    let opts = StandardizeOptions {
        canonical_tautomer: false,
        neutralize_charges: false,
        remove_explicit_h: true,
        largest_fragment_only: false,
        zwitterion_handling: ZwitterionHandling::Keep,
    };

    // Round 0: the probe molecules themselves.
    let mut current_smiles: Vec<String> = probe_smiles;

    let mut stats = CallStats::new();
    let mut canonical_smiles_calls: u64 = 0;
    let mut standardize_calls: u64 = 0;
    let mut fragments_before_filter: u64 = 0;
    let mut fragments_parsed_ok: u64 = 0;
    let mut errors_swallowed: u64 = 0;
    let mut product_molecules_seen: u64 = 0;

    chematic_rxn::perf_counters::reset();
    let wall_t0 = Instant::now();

    for level in 0..depth {
        let mut next_smiles: Vec<String> = Vec::new();
        let mols: Vec<_> = current_smiles
            .iter()
            .filter_map(|s| parse(s).ok())
            .collect();
        eprintln!(
            "round {level}: {} input molecules ({} parsed)",
            current_smiles.len(),
            mols.len()
        );

        for smirks in &templates {
            for mol in &mols {
                let (elapsed, fragments) = apply_retro_like(
                    smirks,
                    mol,
                    &opts,
                    &mut canonical_smiles_calls,
                    &mut standardize_calls,
                    &mut fragments_before_filter,
                    &mut fragments_parsed_ok,
                    &mut errors_swallowed,
                );
                stats.durations.push(elapsed);
                product_molecules_seen += fragments.len() as u64;
                next_smiles.extend(fragments);
            }
        }

        // Dedup + cap: a real search would prune via a visited-set and beam
        // width; this keeps the next round's population bounded without
        // biasing which structures survive (sort+dedup is order-independent).
        next_smiles.sort_unstable();
        next_smiles.dedup();
        next_smiles.truncate(200);
        current_smiles = next_smiles;
    }

    let wall_elapsed = wall_t0.elapsed();
    let counters = chematic_rxn::perf_counters::snapshot();

    println!("=== reaction_transform_perf_report ===");
    println!("total_calls={}", stats.durations.len());
    println!("wall_elapsed_ms={}", wall_elapsed.as_millis());
    println!("elapsed_total_ms={}", stats.total().as_millis());
    println!(
        "elapsed_per_call_ns={:.1}",
        stats.total().as_nanos() as f64 / stats.durations.len().max(1) as f64
    );
    println!("p50_ns={}", stats.percentile(0.50).as_nanos());
    println!("p95_ns={}", stats.percentile(0.95).as_nanos());
    println!("p99_ns={}", stats.percentile(0.99).as_nanos());
    println!("max_ns={}", stats.max().as_nanos());
    println!("product_molecules_seen={product_molecules_seen}");
    println!("canonical_smiles_calls={canonical_smiles_calls}");
    println!("standardize_calls={standardize_calls}");
    println!("fragments_before_filter={fragments_before_filter}");
    println!("fragments_parsed_ok={fragments_parsed_ok}");
    println!("errors_swallowed_by_unwrap_or_default={errors_swallowed}");
    println!(
        "successful_match_rate={:.4}",
        product_molecules_seen as f64 / stats.durations.len().max(1) as f64
    );
    println!("--- chematic-rxn internal work counters (perf-instrumentation feature) ---");
    println!("run_reactants_calls={}", counters.run_reactants_calls);
    println!("reaction_parse_calls={}", counters.reaction_parse_calls);
    println!(
        "reactant_query_match_calls={}",
        counters.reactant_query_match_calls
    );
    println!("vf2_match_count={}", counters.vf2_match_count);
    println!(
        "match_combination_count={}",
        counters.match_combination_count
    );
    println!("build_product_calls={}", counters.build_product_calls);
    println!(
        "product_sets_before_dedup={}",
        counters.product_sets_before_dedup
    );
    println!(
        "product_molecules_built={}",
        counters.product_molecules_built
    );
    println!(
        "atoms_copied_to_products={}",
        counters.atoms_copied_to_products
    );
    println!(
        "bonds_copied_to_products={}",
        counters.bonds_copied_to_products
    );
    if counters.run_reactants_calls == 0 {
        eprintln!(
            "note: all internal work counters are zero -- rebuild with \
             --features perf-instrumentation to populate them."
        );
    }
}
