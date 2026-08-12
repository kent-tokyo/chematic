//! Feature-gated internal work counters for `run_reactants`'s hot path.
//!
//! Enabled only under the `perf-instrumentation` Cargo feature so ordinary
//! builds pay zero cost (every call below compiles to nothing when the
//! feature is off). Counters are process-global `AtomicU64`s (relaxed
//! ordering -- these are cheap diagnostic counters, not a correctness
//! synchronization primitive), safe to bump from any thread, matching how
//! `run_reactants` is actually called by parallel (e.g. rayon) consumers.
//!
//! Added while diagnosing the `run_reactants`/`apply_retro` performance
//! regression between chematic 0.4.25 and 0.4.30 (see
//! `docs/rfcs/reaction_transform_perf.md`): distinguishing "doing more work" from
//! "the same work costing more" was the whole point -- call/product/match
//! *counts* turned out flat across versions, while per-call wall-clock time
//! rose sharply, which is exactly the signature these counters are built to
//! surface (paired with the benchmark's own wall-clock timing).

#[cfg(feature = "perf-instrumentation")]
use std::sync::atomic::{AtomicU64, Ordering};

#[cfg(feature = "perf-instrumentation")]
static RUN_REACTANTS_CALLS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static REACTION_PARSE_CALLS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static REACTANT_QUERY_MATCH_CALLS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static VF2_MATCH_COUNT: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static MATCH_COMBINATION_COUNT: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static BUILD_PRODUCT_CALLS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static PRODUCT_SETS_BEFORE_DEDUP: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static PRODUCT_MOLECULES_BUILT: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static ATOMS_COPIED_TO_PRODUCTS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "perf-instrumentation")]
static BONDS_COPIED_TO_PRODUCTS: AtomicU64 = AtomicU64::new(0);

/// Snapshot of every counter at the moment of the call.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PerfCounters {
    pub run_reactants_calls: u64,
    pub reaction_parse_calls: u64,
    pub reactant_query_match_calls: u64,
    pub vf2_match_count: u64,
    pub match_combination_count: u64,
    pub build_product_calls: u64,
    pub product_sets_before_dedup: u64,
    pub product_molecules_built: u64,
    pub atoms_copied_to_products: u64,
    pub bonds_copied_to_products: u64,
}

/// Read every counter without resetting them. Returns all-zero when the
/// `perf-instrumentation` feature is disabled.
#[cfg(feature = "perf-instrumentation")]
pub fn snapshot() -> PerfCounters {
    PerfCounters {
        run_reactants_calls: RUN_REACTANTS_CALLS.load(Ordering::Relaxed),
        reaction_parse_calls: REACTION_PARSE_CALLS.load(Ordering::Relaxed),
        reactant_query_match_calls: REACTANT_QUERY_MATCH_CALLS.load(Ordering::Relaxed),
        vf2_match_count: VF2_MATCH_COUNT.load(Ordering::Relaxed),
        match_combination_count: MATCH_COMBINATION_COUNT.load(Ordering::Relaxed),
        build_product_calls: BUILD_PRODUCT_CALLS.load(Ordering::Relaxed),
        product_sets_before_dedup: PRODUCT_SETS_BEFORE_DEDUP.load(Ordering::Relaxed),
        product_molecules_built: PRODUCT_MOLECULES_BUILT.load(Ordering::Relaxed),
        atoms_copied_to_products: ATOMS_COPIED_TO_PRODUCTS.load(Ordering::Relaxed),
        bonds_copied_to_products: BONDS_COPIED_TO_PRODUCTS.load(Ordering::Relaxed),
    }
}

#[cfg(not(feature = "perf-instrumentation"))]
pub fn snapshot() -> PerfCounters {
    PerfCounters::default()
}

/// Reset every counter to zero (e.g. between benchmark segments). A no-op
/// when the `perf-instrumentation` feature is disabled.
#[cfg(feature = "perf-instrumentation")]
pub fn reset() {
    RUN_REACTANTS_CALLS.store(0, Ordering::Relaxed);
    REACTION_PARSE_CALLS.store(0, Ordering::Relaxed);
    REACTANT_QUERY_MATCH_CALLS.store(0, Ordering::Relaxed);
    VF2_MATCH_COUNT.store(0, Ordering::Relaxed);
    MATCH_COMBINATION_COUNT.store(0, Ordering::Relaxed);
    BUILD_PRODUCT_CALLS.store(0, Ordering::Relaxed);
    PRODUCT_SETS_BEFORE_DEDUP.store(0, Ordering::Relaxed);
    PRODUCT_MOLECULES_BUILT.store(0, Ordering::Relaxed);
    ATOMS_COPIED_TO_PRODUCTS.store(0, Ordering::Relaxed);
    BONDS_COPIED_TO_PRODUCTS.store(0, Ordering::Relaxed);
}

#[cfg(not(feature = "perf-instrumentation"))]
pub fn reset() {}

// Internal increment helpers used from `transform.rs`. All are no-ops
// (inlined away entirely) when `perf-instrumentation` is off.

#[cfg(feature = "perf-instrumentation")]
pub(crate) fn record_run_reactants_call() {
    RUN_REACTANTS_CALLS.fetch_add(1, Ordering::Relaxed);
}
#[cfg(not(feature = "perf-instrumentation"))]
#[inline(always)]
pub(crate) fn record_run_reactants_call() {}

#[cfg(feature = "perf-instrumentation")]
pub(crate) fn record_reaction_parse_call() {
    REACTION_PARSE_CALLS.fetch_add(1, Ordering::Relaxed);
}
#[cfg(not(feature = "perf-instrumentation"))]
#[inline(always)]
pub(crate) fn record_reaction_parse_call() {}

#[cfg(feature = "perf-instrumentation")]
pub(crate) fn record_reactant_query_match_call(match_count: usize) {
    REACTANT_QUERY_MATCH_CALLS.fetch_add(1, Ordering::Relaxed);
    VF2_MATCH_COUNT.fetch_add(match_count as u64, Ordering::Relaxed);
}
#[cfg(not(feature = "perf-instrumentation"))]
#[inline(always)]
pub(crate) fn record_reactant_query_match_call(_match_count: usize) {}

#[cfg(feature = "perf-instrumentation")]
pub(crate) fn record_match_combination() {
    MATCH_COMBINATION_COUNT.fetch_add(1, Ordering::Relaxed);
}
#[cfg(not(feature = "perf-instrumentation"))]
#[inline(always)]
pub(crate) fn record_match_combination() {}

#[cfg(feature = "perf-instrumentation")]
pub(crate) fn record_build_product_call(atoms: usize, bonds: usize) {
    BUILD_PRODUCT_CALLS.fetch_add(1, Ordering::Relaxed);
    PRODUCT_MOLECULES_BUILT.fetch_add(1, Ordering::Relaxed);
    ATOMS_COPIED_TO_PRODUCTS.fetch_add(atoms as u64, Ordering::Relaxed);
    BONDS_COPIED_TO_PRODUCTS.fetch_add(bonds as u64, Ordering::Relaxed);
}
#[cfg(not(feature = "perf-instrumentation"))]
#[inline(always)]
pub(crate) fn record_build_product_call(_atoms: usize, _bonds: usize) {}

#[cfg(feature = "perf-instrumentation")]
pub(crate) fn record_product_set() {
    PRODUCT_SETS_BEFORE_DEDUP.fetch_add(1, Ordering::Relaxed);
}
#[cfg(not(feature = "perf-instrumentation"))]
#[inline(always)]
pub(crate) fn record_product_set() {}
