//! Conformer ensemble generation on top of [`crate::pipeline_v2::embed_pipeline_v2`]
//! — Track A, item A2 of `docs/rfcs/openeye_materials_advantage_rfc.md`.
//!
//! [`crate::generate_conformer_ensemble_with_config`] is the existing public
//! multi-conformer API, but it routes through the legacy
//! [`crate::etkdg::generate_coords_etkdg_with_noise`] path: no seed parameter
//! (backed by a process-global counter, not reproducible), no energy ranking,
//! and — per the fresh 265-molecule corpus re-audit (A1, PR #370) — a live
//! soundness defect (MMFF94 silently contributing zero energy/gradient for
//! atom-type pairs its parameter tables don't cover). `embed_pipeline_v2` has
//! none of these problems for a single conformer, but has no multi-conformer
//! API of its own.
//!
//! This module is a **new outer loop**, not an in-place extension of either
//! path: `embed_pipeline_v2`'s own `EmbedParameters::max_attempts`/
//! `timeout_ms` already retry the *same* target conformer on a bad stochastic
//! draw — they do not generate genuinely different conformers. [`embed_ensemble_v2`]
//! instead calls `embed_pipeline_v2` `count` times, once per deterministically
//! derived seed, and accumulates the results into the existing
//! [`crate::conformer::ConformerEnsemble`] storage type — no new
//! ensemble/storage type is introduced here.
//!
//! Leaves the legacy `conformer_ensemble()`/
//! `generate_conformer_ensemble_with_config` API completely untouched;
//! rewiring or deprecating it is a separate, future decision.
//!
//! # Two phases: generate, then select — never prune while generating
//!
//! An earlier version of this module pruned near-duplicates as each attempt
//! arrived, in attempt order, then only sorted the *survivors* by energy at
//! the end. That is wrong: if attempt 0 (high energy) and attempt 1 (lower
//! energy, a near-duplicate of attempt 0) both succeed, attempt 0 gets kept
//! first (nothing to compare it against yet) and attempt 1 is discarded as a
//! duplicate of it — the final energy sort can only reorder what already
//! survived, it cannot resurrect a better candidate that was discarded before
//! it ever got a chance to compete. The result was "first-generated
//! representative per RMSD cluster," not "lowest-energy representative per
//! RMSD cluster."
//!
//! [`embed_ensemble_v2`] therefore runs in two strict phases:
//! 1. **Generate**: run every attempt (budget permitting), recording every
//!    outcome — success or typed failure — with no pruning at all.
//! 2. **Select**: among the successes, group by `actual_force_field_used`
//!    (see below for why grouping matters), sort each group by ascending
//!    energy (attempts with no usable energy sort after ones that have it,
//!    tie-broken by attempt index), then greedily walk that order, keeping a
//!    candidate unless it is a near-duplicate of an *already-kept* candidate
//!    **in the same group**. This always selects the lowest-energy
//!    representative of every RMSD cluster, independent of generation order.
//!
//! # Energy comparability and pruning across `ForceFieldPolicy::Mmff94WithUffFallback`
//!
//! `PolicyMinimizeResult::energy_after` is on a physically meaningful scale
//! only *within* a single force field's own parameterization — MMFF94 and UFF
//! energies are not on a common reference zero. Under
//! `ForceFieldPolicy::Mmff94WithUffFallback`, individual attempts in the same
//! ensemble can resolve to different `actual_force_field_used` values (some
//! stay on MMFF94, some fall back to UFF) — see that field's own doc comment
//! in `minimize.rs` for why `fallback_reason.is_some()`, not
//! `actual_force_field_used != requested_force_field`, is the correct
//! "did this one fall back" check.
//!
//! This module never compares or prunes across that boundary: kept
//! conformers are grouped by `actual_force_field_used` (first-seen order),
//! each group is independently energy-ranked and independently
//! deduplicated, and [`EnsembleV2Result::mixed_force_field`] discloses when
//! more than one group has a kept member — a caller must not read
//! cross-group adjacency in `ensemble` as an energy comparison. Pruning
//! never crosses a group boundary either: a UFF-fallback candidate and an
//! MMFF94 candidate that happen to be near-duplicates in geometry are BOTH
//! kept, not resolved by an implicit "whichever came first" or "whichever
//! force field is generally more trustworthy" rule — this module makes no
//! such trustworthiness ranking. A future caller that wants one can filter
//! `ensemble`/`attempts` by `actual_force_field_used` itself.
//!
//! Energy is available at all only when the underlying [`EnergyReport`] is
//! not [`EnergyReport::None`] AND is finite — a non-finite (`NaN`/`±inf`)
//! energy is treated the same as "no usable energy," never fed into the
//! ascending-energy sort (which would otherwise place it in an
//! easy-to-misread pathological position). There is no separate "rank by
//! energy" flag to misconfigure: ranking happens whenever energy is actually
//! present and finite, keyed off the same `EnergyReport` the caller already
//! gets back, never off the *requested* policy.

use crate::clock::Instant;
use crate::conformer::ConformerEnsemble;
use crate::coords::Coords3D;
use crate::minimize::{EnergyReport, ForceFieldBridgeError, ForceFieldPolicy};
use crate::pipeline_v2::{PipelineV2Config, PipelineV2Failure, embed_pipeline_v2};
use chematic_core::Molecule;

use crate::distance_geometry_v2::derive_attempt_seed;

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Configuration for [`embed_ensemble_v2`].
#[derive(Debug, Clone)]
pub struct EnsembleV2Config {
    /// Per-conformer pipeline configuration, reused verbatim for every
    /// attempt except `embed.random_seed`, which this module overrides with
    /// a value derived from `base_seed` (see [`derive_attempt_seed`]).
    pub per_conformer: PipelineV2Config,
    /// Number of independent embedding attempts (before RMSD pruning). The
    /// kept ensemble may have fewer conformers than this, both from typed
    /// failures and from near-duplicate pruning. Also recorded verbatim on
    /// the result as [`EnsembleV2Result::requested_count`].
    pub count: usize,
    /// Base seed; attempt `i`'s seed is `derive_attempt_seed(base_seed, i)`
    /// — the same derivation `embed_pipeline_v2`'s own internal retry loop
    /// uses, reused here rather than inventing a second scheme. The same
    /// `base_seed` always reproduces the same ensemble.
    pub base_seed: u64,
    /// Minimum RMSD (Å) between kept conformers *within the same
    /// `actual_force_field_used` group* (see the module docs). Must be
    /// `0.0` (disables pruning, matching `ConformerConfig::rmsd_threshold`'s
    /// existing convention) or a positive, finite value — `embed_ensemble_v2`
    /// rejects `NaN`, infinite, or negative values via
    /// [`EnsembleV2ConfigError::InvalidRmsdThreshold`].
    pub rmsd_threshold: f64,
    /// When `true` (default), duplicate-checking uses
    /// [`crate::conformer::rmsd_symmetric`] (automorphism-aware — correct on
    /// molecules with interchangeable substituents, e.g. `-CF3`, at the cost
    /// of enumerating automorphisms per comparison). When `false`, uses the
    /// cheaper plain Kabsch RMSD ([`ConformerEnsemble::is_duplicate`]),
    /// which can treat truly-identical conformers as distinct on such
    /// molecules.
    pub use_symmetric_rmsd_pruning: bool,
    /// Wall-clock budget (milliseconds) across all `count` attempts
    /// combined. `None` = no limit. Distinct from
    /// `per_conformer.total_timeout_ms`, which budgets a *single* attempt.
    /// Checked between attempts only, matching this crate's established
    /// convention (see `PipelineV2Config::total_timeout_ms`'s own doc) —
    /// including failing closed on exactly `Some(0)` regardless of
    /// millisecond-resolution elapsed-time rounding. When the budget is
    /// exhausted before `count` attempts run,
    /// [`EnsembleV2Result::termination`] is
    /// [`EnsembleTermination::TimedOut`], never silently inferred from
    /// `attempts.len()`.
    pub ensemble_timeout_ms: Option<u64>,
}

impl EnsembleV2Config {
    /// Convenience constructor: `rmsd_threshold: 0.5` (matches
    /// `ConformerConfig`'s existing default), symmetric pruning on, no
    /// ensemble-level timeout. `per_conformer` has no `Default` of its own
    /// (see `PipelineV2Config`'s doc) so it is always an explicit argument.
    pub fn new(per_conformer: PipelineV2Config, count: usize, base_seed: u64) -> Self {
        Self {
            per_conformer,
            count,
            base_seed,
            rmsd_threshold: 0.5,
            use_symmetric_rmsd_pruning: true,
            ensemble_timeout_ms: None,
        }
    }
}

/// Rejected [`EnsembleV2Config`]. Mirrors this crate's established convention
/// (`PipelineV2FailureCause::InvalidConfiguration`) of returning a typed
/// error for a config that can never succeed, rather than panicking or
/// silently clamping.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum EnsembleV2ConfigError {
    /// `rmsd_threshold` must be `0.0` (pruning disabled) or a positive,
    /// finite value — never `NaN`, infinite, or negative. Carries the
    /// rejected value.
    InvalidRmsdThreshold(f64),
}

impl std::fmt::Display for EnsembleV2ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EnsembleV2ConfigError::InvalidRmsdThreshold(v) => {
                write!(
                    f,
                    "rmsd_threshold must be 0.0 or a positive finite value, got {v}"
                )
            }
        }
    }
}

impl std::error::Error for EnsembleV2ConfigError {}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// Whether [`embed_ensemble_v2`] ran every requested attempt or stopped early
/// because `ensemble_timeout_ms` was exhausted. Never inferred by a caller
/// comparing `attempts.len()` against `requested_count` — always read this
/// field directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum EnsembleTermination {
    /// Every one of `requested_count` attempts ran (`count == 0` also
    /// reports `Completed`, with zero attempts).
    Completed,
    /// Stopped before `requested_count` attempts ran because
    /// `ensemble_timeout_ms` was exhausted.
    TimedOut,
}

/// What became of one successful attempt's coordinates.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum ConformerDisposition {
    /// Retained in [`EnsembleV2Result::ensemble`] at this index.
    Kept { conformer_index: usize },
    /// Discarded as a near-duplicate of another kept conformer **from the
    /// same `actual_force_field_used` group** — pruning never crosses a
    /// force-field group boundary (see the module docs).
    PrunedAsDuplicate {
        /// `attempt_index` of the kept conformer this one duplicated.
        representative_attempt_index: usize,
        /// The RMSD (Å) that triggered pruning — always `< rmsd_threshold`.
        rmsd: f64,
        /// Whether symmetric (automorphism-aware) or plain Kabsch RMSD was
        /// used for this comparison (`config.use_symmetric_rmsd_pruning`).
        symmetric: bool,
    },
}

/// Per-attempt evidence for a successful embed, beyond the coordinates
/// themselves (which live in [`EnsembleV2Result::ensemble`] when
/// `disposition` is `Kept`).
#[derive(Debug, Clone)]
pub struct ConformerSuccess {
    /// `Some(energy)` iff the underlying `EnergyReport` is not
    /// `EnergyReport::None` AND the value is finite (see the module docs'
    /// energy-comparability section). `None` whenever no real force field
    /// ran, or the reported energy was `NaN`/infinite.
    pub energy: Option<f64>,
    pub actual_force_field_used: ForceFieldPolicy,
    /// `Some(reason)` iff this attempt actually fell back from MMFF94 to UFF
    /// under `ForceFieldPolicy::Mmff94WithUffFallback` (see
    /// `PolicyMinimizeResult::fallback_reason`'s own doc for why this, not
    /// `actual_force_field_used != requested`, is the correct check) —
    /// connects a kept/pruned conformer back to the same typed reason A1's
    /// failure ledger already reports for outright MMFF94 failures.
    pub fallback_reason: Option<ForceFieldBridgeError>,
    pub disposition: ConformerDisposition,
}

/// One embedding attempt's full outcome — nothing is dropped, success or
/// failure, matching `embed_pipeline_v2`'s own "carry as much diagnostic
/// evidence as possible" standard one level up. A pruned success is still
/// `Ok(..)` here — its `disposition` records what happened, never silently
/// merged away.
#[derive(Debug)]
pub struct ConformerAttempt {
    pub attempt_index: usize,
    pub seed: u64,
    pub outcome: Result<ConformerSuccess, PipelineV2Failure>,
}

/// Full result of [`embed_ensemble_v2`]. Never `Result`-wrapped for
/// per-molecule/per-attempt outcomes — an ensemble with zero kept conformers
/// (every attempt failed, was pruned, or `count == 0`) is a valid,
/// fully-diagnosable outcome, matching [`ConformerEnsemble::new`]'s own
/// "zero conformers is a normal state" convention; per-attempt typed
/// failures and dispositions already carry the diagnostic detail in
/// `attempts`. `embed_ensemble_v2` itself still returns `Result` — only for
/// a config that can never succeed regardless of input, never for a
/// per-molecule outcome.
pub struct EnsembleV2Result {
    /// Kept conformers only. Ordered group-by-group (first-seen
    /// `actual_force_field_used` order), ascending energy *within* each
    /// group — never a single cross-group energy sort. Insertion order
    /// within a group whenever no member of that group has usable energy.
    pub ensemble: ConformerEnsemble,
    /// Every attempt, success or failure, in `attempt_index` order.
    pub attempts: Vec<ConformerAttempt>,
    /// `true` iff kept conformers span more than one distinct
    /// `actual_force_field_used` value — see the module docs' energy-
    /// comparability section.
    pub mixed_force_field: bool,
    /// Whether every requested attempt ran, or `ensemble_timeout_ms` cut the
    /// run short. See [`EnsembleTermination`].
    pub termination: EnsembleTermination,
    /// `config.count` at call time, for comparing against `attempts.len()`
    /// without needing to have kept a copy of the config around.
    pub requested_count: usize,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Generate a conformer ensemble by calling [`embed_pipeline_v2`]
/// `config.count` times, once per deterministically derived seed, then
/// selecting kept representatives by ascending energy within each
/// `actual_force_field_used` group (see the module docs for why generation
/// and selection are strictly separate phases).
///
/// Returns `Err(EnsembleV2ConfigError)` only for a config that can never
/// succeed regardless of `mol` or how many attempts run (currently: an
/// invalid `rmsd_threshold`) — never for a per-molecule outcome, which is
/// always a `EnsembleV2Result` with per-attempt detail in `attempts`.
///
/// `count == 0` returns an empty, `Completed` ensemble with zero attempts
/// immediately, matching
/// [`crate::generate_conformer_ensemble_with_config`]'s existing
/// early-return convention.
pub fn embed_ensemble_v2(
    mol: &Molecule,
    config: &EnsembleV2Config,
) -> Result<EnsembleV2Result, EnsembleV2ConfigError> {
    let valid_rmsd_threshold = config.rmsd_threshold == 0.0
        || (config.rmsd_threshold.is_finite() && config.rmsd_threshold > 0.0);
    if !valid_rmsd_threshold {
        return Err(EnsembleV2ConfigError::InvalidRmsdThreshold(
            config.rmsd_threshold,
        ));
    }

    if config.count == 0 {
        return Ok(EnsembleV2Result {
            ensemble: ConformerEnsemble::new(mol.clone()),
            attempts: Vec::new(),
            mixed_force_field: false,
            termination: EnsembleTermination::Completed,
            requested_count: 0,
        });
    }

    // -----------------------------------------------------------------
    // Phase 1: generate every attempt. No pruning, no ordering decisions —
    // just run the pipeline and record what happened.
    // -----------------------------------------------------------------
    let overall_start = Instant::now();
    let mut raw: Vec<(u64, Result<RawSuccess, PipelineV2Failure>)> =
        Vec::with_capacity(config.count);
    let mut termination = EnsembleTermination::Completed;

    for i in 0..config.count {
        if let Some(budget) = config.ensemble_timeout_ms {
            // Fail closed on exactly `Some(0)` regardless of millisecond-resolution
            // rounding — mirrors `PipelineV2Config::total_timeout_ms`'s own
            // `check_timeout!` convention (see that macro's doc comment for the
            // real, intermittent-CI-failure history behind this exact check).
            if budget == 0 || overall_start.elapsed().as_millis() as u64 > budget {
                termination = EnsembleTermination::TimedOut;
                break;
            }
        }

        let seed = derive_attempt_seed(config.base_seed, i);
        let mut per_call = config.per_conformer.clone();
        per_call.embed.random_seed = seed;

        let outcome = embed_pipeline_v2(mol, &per_call).map(|result| {
            let energy = match &result.force_field.energy_after {
                EnergyReport::None => None,
                report => {
                    let total = report.total();
                    total.is_finite().then_some(total)
                }
            };
            RawSuccess {
                attempt_index: i,
                coords: result.coords,
                energy,
                actual_force_field_used: result.force_field.actual_force_field_used,
                fallback_reason: result.force_field.fallback_reason.clone(),
            }
        });

        raw.push((seed, outcome));
    }

    // -----------------------------------------------------------------
    // Phase 2: select. Group successes by actual_force_field_used
    // (first-seen order), sort each group by ascending energy, then
    // greedily dedup WITHIN each group only.
    // -----------------------------------------------------------------
    let successes: Vec<&RawSuccess> = raw.iter().filter_map(|(_, o)| o.as_ref().ok()).collect();
    let selection = select_ensemble(
        mol,
        &successes,
        config.rmsd_threshold,
        config.use_symmetric_rmsd_pruning,
    );

    let mut dispositions = selection.dispositions;
    let attempts = raw
        .into_iter()
        .enumerate()
        .map(|(attempt_index, (seed, outcome))| {
            let outcome = outcome.map(|raw_success| ConformerSuccess {
                energy: raw_success.energy,
                actual_force_field_used: raw_success.actual_force_field_used,
                fallback_reason: raw_success.fallback_reason,
                disposition: dispositions
                    .remove(&attempt_index)
                    .expect("every successful attempt has a disposition from select_ensemble"),
            });
            ConformerAttempt {
                attempt_index,
                seed,
                outcome,
            }
        })
        .collect();

    Ok(EnsembleV2Result {
        ensemble: selection.ensemble,
        attempts,
        mixed_force_field: selection.mixed_force_field,
        termination,
        requested_count: config.count,
    })
}

/// A successful attempt's data, before disposition is decided (phase 1
/// output / phase 2 input). Not part of the public API — [`ConformerSuccess`]
/// is the public per-attempt type, once `disposition` is known.
struct RawSuccess {
    attempt_index: usize,
    coords: Coords3D,
    energy: Option<f64>,
    actual_force_field_used: ForceFieldPolicy,
    fallback_reason: Option<ForceFieldBridgeError>,
}

struct SelectionResult {
    ensemble: ConformerEnsemble,
    /// `attempt_index -> disposition`, one entry per input success.
    dispositions: std::collections::HashMap<usize, ConformerDisposition>,
    mixed_force_field: bool,
}

/// Phase 2 in isolation: group `successes` by `actual_force_field_used`
/// (first-seen order), sort each group by ascending energy (attempts with no
/// usable energy last, tie-broken by `attempt_index`), then greedily dedup
/// within each group only. Never crosses a group boundary for either ranking
/// or pruning (see the module docs).
///
/// Factored out from [`embed_ensemble_v2`] specifically so this selection
/// logic is unit-testable against synthetic `RawSuccess`-shaped data,
/// without needing a real `embed_pipeline_v2` call per test case.
fn select_ensemble(
    mol: &Molecule,
    successes: &[&RawSuccess],
    rmsd_threshold: f64,
    use_symmetric_rmsd_pruning: bool,
) -> SelectionResult {
    let mut group_order: Vec<ForceFieldPolicy> = Vec::new();
    for s in successes {
        if !group_order.contains(&s.actual_force_field_used) {
            group_order.push(s.actual_force_field_used);
        }
    }

    let mut dispositions: std::collections::HashMap<usize, ConformerDisposition> =
        std::collections::HashMap::new();
    // (attempt_index, coords), already in final group-then-energy order —
    // built incrementally, one group at a time, in `group_order`.
    let mut kept_in_order: Vec<(usize, ForceFieldPolicy, Coords3D)> = Vec::new();

    for &policy in &group_order {
        let mut candidates: Vec<&RawSuccess> = successes
            .iter()
            .copied()
            .filter(|s| s.actual_force_field_used == policy)
            .collect();
        candidates.sort_by(|a, b| {
            a.energy
                .is_none()
                .cmp(&b.energy.is_none())
                .then_with(|| {
                    a.energy
                        .unwrap_or(f64::INFINITY)
                        .total_cmp(&b.energy.unwrap_or(f64::INFINITY))
                })
                .then_with(|| a.attempt_index.cmp(&b.attempt_index))
        });

        // A scratch ensemble holding only this group's kept representatives
        // so far, so `find_duplicate`/`find_duplicate_symmetric` can be
        // reused as-is; `group_attempt_indices[i]` is the attempt_index
        // behind `group_scratch`'s conformer `i`.
        let mut group_scratch = ConformerEnsemble::new(mol.clone());
        let mut group_attempt_indices: Vec<usize> = Vec::new();

        for candidate in candidates {
            let dup = if use_symmetric_rmsd_pruning {
                group_scratch.find_duplicate_symmetric(&candidate.coords, rmsd_threshold)
            } else {
                group_scratch.find_duplicate(&candidate.coords, rmsd_threshold)
            };

            match dup {
                Some((existing_idx, rmsd)) => {
                    dispositions.insert(
                        candidate.attempt_index,
                        ConformerDisposition::PrunedAsDuplicate {
                            representative_attempt_index: group_attempt_indices[existing_idx],
                            rmsd,
                            symmetric: use_symmetric_rmsd_pruning,
                        },
                    );
                }
                None => {
                    group_scratch
                        .add_conformer(candidate.coords.clone())
                        .expect("candidate coords come from the same mol, atom count must match");
                    group_attempt_indices.push(candidate.attempt_index);
                    kept_in_order.push((candidate.attempt_index, policy, candidate.coords.clone()));
                }
            }
        }
    }

    let mut ensemble = ConformerEnsemble::new(mol.clone());
    for (attempt_index, _, coords) in &kept_in_order {
        let idx = ensemble
            .add_conformer(coords.clone())
            .expect("kept coords come from the same mol, atom count must match");
        dispositions.insert(
            *attempt_index,
            ConformerDisposition::Kept {
                conformer_index: idx,
            },
        );
    }

    let mixed_force_field = {
        let mut seen: Vec<ForceFieldPolicy> = Vec::new();
        for (_, policy, _) in &kept_in_order {
            if !seen.contains(policy) {
                seen.push(*policy);
            }
        }
        seen.len() > 1
    };

    SelectionResult {
        ensemble,
        dispositions,
        mixed_force_field,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::Point3;
    use crate::minimize::ForceFieldPolicy;
    use chematic_core::AtomIdx;
    use chematic_smiles::parse;

    fn config_none(count: usize, seed: u64) -> EnsembleV2Config {
        EnsembleV2Config::new(
            PipelineV2Config::minimal(ForceFieldPolicy::None),
            count,
            seed,
        )
    }

    // -----------------------------------------------------------------------
    // embed_ensemble_v2 integration tests (real embed_pipeline_v2 calls)
    // -----------------------------------------------------------------------

    #[test]
    fn count_zero_returns_empty_immediately() {
        let mol = parse("CCC").unwrap();
        let result = embed_ensemble_v2(&mol, &config_none(0, 42)).unwrap();
        assert_eq!(result.ensemble.conformer_count(), 0);
        assert!(result.attempts.is_empty());
        assert!(!result.mixed_force_field);
        assert_eq!(result.termination, EnsembleTermination::Completed);
        assert_eq!(result.requested_count, 0);
    }

    #[test]
    fn same_base_seed_is_deterministic() {
        let mol = parse("CCCCCCCC").unwrap();
        let config = {
            let mut c = config_none(4, 20260801);
            c.per_conformer.embed.max_attempts = 1;
            c.rmsd_threshold = 0.0; // no pruning -- compare raw attempts
            c
        };
        let r1 = embed_ensemble_v2(&mol, &config).unwrap();
        let r2 = embed_ensemble_v2(&mol, &config).unwrap();
        assert_eq!(r1.ensemble.conformer_count(), r2.ensemble.conformer_count());
        for i in 0..r1.ensemble.conformer_count() {
            let c1 = r1.ensemble.get_conformer(i).unwrap();
            let c2 = r2.ensemble.get_conformer(i).unwrap();
            for a in 0..mol.atom_count() {
                let idx = AtomIdx(a as u32);
                assert_eq!(c1.get(idx), c2.get(idx), "conformer {i} atom {a} mismatch");
            }
        }
    }

    #[test]
    fn different_base_seeds_are_not_aliased() {
        let mol = parse("CCCCCCCC").unwrap();
        let mut c1 = config_none(1, 1);
        c1.per_conformer.embed.max_attempts = 1;
        let mut c2 = config_none(1, 2);
        c2.per_conformer.embed.max_attempts = 1;

        let r1 = embed_ensemble_v2(&mol, &c1).unwrap();
        let r2 = embed_ensemble_v2(&mol, &c2).unwrap();
        assert_eq!(r1.ensemble.conformer_count(), 1);
        assert_eq!(r2.ensemble.conformer_count(), 1);
        let a = r1.ensemble.get_conformer(0).unwrap();
        let b = r2.ensemble.get_conformer(0).unwrap();
        let mut any_diff = false;
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if a.get(idx) != b.get(idx) {
                any_diff = true;
            }
        }
        assert!(
            any_diff,
            "different base seeds must not produce aliased output"
        );
    }

    #[test]
    fn force_field_policy_none_never_reports_energy() {
        let mol = parse("CCCCCCCC").unwrap();
        let mut config = config_none(3, 7);
        config.per_conformer.embed.max_attempts = 1;
        config.rmsd_threshold = 0.0;
        let result = embed_ensemble_v2(&mol, &config).unwrap();
        assert!(!result.mixed_force_field);
        for attempt in &result.attempts {
            if let Ok(success) = &attempt.outcome {
                assert!(
                    success.energy.is_none(),
                    "ForceFieldPolicy::None must never report a real energy"
                );
            }
        }
    }

    #[test]
    fn every_attempt_recorded_success_or_failure() {
        let mol = parse("c1ccccc1").unwrap();
        let mut config = config_none(5, 100);
        config.per_conformer.embed.max_attempts = 1;
        let result = embed_ensemble_v2(&mol, &config).unwrap();
        assert_eq!(result.attempts.len(), 5);
        assert_eq!(result.termination, EnsembleTermination::Completed);
        assert_eq!(result.requested_count, 5);
        for (i, attempt) in result.attempts.iter().enumerate() {
            assert_eq!(attempt.attempt_index, i);
        }
    }

    #[test]
    fn zero_timeout_budget_stops_before_first_attempt() {
        let mol = parse("CCC").unwrap();
        let mut config = config_none(5, 1);
        config.ensemble_timeout_ms = Some(0);
        let result = embed_ensemble_v2(&mol, &config).unwrap();
        assert!(
            result.attempts.is_empty(),
            "a zero ensemble timeout must fail closed before any attempt runs"
        );
        assert_eq!(result.termination, EnsembleTermination::TimedOut);
        assert_eq!(result.requested_count, 5);
    }

    #[test]
    fn generous_timeout_completes_and_reports_completed() {
        let mol = parse("CCC").unwrap();
        let mut config = config_none(3, 1);
        config.ensemble_timeout_ms = Some(60_000);
        let result = embed_ensemble_v2(&mol, &config).unwrap();
        assert_eq!(result.termination, EnsembleTermination::Completed);
        assert_eq!(result.attempts.len(), 3);
    }

    /// A budget picked to very likely cut a moderately expensive run short,
    /// without asserting an exact attempt count (that would be a real-clock
    /// flaky test — this project has hit exactly that class of CI flake
    /// before, see `PipelineV2Config::total_timeout_ms`'s own doc history).
    /// Instead this asserts the *invariant* `embed_ensemble_v2` must uphold
    /// regardless of which way the timing happens to land: `TimedOut` implies
    /// fewer attempts than requested, `Completed` implies exactly as many as
    /// requested. The test can never fail from timing variance, only from a
    /// genuine violation of that contract.
    #[test]
    fn timeout_termination_is_consistent_with_attempts_len() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
        let mut per_conformer = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94BondAngleStrict);
        per_conformer.embed.max_attempts = 8;
        let mut config = EnsembleV2Config::new(per_conformer, 50, 1);
        config.ensemble_timeout_ms = Some(20);

        let result = embed_ensemble_v2(&mol, &config).unwrap();
        assert_eq!(result.requested_count, 50);
        match result.termination {
            EnsembleTermination::TimedOut => {
                assert!(result.attempts.len() < result.requested_count);
            }
            EnsembleTermination::Completed => {
                assert_eq!(result.attempts.len(), result.requested_count);
            }
        }
    }

    #[test]
    fn invalid_rmsd_threshold_is_rejected() {
        let mol = parse("CCC").unwrap();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1] {
            let mut config = config_none(1, 1);
            config.rmsd_threshold = bad;
            match embed_ensemble_v2(&mol, &config) {
                Err(EnsembleV2ConfigError::InvalidRmsdThreshold(_)) => {}
                Ok(_) => panic!("rmsd_threshold {bad} must be rejected"),
            }
        }
    }

    #[test]
    fn zero_rmsd_threshold_is_accepted_as_pruning_disabled() {
        let mol = parse("CCC").unwrap();
        let mut config = config_none(2, 1);
        config.rmsd_threshold = 0.0;
        assert!(embed_ensemble_v2(&mol, &config).is_ok());
    }

    #[test]
    fn symmetric_pruning_collapses_automorphism_only_duplicates() {
        // 1,1,1-trifluoroethane: with a fixed embed seed and noise-free DG (no
        // stochastic torsion perturbation available at this level -- ETKDG-style
        // noise lives in the legacy path, not embed_pipeline_v2), repeated
        // attempts at the SAME seed are identical, so this test instead directly
        // exercises `use_symmetric_rmsd_pruning`'s effect via two independently
        // seeded attempts that a real automorphism (F-atom relabelling) can make
        // symmetric-identical even though plain Kabsch may not.
        let mol = parse("FC(F)(F)C").unwrap();
        let mut symmetric_cfg = config_none(6, 555);
        symmetric_cfg.per_conformer.embed.max_attempts = 1;
        symmetric_cfg.rmsd_threshold = 2.0; // generous -- exercises pruning, not exact-zero matching
        symmetric_cfg.use_symmetric_rmsd_pruning = true;

        let mut plain_cfg = symmetric_cfg.clone();
        plain_cfg.use_symmetric_rmsd_pruning = false;

        let symmetric_result = embed_ensemble_v2(&mol, &symmetric_cfg).unwrap();
        let plain_result = embed_ensemble_v2(&mol, &plain_cfg).unwrap();
        // Symmetric pruning can only keep the same number of conformers or
        // fewer than plain Kabsch pruning at the same threshold -- it is a
        // strictly more aggressive (never less aggressive) duplicate check.
        assert!(
            symmetric_result.ensemble.conformer_count() <= plain_result.ensemble.conformer_count(),
            "symmetric pruning ({}) must never keep more conformers than plain Kabsch pruning ({})",
            symmetric_result.ensemble.conformer_count(),
            plain_result.ensemble.conformer_count()
        );
    }

    // -----------------------------------------------------------------------
    // select_ensemble unit tests (synthetic data -- no real embed_pipeline_v2
    // calls, per the code-review recommendation: these pin the exact
    // selection contract without needing to engineer a real pipeline
    // scenario for it).
    // -----------------------------------------------------------------------

    fn synth_coords(mol: &Molecule, first_atom_x: f64) -> Coords3D {
        // Distinct-but-close coordinates: only the first heavy atom moves, by
        // a small amount, so two synthetic coords are near-duplicates of each
        // other under a generous RMSD threshold but never bit-identical.
        let n = mol.atom_count();
        let mut c = Coords3D::new_zeroed(n);
        for i in 0..n {
            c.set(AtomIdx(i as u32), Point3::new(i as f64 * 1.5, 0.0, 0.0));
        }
        c.set(AtomIdx(0), Point3::new(first_atom_x, 0.0, 0.0));
        c
    }

    fn synth_far_coords(mol: &Molecule) -> Coords3D {
        let n = mol.atom_count();
        let mut c = Coords3D::new_zeroed(n);
        for i in 0..n {
            c.set(AtomIdx(i as u32), Point3::new(0.0, i as f64 * 500.0, 0.0));
        }
        c
    }

    #[test]
    fn select_ensemble_keeps_lowest_energy_representative_regardless_of_generation_order() {
        let mol = parse("CCCCC").unwrap();
        // Near-duplicate pair (small first-atom offset), high-energy one
        // listed FIRST (mimicking "generated first") -- the fix under test is
        // that generation order must not matter.
        let high_energy = RawSuccess {
            attempt_index: 0,
            coords: synth_coords(&mol, 0.0),
            energy: Some(100.0),
            actual_force_field_used: ForceFieldPolicy::Dreiding,
            fallback_reason: None,
        };
        let low_energy = RawSuccess {
            attempt_index: 1,
            coords: synth_coords(&mol, 0.05), // tiny offset -- near-duplicate of attempt 0
            energy: Some(-50.0),
            actual_force_field_used: ForceFieldPolicy::Dreiding,
            fallback_reason: None,
        };
        let successes = [&high_energy, &low_energy];

        let selection = select_ensemble(&mol, &successes, 1.0, true);
        assert_eq!(
            selection.ensemble.conformer_count(),
            1,
            "the pair must collapse to one representative"
        );
        assert!(!selection.mixed_force_field);

        match &selection.dispositions[&1] {
            ConformerDisposition::Kept { conformer_index } => assert_eq!(*conformer_index, 0),
            other => panic!("expected attempt 1 (lower energy) to be Kept, got {other:?}"),
        }
        match &selection.dispositions[&0] {
            ConformerDisposition::PrunedAsDuplicate {
                representative_attempt_index,
                ..
            } => {
                assert_eq!(
                    *representative_attempt_index, 1,
                    "attempt 0 (higher energy, generated first) must be pruned in favor of attempt 1"
                );
            }
            other => panic!("expected attempt 0 (higher energy) to be pruned, got {other:?}"),
        }
    }

    #[test]
    fn select_ensemble_never_prunes_across_force_field_groups() {
        let mol = parse("CCCCC").unwrap();
        // Same geometry (well within threshold of each other) but different
        // actual_force_field_used -- both must survive, never resolved by
        // "whichever came first."
        let uff_fallback = RawSuccess {
            attempt_index: 0,
            coords: synth_coords(&mol, 0.0),
            energy: Some(10.0),
            actual_force_field_used: ForceFieldPolicy::UffOnly,
            fallback_reason: Some(ForceFieldBridgeError::UnsupportedAtomType("X".to_string())),
        };
        let mmff94 = RawSuccess {
            attempt_index: 1,
            coords: synth_coords(&mol, 0.01),
            energy: Some(-5.0),
            actual_force_field_used: ForceFieldPolicy::Mmff94BondAngleStrict,
            fallback_reason: None,
        };
        let successes = [&uff_fallback, &mmff94];

        let selection = select_ensemble(&mol, &successes, 1.0, true);
        assert_eq!(
            selection.ensemble.conformer_count(),
            2,
            "geometrically-similar candidates from different force-field groups must both be kept"
        );
        assert!(selection.mixed_force_field);
        assert!(matches!(
            selection.dispositions[&0],
            ConformerDisposition::Kept { .. }
        ));
        assert!(matches!(
            selection.dispositions[&1],
            ConformerDisposition::Kept { .. }
        ));
    }

    #[test]
    fn select_ensemble_reports_pruning_provenance() {
        let mol = parse("CCCCC").unwrap();
        let representative = RawSuccess {
            attempt_index: 0,
            coords: synth_coords(&mol, 0.0),
            energy: Some(-10.0),
            actual_force_field_used: ForceFieldPolicy::Dreiding,
            fallback_reason: None,
        };
        let duplicate = RawSuccess {
            attempt_index: 1,
            coords: synth_coords(&mol, 0.02),
            energy: Some(-1.0), // higher energy -- must be the one pruned
            actual_force_field_used: ForceFieldPolicy::Dreiding,
            fallback_reason: None,
        };
        let successes = [&representative, &duplicate];

        let selection = select_ensemble(&mol, &successes, 1.0, false); // plain Kabsch this time
        match &selection.dispositions[&1] {
            ConformerDisposition::PrunedAsDuplicate {
                representative_attempt_index,
                rmsd,
                symmetric,
            } => {
                assert_eq!(*representative_attempt_index, 0);
                assert!(
                    *rmsd < 1.0 && *rmsd >= 0.0,
                    "rmsd {rmsd} must be within the threshold"
                );
                assert!(!*symmetric, "this call requested plain Kabsch pruning");
            }
            other => panic!("expected PrunedAsDuplicate, got {other:?}"),
        }
    }

    #[test]
    fn select_ensemble_distinct_conformers_both_kept_no_provenance_confusion() {
        let mol = parse("CCCCC").unwrap();
        let a = RawSuccess {
            attempt_index: 0,
            coords: synth_coords(&mol, 0.0),
            energy: Some(-10.0),
            actual_force_field_used: ForceFieldPolicy::Dreiding,
            fallback_reason: None,
        };
        let b = RawSuccess {
            attempt_index: 1,
            coords: synth_far_coords(&mol), // nowhere near `a`
            energy: Some(-1.0),
            actual_force_field_used: ForceFieldPolicy::Dreiding,
            fallback_reason: None,
        };
        let successes = [&a, &b];
        let selection = select_ensemble(&mol, &successes, 1.0, true);
        assert_eq!(selection.ensemble.conformer_count(), 2);
        assert!(matches!(
            selection.dispositions[&0],
            ConformerDisposition::Kept { .. }
        ));
        assert!(matches!(
            selection.dispositions[&1],
            ConformerDisposition::Kept { .. }
        ));
    }

    #[test]
    fn select_ensemble_no_energy_group_falls_back_to_attempt_order() {
        let mol = parse("CCCCC").unwrap();
        let first = RawSuccess {
            attempt_index: 0,
            coords: synth_coords(&mol, 0.0),
            energy: None,
            actual_force_field_used: ForceFieldPolicy::None,
            fallback_reason: None,
        };
        let second = RawSuccess {
            attempt_index: 1,
            coords: synth_far_coords(&mol),
            energy: None,
            actual_force_field_used: ForceFieldPolicy::None,
            fallback_reason: None,
        };
        let successes = [&first, &second];
        let selection = select_ensemble(&mol, &successes, 1.0, true);
        assert_eq!(selection.ensemble.conformer_count(), 2);
        // Both kept (distinct geometry); order among them is attempt order
        // since neither has energy.
        assert_eq!(
            selection.ensemble.get_conformer(0).unwrap(),
            &synth_coords(&mol, 0.0),
            "first-processed conformer occupies index 0 when no energy signal exists"
        );
    }
}
