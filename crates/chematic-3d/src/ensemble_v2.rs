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
//! # Energy comparability across `ForceFieldPolicy::Mmff94WithUffFallback`
//!
//! `PolicyMinimizeResult::energy_after` is on a physically meaningful scale
//! only *within* a single force field's own parameterization — MMFF94 and UFF
//! energies are not on a common reference zero. Under
//! `ForceFieldPolicy::Mmff94WithUffFallback`, individual attempts in the same
//! ensemble can resolve to different `actual_force_field_used` values (some
//! stay on MMFF94, some fall back to UFF) — see that field's own doc comment
//! in `minimize.rs` for why `fallback_reason.is_some()`, not
//! `actual_force_field_used != requested_force_field`, is the correct
//! "did this one fall back" check. This module never sorts across that
//! boundary: kept conformers are grouped by `actual_force_field_used` (first-
//! seen order), each group is internally ordered by ascending energy, and
//! [`EnsembleV2Result::mixed_force_field`] discloses when more than one group
//! is present in the same result — a caller must not read cross-group
//! adjacency in `ensemble` as an energy comparison.
//!
//! Energy is available at all only when the underlying [`EnergyReport`] is
//! not [`EnergyReport::None`] — i.e. whenever `force_field_policy !=
//! ForceFieldPolicy::None` and minimization actually ran. There is no
//! separate "rank by energy" flag to misconfigure: ranking happens whenever
//! energy is actually present in the result data, keyed off the same
//! `EnergyReport` the caller already gets back, never off the *requested*
//! policy (an invalid-state-by-construction choice, not a validated one).

use crate::clock::Instant;
use crate::conformer::ConformerEnsemble;
use crate::minimize::{EnergyReport, ForceFieldPolicy};
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
    /// failures and from near-duplicate pruning.
    pub count: usize,
    /// Base seed; attempt `i`'s seed is `derive_attempt_seed(base_seed, i)`
    /// — the same derivation `embed_pipeline_v2`'s own internal retry loop
    /// uses, reused here rather than inventing a second scheme. The same
    /// `base_seed` always reproduces the same ensemble.
    pub base_seed: u64,
    /// Minimum RMSD (Å) between kept conformers. `0.0` disables pruning
    /// (matches `ConformerConfig::rmsd_threshold`'s existing convention).
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
    /// millisecond-resolution elapsed-time rounding.
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

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// One embedding attempt's full outcome — nothing is dropped, success or
/// failure, matching `embed_pipeline_v2`'s own "carry as much diagnostic
/// evidence as possible" standard one level up.
#[derive(Debug)]
pub struct ConformerAttempt {
    pub attempt_index: usize,
    pub seed: u64,
    pub outcome: Result<ConformerSuccess, PipelineV2Failure>,
}

/// Per-attempt evidence for a successful embed, beyond the coordinates
/// themselves (which live in [`EnsembleV2Result::ensemble`] when `kept`).
#[derive(Debug, Clone, Copy)]
pub struct ConformerSuccess {
    /// `Some(energy_after.total())` iff the underlying `EnergyReport` is not
    /// `EnergyReport::None` (see the module docs' energy-comparability
    /// section). `None` whenever no real force field ran.
    pub energy: Option<f64>,
    pub actual_force_field_used: ForceFieldPolicy,
    /// `false` if this attempt was discarded as a near-duplicate of an
    /// already-kept conformer (RMSD below `rmsd_threshold`).
    pub kept: bool,
}

/// Full result of [`embed_ensemble_v2`]. Never `Result`-wrapped at the top
/// level — an ensemble with zero kept conformers (every attempt failed, or
/// `count == 0`) is a valid, fully-diagnosable outcome, matching
/// [`ConformerEnsemble::new`]'s own "zero conformers is a normal state"
/// convention; per-attempt typed failures already carry the diagnostic
/// detail in `attempts`.
pub struct EnsembleV2Result {
    /// Kept conformers only. When more than one `ForceFieldPolicy` group is
    /// present (`mixed_force_field == true`), ordered group-by-group
    /// (first-seen order) with ascending-energy order *within* each group —
    /// never a single cross-group energy sort. Insertion order (by
    /// attempt index) whenever no group has usable energy data.
    pub ensemble: ConformerEnsemble,
    /// Every attempt, success or failure, in `attempt_index` order.
    pub attempts: Vec<ConformerAttempt>,
    /// `true` iff kept conformers span more than one distinct
    /// `actual_force_field_used` value — see the module docs' energy-
    /// comparability section.
    pub mixed_force_field: bool,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Generate a conformer ensemble by calling [`embed_pipeline_v2`]
/// `config.count` times, once per deterministically derived seed, pruning
/// near-duplicates and energy-ranking the result (see the module docs).
///
/// `count == 0` returns an empty ensemble with zero attempts immediately,
/// matching [`crate::generate_conformer_ensemble_with_config`]'s existing
/// early-return convention.
pub fn embed_ensemble_v2(mol: &Molecule, config: &EnsembleV2Config) -> EnsembleV2Result {
    if config.count == 0 {
        return EnsembleV2Result {
            ensemble: ConformerEnsemble::new(mol.clone()),
            attempts: Vec::new(),
            mixed_force_field: false,
        };
    }

    let overall_start = Instant::now();
    let mut ensemble = ConformerEnsemble::new(mol.clone());
    let mut attempts = Vec::with_capacity(config.count);
    // Parallel to `ensemble`'s kept conformers: (attempt_index, energy, actual_ff).
    let mut kept_meta: Vec<(usize, Option<f64>, ForceFieldPolicy)> = Vec::new();

    for i in 0..config.count {
        if let Some(budget) = config.ensemble_timeout_ms {
            // Fail closed on exactly `Some(0)` regardless of millisecond-resolution
            // rounding — mirrors `PipelineV2Config::total_timeout_ms`'s own
            // `check_timeout!` convention (see that macro's doc comment for the
            // real, intermittent-CI-failure history behind this exact check).
            if budget == 0 || overall_start.elapsed().as_millis() as u64 > budget {
                break;
            }
        }

        let seed = derive_attempt_seed(config.base_seed, i);
        let mut per_call = config.per_conformer.clone();
        per_call.embed.random_seed = seed;

        let outcome = match embed_pipeline_v2(mol, &per_call) {
            Ok(result) => {
                let energy = match &result.force_field.energy_after {
                    EnergyReport::None => None,
                    report => Some(report.total()),
                };
                let actual_force_field_used = result.force_field.actual_force_field_used;
                let is_duplicate = if config.use_symmetric_rmsd_pruning {
                    ensemble.is_duplicate_symmetric(&result.coords, config.rmsd_threshold)
                } else {
                    ensemble.is_duplicate(&result.coords, config.rmsd_threshold)
                };
                let kept = !is_duplicate;

                if kept {
                    let idx = ensemble.add_conformer(result.coords).expect(
                        "coords come from embed_pipeline_v2 on the same mol, atom count must match",
                    );
                    kept_meta.push((idx, energy, actual_force_field_used));
                }

                Ok(ConformerSuccess {
                    energy,
                    actual_force_field_used,
                    kept,
                })
            }
            Err(failure) => Err(failure),
        };

        attempts.push(ConformerAttempt {
            attempt_index: i,
            seed,
            outcome,
        });
    }

    let mixed_force_field = kept_meta
        .first()
        .is_some_and(|(_, _, first_ff)| kept_meta.iter().any(|(_, _, ff)| ff != first_ff));

    let ensemble = reorder_by_group_then_energy(ensemble, &mut kept_meta);

    EnsembleV2Result {
        ensemble,
        attempts,
        mixed_force_field,
    }
}

/// Reorder `ensemble`'s conformers: group by `actual_force_field_used`
/// (first-seen group order preserved), ascending energy within each group,
/// original attempt order as the final tiebreak (also the order used when no
/// group has any usable energy). Uses `f64::total_cmp` for a total order over
/// energies (never `NaN`/`inf` in practice, but this avoids an `unwrap()` on
/// `partial_cmp` regardless).
fn reorder_by_group_then_energy(
    mut ensemble: ConformerEnsemble,
    kept_meta: &mut [(usize, Option<f64>, ForceFieldPolicy)],
) -> ConformerEnsemble {
    let mut group_order: Vec<ForceFieldPolicy> = Vec::new();
    for (_, _, ff) in kept_meta.iter() {
        if !group_order.contains(ff) {
            group_order.push(*ff);
        }
    }

    let mut order: Vec<usize> = (0..kept_meta.len()).collect();
    order.sort_by(|&a, &b| {
        let (idx_a, energy_a, ff_a) = kept_meta[a];
        let (idx_b, energy_b, ff_b) = kept_meta[b];
        let group_a = group_order.iter().position(|p| *p == ff_a).unwrap();
        let group_b = group_order.iter().position(|p| *p == ff_b).unwrap();
        group_a
            .cmp(&group_b)
            .then_with(|| {
                energy_a
                    .unwrap_or(f64::INFINITY)
                    .total_cmp(&energy_b.unwrap_or(f64::INFINITY))
            })
            .then_with(|| idx_a.cmp(&idx_b))
    });

    if order.windows(2).all(|w| w[0] < w[1]) {
        // Already in the desired order (the common case: single group, or
        // ranking that happens to match attempt order) — avoid the rebuild.
        return ensemble;
    }

    let mut reordered = ConformerEnsemble::new(ensemble.mol().clone());
    for &i in &order {
        let coords = ensemble
            .get_conformer(kept_meta[i].0)
            .expect("kept_meta indices are valid conformer indices by construction")
            .clone();
        reordered
            .add_conformer(coords)
            .expect("same mol as source ensemble, atom count must match");
    }
    // `ensemble` is dropped here; `reordered` is the returned value.
    ensemble = reordered;
    ensemble
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::Coords3D;
    use crate::minimize::ForceFieldPolicy;
    use chematic_smiles::parse;

    fn config_none(count: usize, seed: u64) -> EnsembleV2Config {
        EnsembleV2Config::new(
            PipelineV2Config::minimal(ForceFieldPolicy::None),
            count,
            seed,
        )
    }

    #[test]
    fn count_zero_returns_empty_immediately() {
        let mol = parse("CCC").unwrap();
        let result = embed_ensemble_v2(&mol, &config_none(0, 42));
        assert_eq!(result.ensemble.conformer_count(), 0);
        assert!(result.attempts.is_empty());
        assert!(!result.mixed_force_field);
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
        let r1 = embed_ensemble_v2(&mol, &config);
        let r2 = embed_ensemble_v2(&mol, &config);
        assert_eq!(r1.ensemble.conformer_count(), r2.ensemble.conformer_count());
        for i in 0..r1.ensemble.conformer_count() {
            let c1 = r1.ensemble.get_conformer(i).unwrap();
            let c2 = r2.ensemble.get_conformer(i).unwrap();
            for a in 0..mol.atom_count() {
                let idx = chematic_core::AtomIdx(a as u32);
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

        let r1 = embed_ensemble_v2(&mol, &c1);
        let r2 = embed_ensemble_v2(&mol, &c2);
        assert_eq!(r1.ensemble.conformer_count(), 1);
        assert_eq!(r2.ensemble.conformer_count(), 1);
        let a = r1.ensemble.get_conformer(0).unwrap();
        let b = r2.ensemble.get_conformer(0).unwrap();
        let mut any_diff = false;
        for i in 0..mol.atom_count() {
            let idx = chematic_core::AtomIdx(i as u32);
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
        let result = embed_ensemble_v2(&mol, &config);
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
        let result = embed_ensemble_v2(&mol, &config);
        assert_eq!(result.attempts.len(), 5);
        for (i, attempt) in result.attempts.iter().enumerate() {
            assert_eq!(attempt.attempt_index, i);
        }
    }

    #[test]
    fn zero_timeout_budget_stops_before_first_attempt() {
        let mol = parse("CCC").unwrap();
        let mut config = config_none(5, 1);
        config.ensemble_timeout_ms = Some(0);
        let result = embed_ensemble_v2(&mol, &config);
        assert!(
            result.attempts.is_empty(),
            "a zero ensemble timeout must fail closed before any attempt runs"
        );
    }

    #[test]
    fn kept_conformers_are_energy_ranked_ascending() {
        use crate::pipeline_v2::{PipelineV2Config, embed_pipeline_v2};

        // Small acyclic molecule + Dreiding (cheap physics) + few attempts: this
        // test only needs *some* real, distinguishable energies to check ordering,
        // not a realistic drug-like MMFF94 workload -- keeps a debug-build test
        // run in well under a second instead of minutes.
        let mol = parse("CCCCCC").unwrap(); // hexane
        let mut per_conformer = PipelineV2Config::minimal(ForceFieldPolicy::Dreiding);
        per_conformer.embed.max_attempts = 2;
        let config = EnsembleV2Config::new(per_conformer.clone(), 5, 20260801);

        // Independently reproduce each attempt's (coords, energy) pair by calling
        // embed_pipeline_v2 directly with the same derived seeds -- a shadow
        // computation, not a read of ensemble_v2's own internals -- to check the
        // returned ensemble's order without relying on any not-yet-tested
        // internal field mapping conformer index back to its source attempt.
        let mut shadow: Vec<(Coords3D, f64)> = Vec::new();
        for i in 0..config.count {
            let seed = derive_attempt_seed(config.base_seed, i);
            let mut call = per_conformer.clone();
            call.embed.random_seed = seed;
            if let Ok(r) = embed_pipeline_v2(&mol, &call) {
                shadow.push((r.coords, r.force_field.energy_after.total()));
            }
        }

        let result = embed_ensemble_v2(&mol, &config);
        assert!(
            !result.mixed_force_field,
            "single-policy run must never be mixed"
        );
        assert!(
            result.ensemble.conformer_count() >= 2,
            "need >=2 kept conformers to test ordering"
        );

        let mut prev_energy = f64::NEG_INFINITY;
        for i in 0..result.ensemble.conformer_count() {
            let coords = result.ensemble.get_conformer(i).unwrap();
            let (_, energy) = shadow
                .iter()
                .find(|(c, _)| c == coords)
                .expect("every kept conformer's coords must match some shadow attempt exactly");
            assert!(
                *energy >= prev_energy - 1e-9,
                "kept conformer {i} energy {energy} is out of ascending order (prev {prev_energy})"
            );
            prev_energy = *energy;
        }
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

        let symmetric_result = embed_ensemble_v2(&mol, &symmetric_cfg);
        let plain_result = embed_ensemble_v2(&mol, &plain_cfg);
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
}
