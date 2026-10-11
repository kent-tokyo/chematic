// --- policy dispatch ---------------------------------------------------------

/// Residual-force threshold (kcal/mol/Å) at which the MMFF94 bridge's
/// minimizer reports convergence; the same absolute test chematic-ff's L-BFGS
/// applies.
const MMFF94_BRIDGE_GRADIENT_TOLERANCE: f64 = 1e-4;

/// Which chematic-ff minimizer the MMFF94 bridge runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Mmff94Minimizer {
    /// Dense BFGS with RDKit's line search (`Mmff94EnergyModel::minimize_bfgs`),
    /// the default since batch 13.
    Bfgs,
    /// The earlier L-BFGS (`minimize_lbfgs_bounded_analytic`); the pipeline
    /// retries with it when the BFGS geometry fails declared stereo.
    Lbfgs,
}

struct Mmff94BridgeRun {
    coords: Coords3D,
    coverage: Mmff94CoverageReport,
    energy_before: EnergyBreakdown,
    energy_after: EnergyBreakdown,
    converged: bool,
    iterations: usize,
    termination: Mmff94TerminationReason,
    max_residual_force: f64,
}

fn run_mmff94_bridge(
    mol: &Molecule,
    coords: &Coords3D,
    max_iter: usize,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
    accept_geometry: Option<&dyn Fn(&Coords3D) -> bool>,
    minimizer: Mmff94Minimizer,
) -> Result<Mmff94BridgeRun, ForceFieldBridgeError> {
    let n = mol.atom_count();
    // Must use the same MMFF-specific re-perceived bond orders chematic-ff's
    // own production energy functions use internally (issue #227 Phase 1,
    // torsion parameter gap root cause) -- otherwise this coverage gate and
    // the energy/gradient functions it gates disagree on which classification
    // code a bond/angle/torsion/stretch-bend term resolves to.
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let coverage = compute_mmff94_coverage(&mmff_mol, &types);
    if coverage.has_gate_failure(include_torsion_oop_in_gate, include_stretch_bend_in_gate) {
        return Err(ForceFieldBridgeError::MissingParameters(Box::new(coverage)));
    }

    let coord_vec = coords_to_vec(coords, n);
    let energy_model = Mmff94EnergyModel::new(mol)?;
    let energy_before = energy_model.energy_breakdown(&coord_vec);

    let mut work = coord_vec.clone();
    let convergence = Mmff94Convergence::MaxGradient(MMFF94_BRIDGE_GRADIENT_TOLERANCE);
    let result = match (minimizer, accept_geometry) {
        (Mmff94Minimizer::Bfgs, Some(accept_geometry)) => energy_model
            .minimize_bfgs_with_constraint(&mut work, max_iter, convergence, |candidate| {
                accept_geometry(&vec_to_coords(candidate))
            })?,
        (Mmff94Minimizer::Bfgs, None) => {
            energy_model.minimize_bfgs(&mut work, max_iter, convergence)?
        }
        (Mmff94Minimizer::Lbfgs, Some(accept_geometry)) => energy_model
            .minimize_lbfgs_bounded_analytic_with_constraint(&mut work, max_iter, |candidate| {
                accept_geometry(&vec_to_coords(candidate))
            })?,
        (Mmff94Minimizer::Lbfgs, None) => {
            energy_model.minimize_lbfgs_bounded_analytic(&mut work, max_iter)?
        }
    };

    let energy_after = energy_model.energy_breakdown(&work);
    // The minimizer and the reported energy use this exact prepared MMFF94
    // objective. Re-running a central-difference gradient here cost 6N full
    // energy evaluations after every minimization and disproportionately
    // penalized the stereo-safe explicit-H lane. The analytic gradient is
    // already the production minimizer gradient and is checked against the
    // same total energy in chematic-ff; one final analytic evaluation preserves
    // the residual-force soundness gate without the redundant O(N) objective
    // sweep.
    let max_residual_force = max_gradient_component(&energy_model.bounded_analytic_gradient(&work));

    check_minimization_soundness(
        mol,
        &work,
        ForceFieldPolicy::Mmff94BondAngleStrict,
        OptimizerStatus {
            converged: result.converged,
            iterations: result.iterations,
            mmff94_termination: Some(result.termination),
        },
        max_residual_force,
        false,
    )?;

    Ok(Mmff94BridgeRun {
        coords: vec_to_coords(&work),
        coverage,
        energy_before,
        energy_after,
        converged: result.converged,
        iterations: result.iterations,
        termination: result.termination,
        max_residual_force,
    })
}

struct UffBridgeRun {
    coords: Coords3D,
    energy_before: f64,
    energy_after: f64,
    converged: bool,
    iterations: usize,
    max_residual_force: f64,
    starting_geometry: UffStartingGeometry,
}

/// Issues #185/#188 (UFF minimizer catastrophic bond-length blowup on some
/// starting geometries, sound convergence on others, at the shared default
/// 200-iteration budget): tries the caller-provided `coords` first. A retry
/// occurs for a catastrophic/non-finite result, a declared-stereo violation, or when
/// chematic-ff reports that it rejected an energy-decreasing but unsound
/// proposal during line search. The latter preserves the distinction between
/// an ordinary residual-force failure and a failure caused by the bounded
/// blowup guard. The retry uses `embed_distance_geometry_v2`'s geometry (see
/// [`rescue_with_distance_geometry_v2`] for why this is a post-hoc retry on
/// the actual outcome, not a pre-minimization heuristic guess) — every other
/// caller
/// caller (including every molecule that already passes today) sees zero
/// behavior change, and this never silently substitutes: which geometry
/// actually produced the returned result is always visible on
/// [`UffBridgeRun::starting_geometry`]/[`PolicyMinimizeResult::starting_geometry`],
/// and a failed rescue attempt is disclosed on
/// [`MinimizationFailureDetail::distance_geometry_v2_retry_attempted`].
fn run_uff_bridge(
    mol: &Molecule,
    coords: &Coords3D,
    max_iter: usize,
) -> Result<UffBridgeRun, ForceFieldBridgeError> {
    let n = mol.atom_count();
    let types = assign_uff_types(mol);
    let coord_vec = coords_to_vec(coords, n);
    let energy_before = uff_total_energy(mol, &types, &coord_vec);
    let result = ff_minimize_uff(mol, &types, coord_vec, max_iter);
    let energy_after = uff_total_energy(mol, &types, &result.coords);
    let max_residual_force =
        fd_max_gradient(&result.coords, |c| uff_total_energy(mol, &types, c), 1e-4);

    let first_coords = vec_to_coords(&result.coords);
    let first_preserves_stereo = !crate::distance_geometry_v2::mol_has_declared_stereo(mol)
        || crate::stereo_constraints::verify_stereo(mol, &first_coords).is_fully_satisfied();
    match check_minimization_soundness(
        mol,
        &result.coords,
        ForceFieldPolicy::UffOnly,
        OptimizerStatus {
            converged: result.converged,
            iterations: result.iterations,
            mmff94_termination: None,
        },
        max_residual_force,
        false,
    ) {
        Ok(()) if first_preserves_stereo => Ok(UffBridgeRun {
            coords: first_coords,
            energy_before,
            energy_after,
            converged: result.converged,
            iterations: result.iterations,
            max_residual_force,
            starting_geometry: UffStartingGeometry::AsProvided,
        }),
        Ok(()) => rescue_with_distance_geometry_v2(
            mol,
            &types,
            max_iter,
            Box::new(MinimizationFailureDetail {
                policy: ForceFieldPolicy::UffOnly,
                reason: MinimizationFailureReason::StereoConstraintViolation,
                converged: result.converged,
                iterations: result.iterations,
                mmff94_termination: None,
                max_residual_force,
                worst_bond_length: worst_bond_length_vec(mol, &result.coords),
                distance_geometry_v2_retry_attempted: false,
            }),
        ),
        Err(ForceFieldBridgeError::MinimizationFailed(detail))
            if !first_preserves_stereo
                || result.rejected_unsound_step
                || matches!(
                    detail.reason,
                    MinimizationFailureReason::CatastrophicBondBlowup
                        | MinimizationFailureReason::NonFiniteCoordinates
                ) =>
        {
            rescue_with_distance_geometry_v2(mol, &types, max_iter, detail)
        }
        Err(other) => Err(other),
    }
}

/// The rescue itself: `original_failure` is `coords`' own (already-computed)
/// failure. Raw starting energy alone cannot decide in advance whether a
/// given start needs this — measured (`docs/rfcs/uff_robustness_diagnosis_185_188.md`):
/// naphthalene and anthracene start from virtually identical raw
/// `generate_coords` UFF energies (~1.256e5 vs ~1.267e5 kcal/mol) and
/// worst starting bonds (2.26 Å, identical), yet naphthalene's raw UFF
/// minimization reaches a sound geometry (worst bond 2.19 Å) while
/// anthracene's plateaus permanently unsound (worst bond 3.39 Å, still
/// above `MAX_SANE_BOND_LENGTH` even at 200,000 steps) — so the
/// decision is made from the ACTUAL post-minimization outcome, never a
/// pre-minimization heuristic guess.
///
/// A retry is only accepted as a real rescue if it is BOTH geometrically
/// sound AND preserves every declared stereocenter/E-Z bond
/// ([`crate::stereo_constraints::verify_stereo`]) — measured directly, not
/// assumed. UFF itself is chirality-blind, so a sound retry is given a bounded
/// existing `repair_stereo` reconciliation loop before the final check. A rescue that fixed
/// a bond-length blowup while silently destroying declared stereochemistry
/// is a worse outcome than the honest failure it replaced, so this bridge
/// still returns the original failure when repair and re-verification do not
/// succeed. The first-attempt path is held to the same stereo gate: a sound
/// but stereo-violating result triggers this rescue instead of being returned.
///
/// **Final issue #210 revision:** declared-stereo retries use an
/// `add_hydrogens` working molecule through embedding, constrained UFF, and
/// repair, then truncate to the caller's original heavy-atom count and verify
/// again. Keeping the real H through minimization is load-bearing for ring-
/// fused testosterone/cholesterol centers, where no heavy-atom bridge can be
/// reflected safely. The public result's coordinates, UFF energies, residual
/// force, and soundness gate are all recomputed on the original molecule after
/// truncation. A fixed three-seed search keeps the operation deterministic and
/// bounded. The focused issue harness verifies all five named molecules;
/// this is not a claim that arbitrary stereochemical topologies always converge.
///
/// If `embed_distance_geometry_v2` itself fails, or the retry is unsound or
/// stereo-violating, the ORIGINAL failure is returned unchanged except for
/// `distance_geometry_v2_retry_attempted` flipping to `true` — the rescue
/// never masks a real failure with a different, potentially-misleading one.
fn rescue_with_distance_geometry_v2(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    max_iter: usize,
    original_failure: Box<MinimizationFailureDetail>,
) -> Result<UffBridgeRun, ForceFieldBridgeError> {
    use crate::distance_geometry_v2::{EmbedParameters, mol_has_declared_stereo};

    let n = mol.atom_count();
    // See this function's final issue #210 revision for why
    // these are set instead of `EmbedParameters::default()`.
    let has_stereo = mol_has_declared_stereo(mol);
    // Try a small fixed set of independent starts. This is deterministic and
    // keeps the rescue cost bounded while avoiding dependence on one unlucky
    // stereo-preserving embedding basin.
    const RESCUE_SEED_OFFSETS: [u64; 3] = [0, 0x9E37_79B9_7F4A_7C15, 0xD1B5_4A32_D192_ED03];
    // Ring-fused tetrahedral centers may have no movable heavy-atom bridge.
    // Keep explicit hydrogens present through both embedding and UFF so the
    // stereo repair/constraint machinery has a one-atom substituent it can
    // move without breaking a ring. Heavy atoms retain their indices and the
    // final accepted geometry is truncated back to the caller's molecule.
    let expanded_mol_storage;
    let working_mol: &Molecule = if has_stereo {
        expanded_mol_storage = chematic_chem::add_hydrogens(mol);
        &expanded_mol_storage
    } else {
        mol
    };
    let expanded_types_storage;
    let working_types: &[(AtomIdx, UffType)] = if has_stereo {
        expanded_types_storage = assign_uff_types(working_mol);
        &expanded_types_storage
    } else {
        types
    };
    let working_n = working_mol.atom_count();
    let base_params = EmbedParameters {
        enforce_chirality: has_stereo,
        // Hydrogens are already materialized above and must remain present
        // through minimization rather than being truncated after embedding.
        materialize_implicit_h_for_chirality: false,
        ..EmbedParameters::default()
    };
    let context = UffRescueContext {
        original_mol: mol,
        original_types: types,
        working_mol,
        working_types,
        original_atom_count: n,
        working_atom_count: working_n,
        has_stereo,
        max_iter,
    };
    for offset in RESCUE_SEED_OFFSETS {
        let embed_params = EmbedParameters {
            random_seed: base_params.random_seed.wrapping_add(offset),
            ..base_params.clone()
        };
        if let Some(run) = try_uff_rescue_seed(&context, &embed_params) {
            return Ok(run);
        }
    }

    let mut detail = original_failure;
    detail.distance_geometry_v2_retry_attempted = true;
    Err(ForceFieldBridgeError::MinimizationFailed(detail))
}

struct UffRescueContext<'a> {
    original_mol: &'a Molecule,
    original_types: &'a [(AtomIdx, UffType)],
    working_mol: &'a Molecule,
    working_types: &'a [(AtomIdx, UffType)],
    original_atom_count: usize,
    working_atom_count: usize,
    has_stereo: bool,
    max_iter: usize,
}

fn try_uff_rescue_seed(
    context: &UffRescueContext<'_>,
    params: &crate::distance_geometry_v2::EmbedParameters,
) -> Option<UffBridgeRun> {
    use crate::distance_geometry_v2::{embed_distance_geometry_v2, truncate_coords};
    use crate::stereo_constraints::verify_stereo;

    let embedded = embed_distance_geometry_v2(context.working_mol, params).ok()?;
    let initial_returned = if context.has_stereo {
        truncate_coords(&embedded, context.original_atom_count)
    } else {
        embedded.clone()
    };
    let energy_before = uff_total_energy(
        context.original_mol,
        context.original_types,
        &coords_to_vec(&initial_returned, context.original_atom_count),
    );
    let retry = minimize_uff_rescue(context, &embedded);
    let retry_coords = vec_to_coords(&retry.coords);
    let accepted = reconcile_rescue_stereo(context.working_mol, retry_coords);
    let working_stereo_ok = verify_stereo(context.working_mol, &accepted).is_fully_satisfied();
    let returned = if context.has_stereo {
        truncate_coords(&accepted, context.original_atom_count)
    } else {
        accepted
    };
    let accepted_vec = coords_to_vec(&returned, context.original_atom_count);
    let residual = fd_max_gradient(
        &accepted_vec,
        |coords| uff_total_energy(context.original_mol, context.original_types, coords),
        1e-4,
    );
    let sound = check_minimization_soundness(
        context.original_mol,
        &accepted_vec,
        ForceFieldPolicy::UffOnly,
        OptimizerStatus {
            converged: retry.converged,
            iterations: retry.iterations,
            mmff94_termination: None,
        },
        residual,
        false,
    )
    .is_ok();
    let stereo_ok =
        working_stereo_ok && verify_stereo(context.original_mol, &returned).is_fully_satisfied();
    (sound && stereo_ok).then(|| UffBridgeRun {
        coords: returned,
        energy_before,
        energy_after: uff_total_energy(
            context.original_mol,
            context.original_types,
            &accepted_vec,
        ),
        converged: retry.converged,
        iterations: retry.iterations,
        max_residual_force: residual,
        starting_geometry: UffStartingGeometry::ReplacedWithDistanceGeometryV2,
    })
}

fn minimize_uff_rescue(
    context: &UffRescueContext<'_>,
    embedded: &Coords3D,
) -> chematic_ff::UffMinimizeResult {
    use crate::stereo_constraints::verify_stereo;
    let coords = coords_to_vec(embedded, context.working_atom_count);
    if context.has_stereo && verify_stereo(context.working_mol, embedded).is_fully_satisfied() {
        ff_minimize_uff_with_constraint(
            context.working_mol,
            context.working_types,
            coords,
            context.max_iter,
            |candidate| {
                verify_stereo(context.working_mol, &vec_to_coords(candidate)).is_fully_satisfied()
            },
        )
    } else {
        ff_minimize_uff(
            context.working_mol,
            context.working_types,
            coords,
            context.max_iter,
        )
    }
}

fn reconcile_rescue_stereo(mol: &Molecule, retry_coords: Coords3D) -> Coords3D {
    use crate::stereo_constraints::{repair_stereo, verify_stereo};
    let mut accepted = retry_coords.clone();
    for _ in 0..3 {
        if verify_stereo(mol, &accepted).is_fully_satisfied() {
            break;
        }
        accepted = match repair_stereo(mol, &accepted) {
            Ok(repaired) => repaired.coords,
            Err(failure) => failure.partial_coords,
        };
    }
    let accepted_vec = coords_to_vec(&accepted, mol.atom_count());
    if accepted_vec.iter().all(|point| point.iter().all(|x| x.is_finite()))
        && worst_bond_length_vec(mol, &accepted_vec) <= MAX_SANE_BOND_LENGTH
    {
        accepted
    } else {
        retry_coords
    }
}

fn finish_mmff94(
    r: Mmff94BridgeRun,
    requested_force_field: ForceFieldPolicy,
    actual_force_field_used: ForceFieldPolicy,
    fallback_reason: Option<ForceFieldBridgeError>,
) -> PolicyMinimizeResult {
    let missing_parameter_classes = r.coverage.all_missing();
    PolicyMinimizeResult {
        coords: r.coords,
        requested_force_field,
        actual_force_field_used,
        fallback_reason,
        missing_parameter_classes,
        coverage: Some(r.coverage),
        energy_before: EnergyReport::Mmff94(r.energy_before),
        energy_after: EnergyReport::Mmff94(r.energy_after),
        converged: r.converged,
        iterations: r.iterations,
        mmff94_termination: Some(r.termination),
        max_residual_force: r.max_residual_force,
        starting_geometry: None,
    }
}

fn finish_uff(
    r: UffBridgeRun,
    requested_force_field: ForceFieldPolicy,
    actual_force_field_used: ForceFieldPolicy,
    fallback_reason: Option<ForceFieldBridgeError>,
    missing_parameter_classes: Vec<Mmff94MissingTerm>,
    coverage: Option<Mmff94CoverageReport>,
) -> PolicyMinimizeResult {
    PolicyMinimizeResult {
        coords: r.coords,
        requested_force_field,
        actual_force_field_used,
        fallback_reason,
        missing_parameter_classes,
        coverage,
        energy_before: EnergyReport::Uff {
            total: r.energy_before,
        },
        energy_after: EnergyReport::Uff {
            total: r.energy_after,
        },
        converged: r.converged,
        iterations: r.iterations,
        mmff94_termination: None,
        max_residual_force: r.max_residual_force,
        starting_geometry: Some(r.starting_geometry),
    }
}

/// Bridge chematic-3d's minimization step to chematic-ff's complete
/// force-field implementations. Opt-in only — no existing public function's
/// default behavior changes because of this.
///
/// `include_torsion_oop_in_gate`: when `true`, `Mmff94BondAngleStrict` (and the MMFF94
/// attempt inside `Mmff94WithUffFallback`) also refuses on a missing torsion
/// or out-of-plane term, not just bond/angle. [`minimize_with_policy`] passes
/// `false` — bond+angle is the RFC's mechanism-3 scope (the "no restoring
/// force, vdW pushes atoms apart" failure mode); torsion/oop absence is a
/// real, always-*reported* accuracy gap (see `coverage`/
/// `missing_parameter_classes` on the result) but not the same
/// structural-blowup pathology. Pre-chematic-ff-#183, gating on it too was
/// measured to fail a large fraction of ordinary organic molecules
/// (`angle_type_for` was hardcoded to 0, permanently stranding ~21% of the
/// angle table, which also fed torsion typing). That measurement is now
/// stale: post-#183, `examples/mmff94_bridge_coverage_report.rs` measures
/// 32/58 passing the bond+angle-only gate vs. 31/58 under the widened gate
/// — nearly identical — so widening the default is now a cheap, small
/// follow-up rather than impractical, just not done in this PR (see
/// [`ForceFieldPolicy::Mmff94BondAngleStrict`]'s doc for why not).
///
/// `include_stretch_bend_in_gate`: same opt-in shape, independent of
/// `include_torsion_oop_in_gate` (Priority 2 / Stage 1B, issue #227) — when
/// `true`, also refuses on a missing stretch-bend cross term. Stretch-bend
/// was historically not even measured toward the gate (no field on
/// `Mmff94CoverageReport` existed for it at all until this parameter was
/// added); `stretch_bend_energy` (chematic-ff's `mmff94_minimizer`) silently
/// contributes zero energy for an uncovered stretch-bend term regardless of
/// this flag — this parameter only controls whether that silence is also a
/// typed refusal up front, it does not change the energy function itself.
/// [`minimize_with_policy`] passes `false`, matching its existing
/// `include_torsion_oop_in_gate = false` default — no existing caller's
/// behavior changes.
#[allow(clippy::too_many_arguments)]
fn minimize_with_policy_gated_impl(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
    accept_geometry: Option<&dyn Fn(&Coords3D) -> bool>,
    minimizer: Mmff94Minimizer,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    if mol.atom_count() <= 1 {
        return Ok(trivial_result(coords, policy));
    }

    match policy {
        ForceFieldPolicy::None => Ok(trivial_result(coords, ForceFieldPolicy::None)),
        ForceFieldPolicy::Dreiding => run_dreiding_policy(mol, coords, config),
        ForceFieldPolicy::UffOnly => {
            let r = run_uff_bridge(mol, &coords, config.max_steps)?;
            Ok(finish_uff(
                r,
                ForceFieldPolicy::UffOnly,
                ForceFieldPolicy::UffOnly,
                None,
                Vec::new(),
                None,
            ))
        }

        ForceFieldPolicy::Mmff94BondAngleStrict => run_strict_mmff_policy(
            mol,
            coords,
            config,
            include_torsion_oop_in_gate,
            include_stretch_bend_in_gate,
            accept_geometry,
            minimizer,
        ),
        ForceFieldPolicy::Mmff94WithUffFallback => run_fallback_mmff_policy(
            mol,
            coords,
            config,
            include_torsion_oop_in_gate,
            include_stretch_bend_in_gate,
            accept_geometry,
            minimizer,
        ),
    }
}

fn run_dreiding_policy(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    let dreiding_types = assign_dreiding_types(mol);
    let e_before = total_energy_dreiding(mol, &coords, &dreiding_types);
    let report = minimize_gradient_descent_reporting(mol, coords, config, |candidate| {
        total_energy_dreiding(mol, candidate, &dreiding_types)
    });
    let e_after = total_energy_dreiding(mol, &report.coords, &dreiding_types);
    check_minimization_soundness(
        mol,
        &coords_to_vec(&report.coords, mol.atom_count()),
        ForceFieldPolicy::Dreiding,
        OptimizerStatus {
            converged: report.converged,
            iterations: report.iterations,
            mmff94_termination: None,
        },
        report.final_max_grad,
        false,
    )?;
    Ok(PolicyMinimizeResult {
        coords: report.coords,
        requested_force_field: ForceFieldPolicy::Dreiding,
        actual_force_field_used: ForceFieldPolicy::Dreiding,
        fallback_reason: None,
        missing_parameter_classes: Vec::new(),
        coverage: None,
        starting_geometry: None,
        energy_before: EnergyReport::Dreiding { total: e_before },
        energy_after: EnergyReport::Dreiding { total: e_after },
        converged: report.converged,
        iterations: report.iterations,
        mmff94_termination: None,
        max_residual_force: report.final_max_grad,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_strict_mmff_policy(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
    include_torsion_oop: bool,
    include_stretch_bend: bool,
    accept_geometry: Option<&dyn Fn(&Coords3D) -> bool>,
    minimizer: Mmff94Minimizer,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    let run = run_mmff94_bridge(
        mol,
        &coords,
        config.max_steps,
        include_torsion_oop,
        include_stretch_bend,
        accept_geometry,
        minimizer,
    )?;
    Ok(finish_mmff94(
        run,
        ForceFieldPolicy::Mmff94BondAngleStrict,
        ForceFieldPolicy::Mmff94BondAngleStrict,
        None,
    ))
}

#[allow(clippy::too_many_arguments)]
fn run_fallback_mmff_policy(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
    include_torsion_oop: bool,
    include_stretch_bend: bool,
    accept_geometry: Option<&dyn Fn(&Coords3D) -> bool>,
    minimizer: Mmff94Minimizer,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    match run_mmff94_bridge(
        mol,
        &coords,
        config.max_steps,
        include_torsion_oop,
        include_stretch_bend,
        accept_geometry,
        minimizer,
    ) {
        Ok(run) => Ok(finish_mmff94(
            run,
            ForceFieldPolicy::Mmff94WithUffFallback,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            None,
        )),
        Err(reason) => finish_uff_fallback(mol, coords, config, reason),
    }
}

fn finish_uff_fallback(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
    reason: ForceFieldBridgeError,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    let (coverage, missing) = match &reason {
        ForceFieldBridgeError::MissingParameters(report) => {
            (Some((**report).clone()), report.all_missing())
        }
        ForceFieldBridgeError::UnsupportedAtomType(_)
        | ForceFieldBridgeError::MinimizationFailed(_) => (None, Vec::new()),
    };
    let run = run_uff_bridge(mol, &coords, config.max_steps)?;
    Ok(finish_uff(
        run,
        ForceFieldPolicy::Mmff94WithUffFallback,
        ForceFieldPolicy::UffOnly,
        Some(reason),
        missing,
        coverage,
    ))
}

/// Run the selected force field with optional MMFF94 parameter-coverage
/// gates. This low-level API preserves its historical unconstrained MMFF94
/// line search; the stereo-aware pipeline uses an internal constrained path.
pub fn minimize_with_policy_gated(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    minimize_with_policy_gated_impl(
        mol,
        coords,
        policy,
        config,
        include_torsion_oop_in_gate,
        include_stretch_bend_in_gate,
        None,
        Mmff94Minimizer::Bfgs,
    )
}

/// [`minimize_with_policy_gated`] with the MMFF94 minimizer chosen by the
/// pipeline.
pub(crate) fn minimize_with_policy_gated_using(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
    minimizer: Mmff94Minimizer,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    minimize_with_policy_gated_impl(
        mol,
        coords,
        policy,
        config,
        include_torsion_oop_in_gate,
        include_stretch_bend_in_gate,
        None,
        minimizer,
    )
}

/// Pipeline-only force-field dispatch with a caller-owned MMFF94 line-search
/// acceptance predicate. Public low-level callers retain the historical
/// unconstrained minimizer through [`minimize_with_policy_gated`].
#[allow(clippy::too_many_arguments)]
pub(crate) fn minimize_with_policy_gated_with_constraint<F>(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
    accept_geometry: &F,
    minimizer: Mmff94Minimizer,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError>
where
    F: Fn(&Coords3D) -> bool,
{
    minimize_with_policy_gated_impl(
        mol,
        coords,
        policy,
        config,
        include_torsion_oop_in_gate,
        include_stretch_bend_in_gate,
        Some(accept_geometry),
        minimizer,
    )
}

/// Convenience wrapper over [`minimize_with_policy_gated`] with
/// `include_torsion_oop_in_gate = false` and `include_stretch_bend_in_gate =
/// false` (bond+angle only — mechanism-3's exact scope). Use
/// `minimize_with_policy_gated` directly to widen the strict gate to
/// torsion/out-of-plane and/or stretch-bend.
pub fn minimize_with_policy(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    minimize_with_policy_gated(mol, coords, policy, config, false, false)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
