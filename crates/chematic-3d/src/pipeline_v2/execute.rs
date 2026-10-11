/// Run the opt-in v2 embedding pipeline. Constraint-stalled and final-stereo
/// failures may be retried from deterministic alternate seeds.
#[allow(clippy::result_large_err)]
pub fn embed_pipeline_v2(
    mol: &Molecule,
    config: &PipelineV2Config,
) -> Result<PipelineV2Result, PipelineV2Failure> {
    let first = match embed_pipeline_v2_from_seed(mol, config) {
        Ok(first) => first,
        Err(failure) if retries_final_stereo_violation(config, &failure) => {
            return reseed_after_final_stereo_violation(mol, config, failure);
        }
        Err(failure) => return Err(failure),
    };
    if !minimization_stalled_on_constraint(&first) {
        return Ok(first);
    }
    let mut used_ms = first.elapsed_ms_by_stage.total_ms;
    for k in 1..=STALLED_MINIMIZATION_RESEEDS {
        let Some(retry) = retry_config(config, k, used_ms) else {
            break;
        };
        match embed_pipeline_v2_from_seed(mol, &retry) {
            Ok(result) => {
                used_ms += result.elapsed_ms_by_stage.total_ms;
                if !minimization_stalled_on_constraint(&result) && result.final_validation.sound {
                    return Ok(result);
                }
            }
            Err(failure) => used_ms += failure.elapsed_ms_by_stage.total_ms,
        }
    }
    Ok(first)
}

const STALLED_MINIMIZATION_RESEEDS: u64 = 3;
const FINAL_STEREO_VIOLATION_RESEEDS: u64 = 7;

fn retry_config(config: &PipelineV2Config, k: u64, used_ms: u64) -> Option<PipelineV2Config> {
    let mut retry = config.clone();
    retry.embed.random_seed = config
        .embed
        .random_seed
        .wrapping_add(k.wrapping_mul(0x9E37_79B9_7F4A_7C15));
    if let Some(budget) = config.total_timeout_ms {
        retry.total_timeout_ms = Some(budget.checked_sub(used_ms)?);
    }
    Some(retry)
}

fn retries_final_stereo_violation(config: &PipelineV2Config, failure: &PipelineV2Failure) -> bool {
    config.stereo_policy == StereoPolicy::RepairAndVerify
        && matches!(failure.cause, PipelineV2FailureCause::FinalStereoViolation)
}

#[allow(clippy::result_large_err)]
fn reseed_after_final_stereo_violation(
    mol: &Molecule,
    config: &PipelineV2Config,
    first: PipelineV2Failure,
) -> Result<PipelineV2Result, PipelineV2Failure> {
    let mut used_ms = first.elapsed_ms_by_stage.total_ms;
    for k in 1..=FINAL_STEREO_VIOLATION_RESEEDS {
        let Some(retry) = retry_config(config, k, used_ms) else {
            break;
        };
        match embed_pipeline_v2_from_seed(mol, &retry) {
            Ok(result) => return Ok(result),
            Err(failure) => used_ms += failure.elapsed_ms_by_stage.total_ms,
        }
    }
    Err(first)
}

fn minimization_stalled_on_constraint(result: &PipelineV2Result) -> bool {
    result.force_field.mmff94_termination
        == Some(chematic_ff::Mmff94TerminationReason::ConstraintRejectedFallback)
}

struct PipelineRun {
    started: Instant,
    budget_ms: Option<u64>,
    timings: StageTimings,
    evidence: Evidence,
}

impl PipelineRun {
    fn new(config: &PipelineV2Config) -> Self {
        Self {
            started: Instant::now(),
            budget_ms: config.total_timeout_ms,
            timings: StageTimings::default(),
            evidence: Evidence::default(),
        }
    }

    fn record_stage(&mut self, stage: PipelineStage, elapsed_ms: u64) {
        match stage {
            PipelineStage::ValidateConfig => {}
            PipelineStage::TorsionKnowledge => self.timings.torsion_knowledge_ms = elapsed_ms,
            PipelineStage::MacrocycleBoundAdjustment => {
                self.timings.bound_adjustment_ms = elapsed_ms;
            }
            PipelineStage::DistanceGeometry => self.timings.distance_geometry_ms = elapsed_ms,
            PipelineStage::TorsionEnergyEvaluation => {
                self.timings.torsion_energy_eval_ms = elapsed_ms;
            }
            PipelineStage::TorsionOptimization => {
                self.timings.torsion_optimization_ms = elapsed_ms;
            }
            PipelineStage::StereoVerifyBefore => {
                self.timings.stereo_verify_before_ms = elapsed_ms;
            }
            PipelineStage::StereoRepair => self.timings.stereo_repair_ms = elapsed_ms,
            PipelineStage::StereoVerifyAfterRepair => {
                self.timings.stereo_verify_after_repair_ms = elapsed_ms;
            }
            PipelineStage::ForceFieldMinimization => self.timings.force_field_ms = elapsed_ms,
            PipelineStage::FinalStereoVerify => {
                self.timings.final_stereo_verify_ms = elapsed_ms;
            }
            PipelineStage::FinalGeometryValidationStage => {
                self.timings.final_validation_ms = elapsed_ms;
            }
        }
    }

    #[allow(clippy::result_large_err)]
    fn complete_stage(
        &mut self,
        stage: PipelineStage,
        stage_started: Instant,
    ) -> Result<(), PipelineV2Failure> {
        self.record_stage(stage, stage_started.elapsed().as_millis() as u64);
        self.check_timeout(stage)
    }

    #[allow(clippy::result_large_err)]
    fn check_timeout(&mut self, stage: PipelineStage) -> Result<(), PipelineV2Failure> {
        let elapsed = self.started.elapsed().as_millis() as u64;
        if self
            .budget_ms
            .is_some_and(|budget| budget == 0 || elapsed > budget)
        {
            return Err(self.fail(PipelineV2FailureCause::Timeout, stage));
        }
        Ok(())
    }

    fn fail(&mut self, cause: PipelineV2FailureCause, stage: PipelineStage) -> PipelineV2Failure {
        self.timings.total_ms = self.started.elapsed().as_millis() as u64;
        self.evidence.fail(cause, stage, self.timings)
    }
}

struct KnowledgeStage {
    torsion: TorsionKnowledgeReport,
    bounds: Option<Vec<PairBoundAdjustment>>,
}

struct GeometryStage {
    coords: Coords3D,
    embed_stats: EmbedStats,
    ring_evidence: RingTorsionEvidence,
    torsion_optimization: Option<TorsionOptimizationReport>,
}

struct StereoStage {
    coords: Coords3D,
    before: StereoVerification,
    repair: Option<StereoRepairSummary>,
    after_repair: StereoVerification,
}

struct ForceFieldStage {
    result: PolicyMinimizeResult,
    minimizer: Mmff94Minimizer,
    reconcile_stereo: bool,
}

struct FinalStereoStage {
    coords: Coords3D,
    verification: StereoVerification,
    repair: Option<StereoRepairSummary>,
}

struct MoleculeContext<'a> {
    original: &'a Molecule,
    working: &'a Molecule,
    expanded: bool,
    original_atom_count: usize,
}

#[allow(clippy::result_large_err)]
fn embed_pipeline_v2_from_seed(
    mol: &Molecule,
    config: &PipelineV2Config,
) -> Result<PipelineV2Result, PipelineV2Failure> {
    let mut run = PipelineRun::new(config);
    validate_pipeline_config(config, &mut run)?;
    let knowledge = prepare_knowledge(mol, config, &mut run)?;
    let use_expanded = config.expand_implicit_h_through_pipeline
        && mol_needs_expanded_stereo_geometry(mol);
    let expanded_storage = use_expanded.then(|| expand_pipeline_molecule(mol));
    let molecules = MoleculeContext {
        original: mol,
        working: expanded_storage.as_ref().unwrap_or(mol),
        expanded: use_expanded,
        original_atom_count: mol.atom_count(),
    };
    let geometry = embed_and_optimize_geometry(&molecules, config, &knowledge, &mut run)?;
    let stereo = verify_and_repair_stereo(molecules.working, config, geometry.coords, &mut run)?;
    let force_field = run_force_field(&molecules, config, &stereo, &mut run)?;
    let final_stereo = finalize_stereo(&molecules, config, &force_field, &mut run)?;
    let final_validation = validate_final_geometry(
        molecules.working,
        &knowledge.torsion,
        geometry.torsion_optimization.as_ref(),
        &final_stereo,
        &mut run,
    )?;
    run.timings.total_ms = run.started.elapsed().as_millis() as u64;
    let coords = if molecules.expanded {
        truncate_coords(&final_stereo.coords, molecules.original_atom_count)
    } else {
        final_stereo.coords
    };
    Ok(PipelineV2Result {
        coords,
        embed_stats: geometry.embed_stats,
        bound_adjustment_report: knowledge.bounds,
        torsion_knowledge_report: knowledge.torsion,
        ring_torsion_evidence: geometry.ring_evidence,
        torsion_optimization_report: geometry.torsion_optimization,
        stereo_before: stereo.before,
        stereo_repair: stereo.repair,
        stereo_after_repair: stereo.after_repair,
        force_field: force_field.result,
        final_stereo: final_stereo.verification,
        post_minimization_stereo_repair: final_stereo.repair,
        final_validation,
        elapsed_ms_by_stage: run.timings,
    })
}

#[allow(clippy::result_large_err)]
fn validate_pipeline_config(
    config: &PipelineV2Config,
    run: &mut PipelineRun,
) -> Result<(), PipelineV2Failure> {
    if pipeline_config_is_invalid(config) {
        return Err(run.fail(
            PipelineV2FailureCause::InvalidConfiguration,
            PipelineStage::ValidateConfig,
        ));
    }
    Ok(())
}

fn expand_pipeline_molecule(mol: &Molecule) -> Molecule {
    if mol_needs_full_hydrogen_expansion(mol) {
        chematic_chem::add_hydrogens(mol)
    } else {
        chematic_chem::add_stereocenter_hydrogens(mol)
    }
}

#[allow(clippy::result_large_err)]
fn prepare_knowledge(
    mol: &Molecule,
    config: &PipelineV2Config,
    run: &mut PipelineRun,
) -> Result<KnowledgeStage, PipelineV2Failure> {
    let torsion_config = torsion_knowledge_config(config);
    let started = Instant::now();
    let torsion = build_torsion_knowledge(mol, &torsion_config);
    run.evidence.torsion_knowledge_report = Some(torsion.clone());
    run.complete_stage(PipelineStage::TorsionKnowledge, started)?;
    let started = Instant::now();
    let bounds = if config.embed.use_macrocycle_14_bounds {
        match macrocycle_14_bound_adjustments(mol, &torsion_config) {
            Ok(bounds) => Some(bounds),
            Err(error) => {
                run.record_stage(
                    PipelineStage::MacrocycleBoundAdjustment,
                    started.elapsed().as_millis() as u64,
                );
                return Err(run.fail(
                    PipelineV2FailureCause::TorsionKnowledge(error),
                    PipelineStage::MacrocycleBoundAdjustment,
                ));
            }
        }
    } else {
        None
    };
    run.evidence.bound_adjustment_report = bounds.clone();
    run.complete_stage(PipelineStage::MacrocycleBoundAdjustment, started)?;
    Ok(KnowledgeStage { torsion, bounds })
}

#[allow(clippy::result_large_err)]
fn embed_and_optimize_geometry(
    molecules: &MoleculeContext<'_>,
    config: &PipelineV2Config,
    knowledge: &KnowledgeStage,
    run: &mut PipelineRun,
) -> Result<GeometryStage, PipelineV2Failure> {
    let adjustments = knowledge
        .bounds
        .as_deref()
        .unwrap_or(&[])
        .iter()
        .map(|a| DistanceBoundAdjustment {
            atom1: a.atom_pair.0,
            atom2: a.atom_pair.1,
            lower: a.new_lower,
            upper: a.new_upper,
        })
        .collect::<Vec<_>>();
    let (coords, embed_stats) =
        run_distance_geometry(molecules.working, config, &adjustments, run)?;
    evaluate_initial_torsions(molecules.working, &coords, knowledge, run)?;
    let ring_evidence = build_ring_torsion_evidence(molecules.working, config, knowledge);
    run.evidence.ring_torsion_evidence = Some(ring_evidence.clone());
    enforce_ring_torsion_policy(config, knowledge, &ring_evidence, run)?;
    let (coords, torsion_optimization) =
        optimize_pipeline_torsions(molecules.working, config, knowledge, coords, run)?;
    Ok(GeometryStage {
        coords,
        embed_stats,
        ring_evidence,
        torsion_optimization,
    })
}

#[allow(clippy::result_large_err)]
fn run_distance_geometry(
    mol: &Molecule,
    config: &PipelineV2Config,
    adjustments: &[DistanceBoundAdjustment],
    run: &mut PipelineRun,
) -> Result<(Coords3D, EmbedStats), PipelineV2Failure> {
    let started = Instant::now();
    let result = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
        mol,
        &config.embed,
        adjustments,
    );
    let (coords, stats) = match result {
        Ok(result) => result,
        Err((kind, stats)) => {
            run.record_stage(
                PipelineStage::DistanceGeometry,
                started.elapsed().as_millis() as u64,
            );
            run.evidence.embed_stats = Some(stats);
            let (cause, stage) = match kind {
                EmbedWithAdjustmentsFailure::InvalidAdjustment => (
                    PipelineV2FailureCause::BoundAdjustmentFailed,
                    PipelineStage::MacrocycleBoundAdjustment,
                ),
                EmbedWithAdjustmentsFailure::Embed(cause) => (
                    PipelineV2FailureCause::DistanceGeometry(cause),
                    PipelineStage::DistanceGeometry,
                ),
            };
            return Err(run.fail(cause, stage));
        }
    };
    run.evidence.embed_stats = Some(stats.clone());
    run.evidence.last_known_coords = Some(coords.clone());
    run.complete_stage(PipelineStage::DistanceGeometry, started)?;
    Ok((coords, stats))
}

#[allow(clippy::result_large_err)]
fn evaluate_initial_torsions(
    mol: &Molecule,
    coords: &Coords3D,
    knowledge: &KnowledgeStage,
    run: &mut PipelineRun,
) -> Result<(), PipelineV2Failure> {
    let started = Instant::now();
    if let Err(error) = evaluate_torsion_energy(mol, coords, &knowledge.torsion.potentials) {
        run.record_stage(
            PipelineStage::TorsionEnergyEvaluation,
            started.elapsed().as_millis() as u64,
        );
        return Err(run.fail(
            PipelineV2FailureCause::TorsionKnowledge(error),
            PipelineStage::TorsionEnergyEvaluation,
        ));
    }
    run.complete_stage(PipelineStage::TorsionEnergyEvaluation, started)
}

fn build_ring_torsion_evidence(
    mol: &Molecule,
    config: &PipelineV2Config,
    knowledge: &KnowledgeStage,
) -> RingTorsionEvidence {
    RingTorsionEvidence {
        potentials: knowledge
            .torsion
            .potentials
            .iter()
            .map(|potential| PotentialApplicationEvidence {
                rule_id: potential.rule_id.clone(),
                central_bond: potential.central_bond,
                source: potential.source,
                applied_to_geometry: ring_torsion_is_applicable(mol, potential),
            })
            .collect(),
        diagnostic_only: config.ring_torsion_policy
            == RingTorsionApplicationPolicy::DiagnosticOnly,
    }
}

#[allow(clippy::result_large_err)]
fn enforce_ring_torsion_policy(
    config: &PipelineV2Config,
    knowledge: &KnowledgeStage,
    evidence: &RingTorsionEvidence,
    run: &mut PipelineRun,
) -> Result<(), PipelineV2Failure> {
    let requested = config.embed.use_small_ring_torsions || config.embed.use_macrocycle_torsions;
    let unsupported = knowledge
        .torsion
        .potentials
        .iter()
        .zip(&evidence.potentials)
        .any(|(potential, application)| {
            matches!(
                potential.source,
                TorsionKnowledgeSource::SmallRingExperimental
                    | TorsionKnowledgeSource::MacrocycleAdaptation
            ) && !application.applied_to_geometry
        });
    if requested
        && unsupported
        && config.ring_torsion_policy == RingTorsionApplicationPolicy::FailClosed
    {
        return Err(run.fail(
            PipelineV2FailureCause::RingTorsionApplicationUnsupported,
            PipelineStage::TorsionOptimization,
        ));
    }
    Ok(())
}

#[allow(clippy::result_large_err)]
fn optimize_pipeline_torsions(
    mol: &Molecule,
    config: &PipelineV2Config,
    knowledge: &KnowledgeStage,
    coords: Coords3D,
    run: &mut PipelineRun,
) -> Result<(Coords3D, Option<TorsionOptimizationReport>), PipelineV2Failure> {
    let started = Instant::now();
    let result = if knowledge.torsion.potentials.is_empty() {
        Ok((coords, None))
    } else {
        optimize_torsions(
            mol,
            &coords,
            &knowledge.torsion.potentials,
            &config.torsion_optimization,
        )
        .map(|(coords, report)| (coords, Some(report)))
    };
    let (coords, report) = match result {
        Ok(result) => result,
        Err(error) => {
            run.record_stage(
                PipelineStage::TorsionOptimization,
                started.elapsed().as_millis() as u64,
            );
            return Err(run.fail(
                PipelineV2FailureCause::TorsionKnowledge(error),
                PipelineStage::TorsionOptimization,
            ));
        }
    };
    run.evidence.torsion_optimization_report = report.clone();
    run.evidence.last_known_coords = Some(coords.clone());
    run.complete_stage(PipelineStage::TorsionOptimization, started)?;
    Ok((coords, report))
}

#[allow(clippy::result_large_err)]
fn verify_and_repair_stereo(
    mol: &Molecule,
    config: &PipelineV2Config,
    coords: Coords3D,
    run: &mut PipelineRun,
) -> Result<StereoStage, PipelineV2Failure> {
    let started = Instant::now();
    let before = verify_stereo(mol, &coords);
    run.evidence.stereo_before = Some(before.clone());
    run.complete_stage(PipelineStage::StereoVerifyBefore, started)?;
    let started = Instant::now();
    let (coords, repair) = repair_initial_stereo(mol, config, coords, run, started)?;
    run.evidence.stereo_repair = repair.clone();
    run.evidence.last_known_coords = Some(coords.clone());
    run.complete_stage(PipelineStage::StereoRepair, started)?;
    let started = Instant::now();
    let after_repair = repair
        .as_ref()
        .map_or_else(|| before.clone(), |_| verify_stereo(mol, &coords));
    run.evidence.stereo_after_repair = Some(after_repair.clone());
    run.complete_stage(PipelineStage::StereoVerifyAfterRepair, started)?;
    Ok(StereoStage {
        coords,
        before,
        repair,
        after_repair,
    })
}

#[allow(clippy::result_large_err)]
fn repair_initial_stereo(
    mol: &Molecule,
    config: &PipelineV2Config,
    coords: Coords3D,
    run: &mut PipelineRun,
    started: Instant,
) -> Result<(Coords3D, Option<StereoRepairSummary>), PipelineV2Failure> {
    if config.stereo_policy != StereoPolicy::RepairAndVerify {
        return Ok((coords, None));
    }
    match repair_stereo(mol, &coords) {
        Ok(outcome) => Ok((
            outcome.coords,
            Some(StereoRepairSummary {
                repaired: outcome.repaired,
                failures: Vec::new(),
            }),
        )),
        Err(report) => {
            run.record_stage(
                PipelineStage::StereoRepair,
                started.elapsed().as_millis() as u64,
            );
            run.evidence.last_known_coords = Some(report.partial_coords);
            run.evidence.stereo_repair = Some(StereoRepairSummary {
                repaired: report.repaired,
                failures: report.failures,
            });
            Err(run.fail(
                PipelineV2FailureCause::StereoRepairFailed,
                PipelineStage::StereoRepair,
            ))
        }
    }
}

fn authoritative_stereo(
    molecules: &MoleculeContext<'_>,
    coords: &Coords3D,
) -> StereoVerification {
    verify_authoritative_final_stereo(
        molecules.original,
        molecules.working,
        coords,
        molecules.expanded,
        molecules.original_atom_count,
    )
}

#[allow(clippy::too_many_arguments)]
fn minimize_force_field_once(
    molecules: &MoleculeContext<'_>,
    config: &PipelineV2Config,
    minimize_config: &MinimizeConfig,
    reconcile_stereo: bool,
    coords: Coords3D,
    minimizer: Mmff94Minimizer,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    if !reconcile_stereo {
        return minimize_with_policy_gated_using(
            molecules.working,
            coords,
            config.force_field_policy,
            minimize_config,
            config.gate_mmff94_torsion_oop,
            config.gate_mmff94_stretch_bend,
            minimizer,
        );
    }
    let preserves_stereo = |candidate: &Coords3D| {
        verify_stereo(molecules.working, candidate).is_fully_satisfied()
            && authoritative_stereo(molecules, candidate).n_unevaluable() == 0
    };
    minimize_with_policy_gated_with_constraint(
        molecules.working,
        coords,
        config.force_field_policy,
        minimize_config,
        config.gate_mmff94_torsion_oop,
        config.gate_mmff94_stretch_bend,
        &preserves_stereo,
        minimizer,
    )
}

#[allow(clippy::result_large_err)]
fn run_force_field(
    molecules: &MoleculeContext<'_>,
    config: &PipelineV2Config,
    stereo: &StereoStage,
    run: &mut PipelineRun,
) -> Result<ForceFieldStage, PipelineV2Failure> {
    let minimize_config = MinimizeConfig {
        max_steps: config.force_field_max_iterations,
        ..MinimizeConfig::default()
    };
    let reconcile_stereo = molecules.expanded
        && config.stereo_policy != StereoPolicy::Ignore
        && stereo.after_repair.is_fully_satisfied()
        && authoritative_stereo(molecules, &stereo.coords).n_violations() > 0;
    let started = Instant::now();
    let mut minimizer = Mmff94Minimizer::Bfgs;
    let mut result = minimize_force_field_once(
        molecules,
        config,
        &minimize_config,
        reconcile_stereo,
        stereo.coords.clone(),
        minimizer,
    );
    let mmff94_policy = matches!(
        config.force_field_policy,
        ForceFieldPolicy::Mmff94BondAngleStrict | ForceFieldPolicy::Mmff94WithUffFallback
    );
    if mmff94_policy
        && result.as_ref().map_or(true, |value| {
            config.stereo_policy != StereoPolicy::Ignore
                && authoritative_stereo(molecules, &value.coords).n_violations() > 0
        })
        && let Ok(retry) = minimize_force_field_once(
            molecules,
            config,
            &minimize_config,
            reconcile_stereo,
            stereo.coords.clone(),
            Mmff94Minimizer::Lbfgs,
        )
        && (result.is_err()
            || config.stereo_policy == StereoPolicy::Ignore
            || authoritative_stereo(molecules, &retry.coords).n_violations() == 0)
    {
        result = Ok(retry);
        minimizer = Mmff94Minimizer::Lbfgs;
    }
    let result = match result {
        Ok(result) => result,
        Err(error) => {
            run.record_stage(
                PipelineStage::ForceFieldMinimization,
                started.elapsed().as_millis() as u64,
            );
            return Err(run.fail(
                PipelineV2FailureCause::ForceField(error),
                PipelineStage::ForceFieldMinimization,
            ));
        }
    };
    run.evidence.force_field = Some(result.clone());
    run.evidence.last_known_coords = Some(result.coords.clone());
    run.complete_stage(PipelineStage::ForceFieldMinimization, started)?;
    Ok(ForceFieldStage {
        result,
        minimizer,
        reconcile_stereo,
    })
}

#[allow(clippy::result_large_err)]
fn finalize_stereo(
    molecules: &MoleculeContext<'_>,
    config: &PipelineV2Config,
    force_field: &ForceFieldStage,
    run: &mut PipelineRun,
) -> Result<FinalStereoStage, PipelineV2Failure> {
    let started = Instant::now();
    let mut verification = authoritative_stereo(molecules, &force_field.result.coords);
    run.record_stage(
        PipelineStage::FinalStereoVerify,
        started.elapsed().as_millis() as u64,
    );
    run.evidence.final_stereo = Some(verification.clone());
    let mut coords = force_field.result.coords.clone();
    let mut repair = None;
    if config.stereo_policy != StereoPolicy::Ignore && verification.n_violations() > 0 {
        let repair_started = Instant::now();
        if config.stereo_policy == StereoPolicy::RepairAndVerify
            && let Some(repaired) = try_post_minimization_repair(
                molecules,
                config,
                force_field,
                &coords,
            )
        {
            coords = repaired.0;
            repair = Some(repaired.1);
            verification = authoritative_stereo(molecules, &coords);
        }
        run.timings.post_min_stereo_repair_ms = repair_started.elapsed().as_millis() as u64;
        if repair.is_none() || verification.n_violations() > 0 {
            run.evidence.final_stereo = Some(verification);
            return Err(run.fail(
                PipelineV2FailureCause::FinalStereoViolation,
                PipelineStage::FinalStereoVerify,
            ));
        }
        run.evidence.final_stereo = Some(verification.clone());
        run.evidence.post_minimization_stereo_repair = repair.clone();
        run.evidence.last_known_coords = Some(coords.clone());
    }
    if config.stereo_policy != StereoPolicy::Ignore
        && config.fail_on_unevaluable_stereo
        && verification.n_unevaluable() > 0
    {
        return Err(run.fail(
            PipelineV2FailureCause::StereoUnevaluableUnderStrictPolicy,
            PipelineStage::FinalStereoVerify,
        ));
    }
    run.check_timeout(PipelineStage::FinalStereoVerify)?;
    Ok(FinalStereoStage {
        coords,
        verification,
        repair,
    })
}

fn try_post_minimization_repair(
    molecules: &MoleculeContext<'_>,
    config: &PipelineV2Config,
    force_field: &ForceFieldStage,
    coords: &Coords3D,
) -> Option<(Coords3D, StereoRepairSummary)> {
    let outcome = repair_stereo(molecules.working, coords).ok()?;
    let verified = verify_stereo(molecules.working, &outcome.coords);
    let sound = outcome.coords.is_finite()
        && worst_bond_length(molecules.working, &outcome.coords) <= MAX_SANE_BOND_LENGTH;
    if verified.n_violations() > 0 || !sound {
        return None;
    }
    let mut repaired_coords = outcome.coords;
    if molecules.expanded {
        let minimize_config = MinimizeConfig {
            max_steps: config.force_field_max_iterations,
            ..MinimizeConfig::default()
        };
        if let Ok(relaxed) = minimize_force_field_once(
            molecules,
            config,
            &minimize_config,
            force_field.reconcile_stereo,
            repaired_coords.clone(),
            force_field.minimizer,
        ) && verify_stereo(molecules.working, &relaxed.coords).n_violations() == 0
        {
            repaired_coords = relaxed.coords;
        }
    }
    Some((
        repaired_coords,
        StereoRepairSummary {
            repaired: outcome.repaired,
            failures: Vec::new(),
        },
    ))
}

#[allow(clippy::result_large_err)]
fn validate_final_geometry(
    mol: &Molecule,
    knowledge: &TorsionKnowledgeReport,
    torsion_optimization: Option<&TorsionOptimizationReport>,
    final_stereo: &FinalStereoStage,
    run: &mut PipelineRun,
) -> Result<FinalGeometryValidation, PipelineV2Failure> {
    let started = Instant::now();
    let torsion_energy = match evaluate_torsion_energy(
        mol,
        &final_stereo.coords,
        &knowledge.potentials,
    ) {
        Ok(report) => report.total_energy,
        Err(_) => {
            run.record_stage(
                PipelineStage::FinalGeometryValidationStage,
                started.elapsed().as_millis() as u64,
            );
            return Err(run.fail(
                PipelineV2FailureCause::FinalGeometryInvalid,
                PipelineStage::FinalGeometryValidationStage,
            ));
        }
    };
    let ring_closure_delta = torsion_optimization
        .map(|report| report.max_ring_closure_delta)
        .unwrap_or(0.0);
    let validation = compute_final_validation(
        mol,
        &final_stereo.coords,
        &final_stereo.verification,
        torsion_energy,
        ring_closure_delta,
    );
    run.evidence.final_validation = Some(validation.clone());
    run.record_stage(
        PipelineStage::FinalGeometryValidationStage,
        started.elapsed().as_millis() as u64,
    );
    if !validation.sound {
        return Err(run.fail(
            PipelineV2FailureCause::FinalGeometryInvalid,
            PipelineStage::FinalGeometryValidationStage,
        ));
    }
    run.check_timeout(PipelineStage::FinalGeometryValidationStage)?;
    Ok(validation)
}
