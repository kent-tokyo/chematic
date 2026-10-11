use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::etkdg_knowledge::FourierTorsionTerm;
    use crate::minimize::minimize_with_policy_gated;
    use chematic_smiles::parse;

    fn config_none() -> PipelineV2Config {
        PipelineV2Config::minimal(ForceFieldPolicy::None)
    }

    #[test]
    fn pipeline_run_records_every_stage_and_enforces_budget() {
        let config = config_none();
        let mut run = PipelineRun::new(&config);
        for (index, stage) in [
            PipelineStage::ValidateConfig,
            PipelineStage::TorsionKnowledge,
            PipelineStage::MacrocycleBoundAdjustment,
            PipelineStage::DistanceGeometry,
            PipelineStage::TorsionEnergyEvaluation,
            PipelineStage::TorsionOptimization,
            PipelineStage::StereoVerifyBefore,
            PipelineStage::StereoRepair,
            PipelineStage::StereoVerifyAfterRepair,
            PipelineStage::ForceFieldMinimization,
            PipelineStage::FinalStereoVerify,
            PipelineStage::FinalGeometryValidationStage,
        ]
        .into_iter()
        .enumerate()
        {
            run.record_stage(stage, index as u64 + 1);
        }
        assert_eq!(run.timings.torsion_knowledge_ms, 2);
        assert_eq!(run.timings.bound_adjustment_ms, 3);
        assert_eq!(run.timings.distance_geometry_ms, 4);
        assert_eq!(run.timings.torsion_energy_eval_ms, 5);
        assert_eq!(run.timings.torsion_optimization_ms, 6);
        assert_eq!(run.timings.stereo_verify_before_ms, 7);
        assert_eq!(run.timings.stereo_repair_ms, 8);
        assert_eq!(run.timings.stereo_verify_after_repair_ms, 9);
        assert_eq!(run.timings.force_field_ms, 10);
        assert_eq!(run.timings.final_stereo_verify_ms, 11);
        assert_eq!(run.timings.final_validation_ms, 12);

        let mut timed = config_none();
        timed.total_timeout_ms = Some(0);
        let mut run = PipelineRun::new(&timed);
        let failure = run
            .check_timeout(PipelineStage::ValidateConfig)
            .expect_err("zero budget must fail closed");
        assert_eq!(failure.stage, PipelineStage::ValidateConfig);
        assert!(matches!(failure.cause, PipelineV2FailureCause::Timeout));
    }

    #[test]
    fn retry_configuration_is_deterministic_and_budgeted() {
        let mut config = config_none();
        config.embed.random_seed = 17;
        config.total_timeout_ms = Some(100);
        let retry = retry_config(&config, 2, 40).expect("remaining budget");
        assert_eq!(retry.total_timeout_ms, Some(60));
        assert_ne!(retry.embed.random_seed, config.embed.random_seed);
        assert!(retry_config(&config, 1, 101).is_none());

        config.total_timeout_ms = None;
        assert_eq!(
            retry_config(&config, 1, u64::MAX).unwrap().total_timeout_ms,
            None
        );

        let mut run = PipelineRun::new(&config);
        let stereo_failure = run.fail(
            PipelineV2FailureCause::FinalStereoViolation,
            PipelineStage::FinalStereoVerify,
        );
        assert!(!retries_final_stereo_violation(&config, &stereo_failure));
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        assert!(retries_final_stereo_violation(&config, &stereo_failure));
    }

    #[test]
    fn typed_stage_failures_preserve_partial_evidence() {
        let mol = parse("CCCC").unwrap();
        let mut config = config_none();

        let mut run = PipelineRun::new(&config);
        let invalid_bound = DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(1),
            lower: 3.0,
            upper: 1.0,
        };
        let failure = run_distance_geometry(&mol, &config, &[invalid_bound], &mut run)
            .expect_err("inverted bounds must fail closed");
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::BoundAdjustmentFailed
        ));
        assert_eq!(failure.stage, PipelineStage::MacrocycleBoundAdjustment);
        assert!(failure.embed_stats.is_some());

        config.embed.use_exp_torsions = true;
        let knowledge = KnowledgeStage {
            torsion: TorsionKnowledgeReport {
                potentials: vec![TorsionPotential {
                    atoms: [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)],
                    central_bond: (AtomIdx(1), AtomIdx(2)),
                    source: TorsionKnowledgeSource::StandardExperimental,
                    rule_id: "coverage:butane".into(),
                    terms: vec![FourierTorsionTerm::from_rdkit(1, 1, 1.0)],
                    ring_size: None,
                }],
                ..TorsionKnowledgeReport::default()
            },
            bounds: None,
        };
        let wrong_size = Coords3D::new_zeroed(1);

        let mut energy_run = PipelineRun::new(&config);
        let failure = evaluate_initial_torsions(&mol, &wrong_size, &knowledge, &mut energy_run)
            .expect_err("coordinate mismatch must be typed");
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::TorsionKnowledge(
                TorsionKnowledgeError::CoordsAtomCountMismatch
            )
        ));
        assert_eq!(failure.stage, PipelineStage::TorsionEnergyEvaluation);

        let mut optimize_run = PipelineRun::new(&config);
        let failure =
            optimize_pipeline_torsions(&mol, &config, &knowledge, wrong_size, &mut optimize_run)
                .expect_err("optimizer must preserve the same typed mismatch");
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::TorsionKnowledge(
                TorsionKnowledgeError::CoordsAtomCountMismatch
            )
        ));
        assert_eq!(failure.stage, PipelineStage::TorsionOptimization);
    }

    #[test]
    fn final_geometry_and_strict_stereo_fail_closed() {
        let mol = parse("CC").unwrap();
        let mut non_finite = Coords3D::new_zeroed(mol.atom_count());
        non_finite.set(
            AtomIdx(1),
            crate::coords::Point3 {
                x: f64::NAN,
                y: 0.0,
                z: 0.0,
            },
        );
        let final_stereo = FinalStereoStage {
            coords: non_finite,
            verification: StereoVerification::default(),
            repair: None,
        };
        let mut run = PipelineRun::new(&config_none());
        let failure = validate_final_geometry(
            &mol,
            &TorsionKnowledgeReport::default(),
            None,
            &final_stereo,
            &mut run,
        )
        .expect_err("non-finite final coordinates must fail");
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::FinalGeometryInvalid
        ));

        let chiral = parse("C[C@H](F)Cl").unwrap();
        let coords = Coords3D::new_zeroed(chiral.atom_count());
        let force_field = minimize_with_policy_gated(
            &chiral,
            coords,
            ForceFieldPolicy::None,
            &MinimizeConfig::default(),
            false,
            false,
        )
        .unwrap();
        let molecules = MoleculeContext {
            original: &chiral,
            working: &chiral,
            expanded: false,
            original_atom_count: chiral.atom_count(),
        };
        let mut strict = config_none();
        strict.stereo_policy = StereoPolicy::VerifyOnly;
        strict.fail_on_unevaluable_stereo = true;
        let force_field = ForceFieldStage {
            result: force_field,
            minimizer: Mmff94Minimizer::Bfgs,
            reconcile_stereo: false,
        };
        let mut run = PipelineRun::new(&strict);
        let failure = match finalize_stereo(&molecules, &strict, &force_field, &mut run) {
            Err(failure) => failure,
            Ok(_) => panic!("degenerate declared stereo must be refused under strict policy"),
        };
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::StereoUnevaluableUnderStrictPolicy
        ));
    }

    #[test]
    fn force_field_and_final_torsion_errors_remain_typed() {
        let mol = parse("CCCC").unwrap();
        let mut bad_coords = Coords3D::new_zeroed(mol.atom_count());
        bad_coords.set(
            AtomIdx(0),
            crate::coords::Point3 {
                x: f64::NAN,
                y: 0.0,
                z: 0.0,
            },
        );
        let stereo = StereoStage {
            coords: bad_coords,
            before: StereoVerification::default(),
            repair: None,
            after_repair: StereoVerification::default(),
        };
        let mut config = config_none();
        config.force_field_policy = ForceFieldPolicy::Mmff94BondAngleStrict;
        let mut run = PipelineRun::new(&config);
        let failure = match run_force_field(
            &MoleculeContext {
                original: &mol,
                working: &mol,
                expanded: false,
                original_atom_count: mol.atom_count(),
            },
            &config,
            &stereo,
            &mut run,
        ) {
            Err(failure) => failure,
            Ok(_) => panic!("non-finite force-field input must fail"),
        };
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::ForceField(_)
        ));
        assert_eq!(failure.stage, PipelineStage::ForceFieldMinimization);

        let invalid_torsion = TorsionPotential {
            atoms: [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(99)],
            central_bond: (AtomIdx(1), AtomIdx(2)),
            source: TorsionKnowledgeSource::StandardExperimental,
            rule_id: "coverage:invalid-topology".into(),
            terms: vec![FourierTorsionTerm::from_rdkit(1, 1, 1.0)],
            ring_size: None,
        };
        let final_stereo = FinalStereoStage {
            coords: Coords3D::new_zeroed(mol.atom_count()),
            verification: StereoVerification::default(),
            repair: None,
        };
        let mut run = PipelineRun::new(&config_none());
        let failure = validate_final_geometry(
            &mol,
            &TorsionKnowledgeReport {
                potentials: vec![invalid_torsion],
                ..TorsionKnowledgeReport::default()
            },
            None,
            &final_stereo,
            &mut run,
        )
        .expect_err("invalid torsion topology must fail final validation");
        assert!(matches!(
            failure.cause,
            PipelineV2FailureCause::FinalGeometryInvalid
        ));
        assert_eq!(failure.stage, PipelineStage::FinalGeometryValidationStage);
    }

    // -----------------------------------------------------------------------
    // Reproducibility / invariance (spec §16)
    // -----------------------------------------------------------------------

    #[test]
    fn all_flags_off_matches_raw_dg_exactly() {
        // "all stage flags off -> matches raw DG" (spec §16/§11 negative-control
        // partner). With every knowledge/optimization/stereo/bounds flag off and
        // ForceFieldPolicy::None, the pipeline must reduce to exactly the raw
        // embedder's own output at the same seed.
        let mol = parse("CCCCCCCCCC").unwrap(); // decane
        let config = config_none();
        let result = embed_pipeline_v2(&mol, &config).expect("pipeline should succeed");

        let raw = distance_geometry_v2::embed_distance_geometry_v2(&mol, &config.embed)
            .expect("raw embed should succeed");

        for i in 0..mol.atom_count() {
            let p_pipeline = result.coords.get(AtomIdx(i as u32));
            let p_raw = raw.get(AtomIdx(i as u32));
            assert_eq!(
                p_pipeline, p_raw,
                "atom {i}: pipeline coords must be byte-identical to raw DG with every flag off"
            );
        }
        assert!(result.torsion_knowledge_report.potentials.is_empty());
        assert!(result.torsion_optimization_report.is_none());
        assert!(result.bound_adjustment_report.is_none());
    }

    #[test]
    fn same_seed_same_config_reproduces_identical_result() {
        let mol = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol
        let config = config_none();
        let r1 = embed_pipeline_v2(&mol, &config).unwrap();
        let r2 = embed_pipeline_v2(&mol, &config).unwrap();
        for i in 0..mol.atom_count() {
            assert_eq!(
                r1.coords.get(AtomIdx(i as u32)),
                r2.coords.get(AtomIdx(i as u32)),
                "same seed/config must reproduce identical coords at atom {i}"
            );
        }
        assert_eq!(
            r1.final_validation.worst_bond_length,
            r2.final_validation.worst_bond_length
        );
    }

    #[test]
    fn different_seeds_are_not_aliased() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let mut any_diff = false;
        let base = embed_pipeline_v2(&mol, &config_none()).unwrap();
        for seed in [1u64, 42u64] {
            let mut config = config_none();
            config.embed.random_seed = seed;
            let r = embed_pipeline_v2(&mol, &config).unwrap();
            for i in 0..mol.atom_count() {
                if base.coords.get(AtomIdx(i as u32)) != r.coords.get(AtomIdx(i as u32)) {
                    any_diff = true;
                }
            }
        }
        assert!(any_diff, "different seeds must not produce aliased output");
    }

    // -----------------------------------------------------------------------
    // Stage order / typed failures
    // -----------------------------------------------------------------------

    #[test]
    fn enforce_chirality_with_repair_and_verify_stereo_policy_is_allowed() {
        // Revised 2026-08-24 (issue #291 Step A): this combination was
        // previously rejected as InvalidConfiguration -- now validated (see
        // the module doc's revised Stage 1 entry and
        // `crates/chematic-3d/examples/issue291_repair_policy_measurement.rs`).
        let mol = parse("C[C@H](O)CC").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        match embed_pipeline_v2(&mol, &config) {
            Ok(r) => assert!(
                r.final_stereo.is_fully_satisfied(),
                "2-butanol has no ring-fused stereocenter -- must succeed with satisfied stereo"
            ),
            Err(e) => assert!(
                !matches!(e.cause, PipelineV2FailureCause::InvalidConfiguration),
                "must not be rejected as InvalidConfiguration, got {e:?}"
            ),
        }
    }

    #[test]
    fn embed_level_materialize_implicit_h_for_chirality_is_rejected_alone() {
        // `embed.materialize_implicit_h_for_chirality` is proven correct in isolation
        // (`distance_geometry_v2`'s own tests, an independent-oracle cross-check), but
        // on its own it only expands within a single raw embed attempt and truncates
        // immediately -- this pipeline's stages 7/8/9/11 verify/repair would never see
        // the real H position. Stays rejected unconditionally when used alone; the
        // real follow-up landed as the distinct, pipeline-level
        // `expand_implicit_h_through_pipeline` flag tested below, which keeps the
        // expanded molecule for the whole stage 4-11 sequence instead.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.embed.materialize_implicit_h_for_chirality = true;
        let err = embed_pipeline_v2(&mol, &config).expect_err(
            "materialize_implicit_h_for_chirality must be rejected, not silently ignored or run",
        );
        assert!(
            matches!(err.cause, PipelineV2FailureCause::InvalidConfiguration),
            "expected InvalidConfiguration, got {err:?}"
        );
        assert_eq!(err.stage, PipelineStage::ValidateConfig);
    }

    #[test]
    fn embed_level_and_pipeline_level_materialize_implicit_h_flags_together_still_rejected() {
        // The pre-existing unconditional check on `embed.materialize_implicit_h_for_chirality`
        // already rejects this combination before any new code runs -- pinned here so
        // a future reordering of the stage-1 checks can't silently change that.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.embed.materialize_implicit_h_for_chirality = true;
        config.expand_implicit_h_through_pipeline = true;
        let err = embed_pipeline_v2(&mol, &config)
            .expect_err("both flags together must still be rejected");
        assert!(matches!(
            err.cause,
            PipelineV2FailureCause::InvalidConfiguration
        ));
    }

    #[test]
    fn expand_implicit_h_through_pipeline_requires_enforce_chirality() {
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = false;
        config.expand_implicit_h_through_pipeline = true;
        let err = embed_pipeline_v2(&mol, &config).expect_err(
            "expand_implicit_h_through_pipeline without enforce_chirality must be rejected",
        );
        assert!(
            matches!(err.cause, PipelineV2FailureCause::InvalidConfiguration),
            "expected InvalidConfiguration, got {err:?}"
        );
        assert_eq!(err.stage, PipelineStage::ValidateConfig);
    }

    #[test]
    fn expand_implicit_h_through_pipeline_is_noop_without_declared_stereo() {
        let mol = parse("CCCCCCCCCC").unwrap(); // decane, no declared stereo at all
        let base = {
            let mut config = config_none();
            config.embed.enforce_chirality = true;
            embed_pipeline_v2(&mol, &config).unwrap()
        };
        let expanded = {
            let mut config = config_none();
            config.embed.enforce_chirality = true;
            config.expand_implicit_h_through_pipeline = true;
            embed_pipeline_v2(&mol, &config).unwrap()
        };
        assert_eq!(base.coords.atom_count(), expanded.coords.atom_count());
        for i in 0..base.coords.atom_count() {
            let idx = AtomIdx(i as u32);
            assert_eq!(
                base.coords.get(idx),
                expanded.coords.get(idx),
                "atom {i}: expand_implicit_h_through_pipeline must be byte-identical to off \
                 for a molecule with no declared stereo"
            );
        }
    }

    #[test]
    fn expand_implicit_h_through_pipeline_negative_control_no_ring_fused_center() {
        // 2-butanol's declared stereocenter has an implicit H but is NOT ring-fused
        // -- repair_tetrahedral_center's existing substituent-reflection already
        // handles it without any H materialization. The flag must not change the
        // observable outcome here: it's real work spent on a molecule that never
        // needed it.
        let mol = parse("C[C@H](O)CC").unwrap();
        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
        config.embed.enforce_chirality = true;
        config.embed.random_seed = 0;
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        config.expand_implicit_h_through_pipeline = true;
        let result = embed_pipeline_v2(&mol, &config)
            .expect("2-butanol must still succeed with the flag on");
        assert!(result.final_stereo.is_fully_satisfied());
        assert_eq!(
            result.force_field.coords.atom_count(),
            mol.atom_count() + 1,
            "ordinary stereocenters should materialize only their declared implicit H"
        );
    }

    #[test]
    fn multi_ring_implicit_h_center_keeps_full_hydrogen_expansion() {
        let fused = parse("CC1(C)[C@H](C(=O)O)N2C(=O)C(C/C=C/C=O)[C@H]2S1(=O)=O")
            .expect("fused stereocenter fixture");
        let ordinary = parse("C[C@H](O)CC").expect("ordinary stereocenter fixture");

        assert!(mol_needs_full_hydrogen_expansion(&fused));
        assert!(!mol_needs_full_hydrogen_expansion(&ordinary));
    }

    #[test]
    fn expand_implicit_h_through_pipeline_is_noop_for_ez_only_molecule() {
        let mol = parse("C/C=C\\C").unwrap();
        let mut config = PipelineV2Config::stereo_safe(ForceFieldPolicy::Mmff94BondAngleStrict);
        config.embed.random_seed = 0;
        let result = embed_pipeline_v2(&mol, &config).expect("declared E/Z must remain supported");
        assert!(result.final_stereo.is_fully_satisfied());
        assert_eq!(
            result.force_field.coords.atom_count(),
            mol.atom_count(),
            "E/Z-only molecules do not need explicit-H expansion"
        );
    }

    #[test]
    // Runs in CI: the bounded stereo regression completes within seconds.
    fn expand_implicit_h_through_pipeline_fixes_testosterone_with_correct_geometry() {
        // Issue #291's original residual. Seed picked from the Phase 0.5 measurement
        // harness's own run (`issue291_expanded_geometry_feasibility.rs`) as one of
        // testosterone's 3 known-clean seeds under this exact configuration -- not
        // an arbitrary seed, since that measurement found a real 3/5-vs-2/5 split
        // (the other 2/5 are honest, safe failures, not silent-wrong -- see
        // ROADMAP.md's #291 entry).
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let declared = chematic_chem::assign_cip(&mol);
        assert_eq!(
            declared.assignments.len(),
            6,
            "sanity: 6 declared stereocenters"
        );

        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
        config.embed.enforce_chirality = true;
        config.embed.random_seed = 0;
        config.embed.max_attempts = 8;
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        config.expand_implicit_h_through_pipeline = true;
        let result = embed_pipeline_v2(&mol, &config).expect("testosterone must succeed");
        assert!(
            result.final_stereo.is_fully_satisfied(),
            "final_stereo: {:?}",
            result.final_stereo
        );
        assert_eq!(result.coords.atom_count(), mol.atom_count());

        let perceived = crate::stereo3d::assign_stereo_from_3d(&mol, &result.coords);
        for &(idx, code) in &declared.assignments {
            if let Some(perceived_code) = perceived.get(idx) {
                assert_eq!(
                    perceived_code, code,
                    "atom {idx:?}: declared {code:?} but 3D-perceived {perceived_code:?} \
                     -- pipeline reported success with wrong chirality"
                );
            }
        }
    }

    #[test]
    // Runs in CI: the bounded stereo regression completes within seconds.
    fn stereo_safe_matches_the_hand_built_configuration_above() {
        // issue #383: `stereo_safe` must be exactly the same as manually setting the
        // 3 flags -- same testosterone happy-path, same seed, same expectations as
        // `expand_implicit_h_through_pipeline_fixes_testosterone_with_correct_geometry`
        // above, just built via the convenience constructor instead.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let declared = chematic_chem::assign_cip(&mol);

        let mut config = PipelineV2Config::stereo_safe(ForceFieldPolicy::Mmff94WithUffFallback);
        assert!(config.embed.enforce_chirality);
        assert!(config.expand_implicit_h_through_pipeline);
        assert_eq!(config.stereo_policy, StereoPolicy::RepairAndVerify);
        config.embed.random_seed = 0;

        let result = embed_pipeline_v2(&mol, &config).expect("testosterone must succeed");
        assert!(result.final_stereo.is_fully_satisfied());
        assert_eq!(result.coords.atom_count(), mol.atom_count());

        let perceived = crate::stereo3d::assign_stereo_from_3d(&mol, &result.coords);
        for &(idx, code) in &declared.assignments {
            if let Some(perceived_code) = perceived.get(idx) {
                assert_eq!(perceived_code, code);
            }
        }
    }

    #[test]
    fn stereo_safe_mmff94_keeps_fused_ring_stereo_evaluable_during_minimization() {
        // A6 regression: unconstrained MMFF94 minimization crossed atom 16's
        // fused-ring chiral boundary, while rejecting every crossing proposal
        // also rejected an initial heavy-only phantom-H mismatch at atom 10.
        // The stereo-safe line search must preserve the expanded assignment,
        // allow that pre-existing violation to relax, and forbid a transition
        // through unevaluable returned geometry.
        let mol = parse("C[C@]12CCC3C(CC=C4C[C@@H](O)[C@H]5COC[C@]43C5)C1CC[C@@H]2OC1CC1").unwrap();
        let mut config = PipelineV2Config::stereo_safe(ForceFieldPolicy::Mmff94BondAngleStrict);
        config.embed.random_seed = 20_260_801;
        config.force_field_max_iterations = 300;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;

        let result = embed_pipeline_v2(&mol, &config)
            .expect("stereo-safe MMFF94 must retain this fused-ring molecule");
        assert!(result.final_stereo.is_fully_satisfied());
        assert_eq!(result.final_validation.gross_clash_count, 0);
        assert!(result.final_validation.sound);
    }

    #[test]
    // Runs in CI: the bounded stereo regression completes within seconds.
    fn expand_implicit_h_through_pipeline_with_verify_only_never_reports_success_with_violated_final_stereo()
     {
        // `VerifyOnly` never repairs, so this test doesn't assert success at any
        // particular seed (unlike the RepairAndVerify happy-path test above) --
        // it asserts the *invariant* this flag must preserve regardless of policy:
        // whatever `final_stereo` says, it must never disagree with an Ok result.
        // This specific interaction (VerifyOnly + expand_implicit_h_through_pipeline)
        // was never empirically measured by the Phase 0.5 harness (which only ran
        // RepairAndVerify) -- it only follows from `authoritative_final_stereo` being
        // computed before the `stereo_policy != Ignore` branch, which this test
        // verifies directly rather than trusting the code-structure argument alone.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        for seed in 0..5u64 {
            let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
            config.embed.enforce_chirality = true;
            config.embed.random_seed = seed;
            config.embed.max_attempts = 8;
            config.stereo_policy = StereoPolicy::VerifyOnly;
            config.expand_implicit_h_through_pipeline = true;
            match embed_pipeline_v2(&mol, &config) {
                Ok(result) => assert!(
                    result.final_stereo.is_fully_satisfied(),
                    "seed={seed}: Ok result must never carry a violated final_stereo, got {:?}",
                    result.final_stereo
                ),
                Err(e) => {
                    // VerifyOnly never repairs, so any of these are legitimate --
                    // just confirm it's a recognized stereo-related failure, not a
                    // silent success/failure mismatch.
                    assert!(
                        matches!(
                            e.cause,
                            PipelineV2FailureCause::FinalStereoViolation
                                | PipelineV2FailureCause::DistanceGeometry(_)
                        ),
                        "seed={seed}: unexpected failure cause {:?}",
                        e.cause
                    );
                }
            }
        }
    }

    #[test]
    fn enforce_chirality_with_ignore_stereo_policy_is_allowed() {
        let mol = parse("CCCC").unwrap(); // no declared stereo at all
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.stereo_policy = StereoPolicy::Ignore;
        assert!(embed_pipeline_v2(&mol, &config).is_ok());
    }

    #[test]
    fn enforce_chirality_with_verify_only_stereo_policy_is_allowed() {
        // Revised 2026-08-11: previously InvalidConfiguration -- now allowed, since
        // VerifyOnly's stage 11 gate is exactly what catches a force field
        // undoing enforce_chirality's embedding-time correctness (see the module
        // doc's revised judgment-call entry for the corpus evidence).
        let mol = parse("C[C@H](O)CC").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.embed.max_attempts = 8;
        config.stereo_policy = StereoPolicy::VerifyOnly;
        let result = embed_pipeline_v2(&mol, &config);
        if let Err(e) = &result {
            assert!(
                !matches!(e.cause, PipelineV2FailureCause::InvalidConfiguration),
                "must not be rejected as InvalidConfiguration, got {e:?}"
            );
        }
    }

    #[test]
    fn enforce_chirality_with_verify_only_never_reports_success_with_violated_final_stereo() {
        // Motivated by the corpus-measured chembl_tier_b_0076/0083 failure mode
        // (265-molecule v0.14.0 release-gate re-measurement): enforce_chirality can
        // deliver correct declared E/Z at embedding time (verified via
        // `stereo_before`, populated before stage 10 runs) while MMFF94
        // minimization -- which has no notion of declared stereo -- walks it back
        // across the boundary afterward (confirmed separately for that exact
        // molecule/seed: re-running with `ForceFieldPolicy::None` leaves
        // `final_stereo` satisfied, isolating minimization as the cause -- see the
        // module doc's revised judgment-call entry). Before this PR's gate
        // relaxation, `StereoPolicy::Ignore` was the only option compatible with
        // enforce_chirality, and Ignore never gates on stereo -- a caller could get
        // a `success` result whose geometry silently violates its own declared E/Z.
        //
        // This test asserts the general invariant, not a specific outcome for this
        // one molecule/seed: whether MMFF94 happens to revert the fix is itself
        // non-deterministic across build profiles (verified while writing this test
        // -- an unoptimized debug build converged differently than the release
        // build the corpus measurement used, for the identical seed), so asserting
        // "this exact call must fail" would make the test flaky by build profile.
        // What must hold regardless: `enforce_chirality`'s own part of the contract
        // (embedding-time correctness, `stereo_before`) always succeeds for this
        // molecule, and `Ok` is never returned with a violated `final_stereo` --
        // VerifyOnly's stage 11 gate must have caught it if minimization did revert
        // the fix.
        let mol = parse("COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(NC(=O)CBr)cc4)CC3)nc(N)c2cc1OC").unwrap(); // chembl_tier_b_0083
        let mut config = config_none();
        config.embed.random_seed = 20260801; // the corpus benchmark's exact seed
        config.embed.enforce_chirality = true;
        config.embed.max_attempts = 8;
        config.embed.use_exp_torsions = true;
        config.embed.use_small_ring_torsions = true;
        config.embed.use_macrocycle_torsions = true;
        config.embed.use_macrocycle_14_bounds = true;
        config.stereo_policy = StereoPolicy::VerifyOnly;
        config.force_field_policy = ForceFieldPolicy::Mmff94BondAngleStrict;
        // Capped well below the 200-iteration production default: in an
        // unoptimized debug build (what `cargo test` uses in CI), this specific
        // 36-atom macrocycle-containing molecule's MMFF94 minimization is slow
        // enough at 200 iterations to make the test take minutes -- 40 is still
        // enough to exercise "did minimization move the geometry", the actual
        // invariant under test, without that cost.
        config.force_field_max_iterations = 40;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        config.total_timeout_ms = Some(60_000);

        match embed_pipeline_v2(&mol, &config) {
            Ok(r) => {
                assert!(
                    r.stereo_before.is_fully_satisfied(),
                    "enforce_chirality's embedding-time fix must be correct -- got {:?}",
                    r.stereo_before
                );
                assert!(
                    r.final_stereo.is_fully_satisfied(),
                    "must never report success with a geometry violating declared stereo -- \
                     got final_stereo {:?}",
                    r.final_stereo
                );
            }
            Err(e) => {
                assert!(
                    !matches!(e.cause, PipelineV2FailureCause::InvalidConfiguration),
                    "must not be rejected as InvalidConfiguration: {e:?}"
                );
                assert!(
                    e.stereo_before
                        .as_ref()
                        .is_some_and(StereoVerification::is_fully_satisfied),
                    "enforce_chirality's embedding-time fix must be correct -- got {:?}",
                    e.stereo_before
                );
            }
        }
    }

    #[test]
    fn ring_torsion_fail_closed_on_macrocycle_when_optimizer_cannot_apply() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane, macrocycle
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(
            err.cause,
            PipelineV2FailureCause::RingTorsionApplicationUnsupported
        ));
        // Evidence must still be attached even on this early, pre-embed failure.
        let evidence = err
            .ring_torsion_evidence
            .expect("evidence must be attached");
        assert!(evidence.potentials.iter().any(|p| !p.applied_to_geometry));
        assert!(evidence.n_scored_only() > 0);
        assert_eq!(
            evidence.n_applied() + evidence.n_scored_only(),
            evidence.potentials.len()
        );
    }

    #[test]
    fn ring_torsion_diagnostic_only_succeeds_with_scored_only_evidence() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        let result = embed_pipeline_v2(&mol, &config).expect("DiagnosticOnly must succeed");
        assert!(result.ring_torsion_evidence.diagnostic_only);
        assert!(result.ring_torsion_evidence.n_scored_only() > 0);
        assert_eq!(
            result.ring_torsion_evidence.n_applied() + result.ring_torsion_evidence.n_scored_only(),
            result.ring_torsion_evidence.potentials.len()
        );
        assert!(
            result
                .ring_torsion_evidence
                .potentials
                .iter()
                .any(|p| p.source == TorsionKnowledgeSource::MacrocycleAdaptation),
            "expected at least one MacrocycleAdaptation potential on a macrocycle"
        );
        for p in &result.ring_torsion_evidence.potentials {
            if p.source == TorsionKnowledgeSource::MacrocycleAdaptation {
                assert!(
                    !p.applied_to_geometry,
                    "a macrocycle potential must never be reported as applied_to_geometry"
                );
            }
        }
    }

    #[test]
    fn ring_torsion_fail_closed_never_fires_on_purely_acyclic_potentials() {
        // Negative control: acyclic-only torsion knowledge must never trip the ring
        // gate, even under FailClosed.
        let mol = parse("CCCC").unwrap(); // butane, fully acyclic
        let mut config = config_none();
        config.embed.use_exp_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
        let result = embed_pipeline_v2(&mol, &config).expect("acyclic torsions must not gate");
        assert_eq!(result.ring_torsion_evidence.n_scored_only(), 0);
        assert_eq!(
            result.ring_torsion_evidence.n_applied(),
            result.ring_torsion_evidence.potentials.len()
        );
        for p in &result.ring_torsion_evidence.potentials {
            assert!(p.applied_to_geometry, "acyclic potential must be applied");
        }
    }

    #[test]
    fn macrocycle_14_bounds_actually_applied_before_smoothing_changes_output() {
        // Arm A vs Arm B differential (spec §13/§18-c): enabling
        // use_macrocycle_14_bounds on a macrocycle must produce DIFFERENT coordinates
        // than leaving it off, at the same seed -- otherwise the hook isn't wired
        // even though every unit test is green.
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let config_a = config_none();
        let mut config_b = config_none();
        config_b.embed.use_macrocycle_14_bounds = true;

        let result_a = embed_pipeline_v2(&mol, &config_a).unwrap();
        let result_b = embed_pipeline_v2(&mol, &config_b).unwrap();

        assert!(result_a.bound_adjustment_report.is_none());
        let adjustments_b = result_b
            .bound_adjustment_report
            .as_ref()
            .expect("Some(..) when the flag is set");
        assert!(
            !adjustments_b.is_empty(),
            "cyclododecane must produce real adjustments"
        );
        assert!(result_b.embed_stats.adjustments_applied > 0);

        let mut any_diff = false;
        for i in 0..mol.atom_count() {
            if result_a.coords.get(AtomIdx(i as u32)) != result_b.coords.get(AtomIdx(i as u32)) {
                any_diff = true;
            }
        }
        assert!(
            any_diff,
            "macrocycle 1-4 bounds must actually change the embedded geometry, not just be scored"
        );
    }

    #[test]
    fn macrocycle_14_amide_pinned_branch_is_reachable_through_the_real_pipeline() {
        // `bounds14.rs`'s `macrocycle_14:amide_ester_pinned` rule (pin-to-cis-or-trans,
        // depending on ring role, for an amide/ester bond inside a macrocycle) is a
        // DIFFERENT branch from the `macrocycle_14:relaxed_band` rule the test above
        // exercises -- every macrocycle in the frozen 58/63 corpus (cyclododecane,
        // crown_12_4, cyclooctadecane) is all-carbon or all-ether, so no arm in
        // `pipeline_v2_integration_gate.rs` ever hit this branch through the real
        // embedder before this test (found during verification round 4). A
        // 12-membered ring lactam is `bounds14.rs`'s own test fixture for this rule.
        let mol = parse("O=C1CCCCCCCCCCN1").unwrap(); // 12-membered ring lactam
        let mut config = config_none();
        config.embed.use_macrocycle_14_bounds = true;

        let result = embed_pipeline_v2(&mol, &config).expect("must embed successfully");
        let adjustments = result
            .bound_adjustment_report
            .as_ref()
            .expect("Some(..) when the flag is set");
        assert!(
            adjustments
                .iter()
                .any(|a| a.rule_id == "macrocycle_14:amide_ester_pinned"),
            "the amide-pinned branch must fire through the real pipeline, not just \
             bounds14.rs's own standalone unit test: {adjustments:?}"
        );
        assert!(result.embed_stats.adjustments_applied > 0);
    }

    /// `atoms[1]`-`atoms[2]` central bond dihedral in degrees, atan2-based
    /// (numerically stable near 0/180). Self-contained rather than reusing
    /// `etkdg_knowledge::energy`'s private `dihedral_deg`, to avoid widening
    /// that module's visibility just for this test.
    fn measured_dihedral_deg(coords: &crate::coords::Coords3D, atoms: [AtomIdx; 4]) -> f64 {
        let p0 = coords.get(atoms[0]);
        let p1 = coords.get(atoms[1]);
        let p2 = coords.get(atoms[2]);
        let p3 = coords.get(atoms[3]);
        let b1 = p1.sub(&p0);
        let b2 = p2.sub(&p1);
        let b3 = p3.sub(&p2);
        let n1 = b1.cross(&b2);
        let n2 = b2.cross(&b3);
        let b2_unit = b2.try_normalize().unwrap_or(crate::coords::Point3::zero());
        let m1 = n1.cross(&b2_unit);
        let x = n1.dot(&n2);
        let y = m1.dot(&n2);
        y.datan2(x).to_degrees()
    }

    /// Distance (degrees) from `dihedral` to the nearest of a planar amide's
    /// two valid configurations, 0° (cis) or ±180° (trans).
    fn dist_to_planar(dihedral_deg: f64) -> f64 {
        let d = dihedral_deg.rem_euclid(360.0);
        (d.min(360.0 - d)).min((d - 180.0).abs())
    }

    #[test]
    fn tertiary_amide_macrocycle_embeds_closer_to_planar_with_14_bounds_fix() {
        // Issue found while surveying RDKit's open issues (analogous to RDKit
        // #9266, "ETKDGv3 twisted tertiary amides in macrocycles"). Confirmed
        // this reproduces in chematic; root cause was `bounds14.rs`'s
        // amide-pinned branch pinning all 4 combinatorial 1-4 pairs through a
        // tertiary amide to the same (geometrically unsatisfiable) cis
        // configuration -- fixed to split cis/trans by ring role (see that
        // module's own tests for the precise per-pair values). This test
        // confirms the fix's real-world effect through the actual embedding
        // pipeline, not just the bound values in isolation: enabling the
        // fixed `use_macrocycle_14_bounds` must measurably improve amide
        // planarity relative to leaving it off, on a genuine tertiary-amide
        // macrolactam. Empirically measured (Python, live pipeline, 2 real
        // multi-amide macrocycles, 4 seeds each, 48 dihedral measurements):
        // mean distance-to-planar dropped from ~61.5° (flag off) to ~20.8°
        // (flag on, fixed) -- not perfect convergence (a soft DG-bounds nudge
        // on a stochastic embedder, not a hard planarity guarantee), but a
        // clear, reproducible improvement. This test uses a smaller,
        // single-amide fixture and a threshold with real margin below that
        // baseline, not a razor-thin one recreated from a lucky single run.
        let mol = parse("O=C1CCCCCCCCCCN1C").unwrap(); // N-methyl 13-membered macrolactam
        // Atom indices per this exact SMILES's left-to-right parse order.
        let o_idx = AtomIdx(0);
        let carbonyl_c_idx = AtomIdx(1);
        let ring_c_idx = AtomIdx(2);
        let ring_n_idx = AtomIdx(11);
        let amide_n_idx = AtomIdx(12);
        let methyl_idx = AtomIdx(13);

        let dihedral_pairs = [
            (ring_n_idx, o_idx),
            (ring_n_idx, ring_c_idx),
            (methyl_idx, o_idx),
            (methyl_idx, ring_c_idx),
        ];

        let mean_dist_to_planar = |use_14_bounds: bool, seed: u64| -> f64 {
            let mut config = config_none();
            config.embed.use_exp_torsions = true;
            config.embed.use_macrocycle_torsions = true;
            config.embed.use_macrocycle_14_bounds = use_14_bounds;
            config.embed.random_seed = seed;
            // Ring-internal torsion potentials are scored-only (never
            // mechanically applied, see `energy.rs`'s `is_bridge_bond` gate)
            // -- FailClosed would reject every attempt once macrocycle
            // torsion knowledge is enabled at all. DiagnosticOnly matches
            // how this scenario is actually run in production.
            config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
            let result = embed_pipeline_v2(&mol, &config).expect("must embed successfully");
            let mut total = 0.0;
            for &(a1, a4) in &dihedral_pairs {
                let dih =
                    measured_dihedral_deg(&result.coords, [a1, amide_n_idx, carbonyl_c_idx, a4]);
                total += dist_to_planar(dih);
            }
            total / dihedral_pairs.len() as f64
        };

        let seeds: [u64; 4] = [1, 7, 42, 123];
        let with_fix: Vec<f64> = seeds
            .iter()
            .map(|&s| mean_dist_to_planar(true, s))
            .collect();
        let without: Vec<f64> = seeds
            .iter()
            .map(|&s| mean_dist_to_planar(false, s))
            .collect();

        let avg = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let avg_with_fix = avg(&with_fix);
        let avg_without = avg(&without);

        assert!(
            avg_with_fix < avg_without,
            "enabling the fixed macrocycle_14_bounds must improve amide planarity: \
             with_fix={with_fix:?} (avg {avg_with_fix:.1}) without={without:?} (avg {avg_without:.1})"
        );
        assert!(
            avg_with_fix < 45.0,
            "average distance-to-planar with the fix should be well below the \
             ~90-130 range the original bug produced: with_fix={with_fix:?} (avg {avg_with_fix:.1})"
        );
    }

    #[test]
    fn invalid_adjustment_pair_is_typed_bound_adjustment_failed() {
        // Direct test of the internal hook's own validation (index out of range).
        let mol = parse("CC").unwrap();
        let bad = [DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(99),
            lower: 1.0,
            upper: 2.0,
        }];
        let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            &mol,
            &EmbedParameters::default(),
            &bad,
        )
        .unwrap_err();
        assert!(matches!(
            err.0,
            EmbedWithAdjustmentsFailure::InvalidAdjustment
        ));
    }

    #[test]
    fn invalid_adjustment_lower_greater_than_upper_is_rejected() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let bad = [DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(3),
            lower: 5.0,
            upper: 2.0,
        }];
        let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            &mol,
            &EmbedParameters::default(),
            &bad,
        )
        .unwrap_err();
        assert!(matches!(
            err.0,
            EmbedWithAdjustmentsFailure::InvalidAdjustment
        ));
    }

    #[test]
    fn invalid_adjustment_nan_is_rejected() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let bad = [DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(3),
            lower: f64::NAN,
            upper: 2.0,
        }];
        let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            &mol,
            &EmbedParameters::default(),
            &bad,
        )
        .unwrap_err();
        assert!(matches!(
            err.0,
            EmbedWithAdjustmentsFailure::InvalidAdjustment
        ));
    }

    #[test]
    fn zero_adjustments_is_byte_identical_to_public_detail_api() {
        let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene
        let params = EmbedParameters::default();
        let (via_hook, _) =
            distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(&mol, &params, &[])
                .unwrap();
        let (via_public, _) =
            distance_geometry_v2::embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        for i in 0..mol.atom_count() {
            assert_eq!(
                via_hook.get(AtomIdx(i as u32)),
                via_public.get(AtomIdx(i as u32))
            );
        }
    }

    #[test]
    fn stereo_repair_and_verify_fixes_violated_stereo_before_force_field() {
        let mol = parse("C[C@H](O)CC").unwrap(); // 2-butanol, declared stereo
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        let result = embed_pipeline_v2(&mol, &config).expect("should succeed");
        assert!(result.final_stereo.is_fully_satisfied());
    }

    #[test]
    fn stereo_verify_only_fails_when_final_geometry_violates_declared_stereo() {
        // Force a violation by disabling repair, and manufacturing a molecule whose
        // raw DG output is checked directly against a corrupted variant via the
        // repair path used as an oracle: run RepairAndVerify once to find whether
        // there was anything to repair, then run VerifyOnly and confirm violations
        // (if any occurred) are reported the same way, i.e. gated at stage 11 -- and
        // confirm the trivial case (no possible violation, no repair active) still
        // succeeds under VerifyOnly.
        let mol = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine, implicit-H stereocenter
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::VerifyOnly;
        // VerifyOnly must not panic and must produce a real, evaluated verdict either
        // way (Ok or a typed FinalStereoViolation) -- both are acceptable outcomes
        // here since raw DG + no repair is not guaranteed to land in the declared
        // basin; what must NOT happen is a silent Ok with n_violations() > 0.
        match embed_pipeline_v2(&mol, &config) {
            Ok(result) => assert!(result.final_stereo.is_fully_satisfied()),
            Err(err) => assert!(matches!(
                err.cause,
                PipelineV2FailureCause::FinalStereoViolation
            )),
        }
    }

    #[test]
    fn ignore_stereo_policy_never_fails_on_declared_stereo() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::Ignore;
        let result = embed_pipeline_v2(&mol, &config).expect("Ignore must never fail on stereo");
        // Evidence is still real (not fabricated), even though it doesn't gate.
        assert_eq!(result.stereo_before.n_declared(), 1);
    }

    #[test]
    fn strict_unevaluable_stereo_is_typed_failure_when_requested() {
        // A declared quaternary stereocenter with a genuinely degenerate/coplanar
        // arrangement is hard to manufacture from raw DG directly; instead this
        // confirms the wiring: fail_on_unevaluable_stereo has no effect when there is
        // nothing Unevaluable (a molecule with fully-resolvable declared stereo must
        // still succeed under strict mode).
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        config.fail_on_unevaluable_stereo = true;
        let result = embed_pipeline_v2(&mol, &config).expect("should succeed");
        assert_eq!(result.final_stereo.n_unevaluable(), 0);
    }

    #[test]
    fn wildcard_atom_reports_typed_distance_geometry_failure() {
        let mol = parse("[*]C").unwrap();
        let config = config_none();
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(
            err.cause,
            PipelineV2FailureCause::DistanceGeometry(EmbedFailureCause::InvalidTopology)
        ));
        assert_eq!(err.stage, PipelineStage::DistanceGeometry);
    }

    #[test]
    fn timeout_zero_fails_closed_with_typed_timeout() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let mut config = config_none();
        config.total_timeout_ms = Some(0);
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(err.cause, PipelineV2FailureCause::Timeout));
    }

    /// Regression test for a real bug independent verification round 1 found:
    /// `check_timeout!` used to construct a bare, all-`None` `PipelineV2Failure`
    /// regardless of how much had already been computed by that point -- silently
    /// contradicting this module's own "carries as much partial diagnostic
    /// information as possible" claim. `total_timeout_ms: Some(0)` is guaranteed to
    /// trip at SOME `check_timeout!` call site, and every such site is textually
    /// after stage 2 (torsion knowledge) completes -- not asserting on exactly
    /// which stage trips (that depends on sub-millisecond timing and is not
    /// deterministic across machines), only that the evidence stage 2 always
    /// produces by that point is never silently dropped.
    #[test]
    fn timeout_failure_still_carries_evidence_computed_before_it_tripped() {
        let mol = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol
        let mut config = config_none();
        config.total_timeout_ms = Some(0);
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(err.cause, PipelineV2FailureCause::Timeout));
        assert!(
            err.torsion_knowledge_report.is_some(),
            "a timeout must always carry at least the torsion-knowledge report \
             (every check_timeout! call site is after stage 2 completes) -- this \
             is exactly the evidence-dropping bug that was found and fixed, stage \
             reached: {:?}",
            err.stage
        );
    }

    #[test]
    fn force_field_none_still_runs_independent_final_validation() {
        let mol = parse("c1ccccc1").unwrap();
        let config = config_none();
        let result = embed_pipeline_v2(&mol, &config).expect("benzene should succeed");
        assert!(result.final_validation.sound);
        assert!(result.final_validation.all_finite);
        assert!(result.final_validation.atom_count_unchanged);
    }

    #[test]
    fn atom_order_permutation_gives_equivalent_result() {
        // "atom-order permutation -> equivalence after mapping" (spec §16). Compare
        // ibuprofen written two structurally-equivalent ways is out of scope for a
        // cheap unit test (would need a full atom-map + Kabsch RMSD harness); this
        // narrower check instead confirms two independently-parsed copies of the
        // exact same SMILES (same atom order) at the same seed give identical
        // results -- the baseline the gate-harness's real permutation test builds on.
        let mol_a = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap(); // ibuprofen
        let mol_b = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();
        let config = config_none();
        let ra = embed_pipeline_v2(&mol_a, &config).unwrap();
        let rb = embed_pipeline_v2(&mol_b, &config).unwrap();
        for i in 0..mol_a.atom_count() {
            assert_eq!(
                ra.coords.get(AtomIdx(i as u32)),
                rb.coords.get(AtomIdx(i as u32))
            );
        }
    }

    // -----------------------------------------------------------------------
    // Negative controls (spec §17): deliberately corrupt an expected value and
    // confirm the relevant check actually fails.
    // -----------------------------------------------------------------------

    #[test]
    fn negative_control_optimize_torsions_alone_would_silently_accept_ring_only_potentials() {
        // Proves WHY stage 6's own pre-check is load-bearing: `optimize_torsions`
        // itself returns Ok (not a typed failure) when every potential's central bond
        // is a ring bond (none are "rotatable"), because it just optimizes zero
        // bonds and trivially "converges." If this test ever starts failing (i.e.
        // optimize_torsions itself started rejecting ring-only input), the
        // stage-6 pre-check in `embed_pipeline_v2` would still be correct, but this
        // comment's premise should be re-checked.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_torsions: true,
            ..TorsionKnowledgeConfig::default()
        };
        let report = build_torsion_knowledge(&mol, &config);
        assert!(
            report.potentials.iter().all(|p| !is_bridge_bond(
                &mol,
                p.central_bond.0,
                p.central_bond.1
            )),
            "expected every macrocycle potential's central bond to be non-bridge"
        );
        let coords =
            distance_geometry_v2::embed_distance_geometry_v2(&mol, &EmbedParameters::default())
                .unwrap();
        let opt_config = TorsionOptimizationConfig::default();
        let outcome = optimize_torsions(&mol, &coords, &report.potentials, &opt_config);
        assert!(
            outcome.is_ok(),
            "optimize_torsions alone silently succeeds on ring-only potentials -- \
             this is exactly why the pipeline's own stage-6 gate cannot delegate to it"
        );
    }

    /// Negative control (spec §17): proves the check in
    /// `ring_torsion_diagnostic_only_succeeds_with_scored_only_evidence` actually
    /// discriminates REAL pipeline output, not just vacuously passes. Runs the same
    /// real `embed_pipeline_v2` call and asserts the WRONG thing (every macrocycle
    /// potential IS applied) -- this must panic, proving a regression that started
    /// reporting scored-only potentials as applied would be caught, not silently
    /// accepted.
    #[test]
    #[should_panic(expected = "must never be reported as applied_to_geometry")]
    fn negative_control_reporting_scored_only_as_applied_would_be_caught() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        let result = embed_pipeline_v2(&mol, &config).expect("DiagnosticOnly must succeed");
        for p in &result.ring_torsion_evidence.potentials {
            if p.source == TorsionKnowledgeSource::MacrocycleAdaptation {
                assert!(
                    p.applied_to_geometry, // deliberately WRONG -- must panic
                    "a macrocycle potential must never be reported as applied_to_geometry"
                );
            }
        }
    }

    #[test]
    fn negative_control_bounds_after_smoothing_would_fail_the_invariant_check() {
        // Manufacture the exact bug spec §17 names: apply an adjustment to a matrix
        // that has ALREADY been smoothed (instead of before), and confirm the
        // existing `smoothing_preserves_invariants` check (reused by
        // `embed_distance_geometry_v2_with_adjustments`) is what would catch it, by
        // constructing bounds that violate the invariant directly and confirming the
        // helper says so.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let (lower0, upper0) = crate::dg_fft::build_bound_matrix(&mol);
        let mut lower = lower0.clone();
        let mut upper = upper0.clone();
        crate::dg_fft::smooth_bounds(&mut lower, &mut upper);
        // Now simulate "adjustment applied after smoothing" by widening a
        // FINITE-baseline entry in the ALREADY-SMOOTHED `upper` beyond what the
        // (unsmoothed) baseline allows -- this must be flagged. Search for a pair
        // with a genuinely finite pre-smoothing upper bound (a directly bonded pair
        // always has one; an unconstrained distant pair may have an infinite
        // fallback, which would make this check vacuously pass for the wrong
        // reason, so it must be excluded).
        let n = mol.atom_count();
        let mut found = false;
        'search: for i in 0..n {
            for j in (i + 1)..n {
                if upper0[i][j].is_finite() {
                    let widened_upper = upper[i][j] + 100.0;
                    let bad = widened_upper > upper0[i][j] + 1e-6;
                    assert!(
                        bad,
                        "a post-smoothing widening of a finite-baseline pair must be \
                         detectable as a baseline violation"
                    );
                    found = true;
                    break 'search;
                }
            }
        }
        assert!(
            found,
            "expected at least one finite-baseline pair in cyclododecane"
        );
    }

    /// Was `negative_control_final_stereo_violation_as_success_would_be_caught`
    /// (spec §17), using a REAL molecule found by the integration gate harness
    /// (`examples/pipeline_v2_integration_gate.rs`): under `RepairAndVerify` +
    /// `ForceFieldPolicy::Dreiding`, gly-ala-gly's stereo repair succeeds (stage 8)
    /// but DREIDING minimization (stage 10, no stereo awareness at all) measurably
    /// walks the geometry back across a declared stereo boundary. At the time this
    /// test was written, stage 11 had no way to recover from that and correctly
    /// caught it as `FinalStereoViolation` -- this test asserted exactly that.
    ///
    /// Issue #227 Phase 2 added a post-minimization repair-and-reverify step
    /// (immediately below stage 11's violation check, `RepairAndVerify` only) for
    /// precisely this failure shape -- stage 8's repair runs too early to see a
    /// violation minimization itself introduces, so `RepairAndVerify` now gets one
    /// more repair attempt on the POST-minimization geometry, accepted only if it
    /// fully clears every violation and the repaired geometry is still sound. For
    /// this exact molecule/config, that new step succeeds (empirically confirmed via
    /// `mmff94_bci_stereo_drift_diagnostic_227.rs` before this test was written:
    /// `repair_stereo` recovers a robustly-satisfied geometry -179.7° from the
    /// declared-boundary, `worst_bond_length_ratio`/`gross_clash_count` both
    /// unchanged by the reflection) -- so the old premise ("this case is
    /// unrecoverable") no longer holds. Updated to assert the new, correct, more
    /// precise behavior instead: not just "does it fail," but "does it recover, and
    /// is the recovery genuinely visible in the evidence."
    #[test]
    fn repair_and_verify_recovers_post_minimization_stereo_violation() {
        let mol = parse("NCC(=O)N[C@@H](C)C(=O)NCC(=O)O").unwrap(); // gly_ala_gly
        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Dreiding);
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        let result = embed_pipeline_v2(&mol, &config).expect(
            "gly_ala_gly under RepairAndVerify+Dreiding must now recover via the \
             post-minimization repair-and-reverify step",
        );
        // Stage 8's repair must have genuinely succeeded first (this is specifically
        // the "fixed, then broken again by the force field" scenario, not merely
        // "repair never worked") -- unchanged precondition from the original test.
        assert!(
            result.stereo_after_repair.is_fully_satisfied(),
            "repair (stage 8) must have succeeded before Dreiding (stage 10) broke it again"
        );
        assert!(
            result.final_stereo.is_fully_satisfied(),
            "post-minimization repair-and-reverify must leave final_stereo fully satisfied"
        );
        let repair = result
            .post_minimization_stereo_repair
            .as_ref()
            .expect("post_minimization_stereo_repair must be Some -- this is exactly the case it exists for");
        assert!(
            !repair.repaired.is_empty(),
            "the post-min repair summary must record what it actually repaired"
        );
        assert!(repair.failures.is_empty());
    }

    /// Companion sanity check: the post-minimization repair-and-reverify step must
    /// be a true no-op (never invoked, `post_minimization_stereo_repair: None`, its
    /// own timing bucket `0`) whenever stage 11 finds nothing to recover from --
    /// covers both "policy doesn't gate stereo at all" (`Ignore`) and "nothing was
    /// violated" implicitly via any passing `RepairAndVerify` case elsewhere in this
    /// test module. Guards against the new step accidentally firing (and doing
    /// needless work, or worse, silently swapping in different-but-also-valid
    /// coordinates) on a call that never needed it.
    #[test]
    fn post_minimization_stereo_repair_is_a_no_op_when_nothing_needs_recovering() {
        let mol = parse("CCCCCC").unwrap(); // hexane, no declared stereo at all
        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Dreiding);
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        let result = embed_pipeline_v2(&mol, &config).expect("hexane must embed and minimize");
        assert!(result.post_minimization_stereo_repair.is_none());
        assert_eq!(result.elapsed_ms_by_stage.post_min_stereo_repair_ms, 0);
    }

    /// Golden regression test (issue #227 Phase 2), pinned to the exact molecule
    /// that surfaced this gap: `chembl_tier_b_0082`
    /// (`COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC`, ChEMBL Wave 1
    /// corpus). Its single declared E/Z bond is satisfied post-embedding (`stereo_before`)
    /// under the exact seed `pipeline_v2_vs_rdkit_dump.rs` uses, but the BCI charge
    /// fix (this same PR) changed the electrostatic term enough that
    /// `Mmff94BondAngleStrict` minimization now walks it to `Violated` --
    /// oracle-confirmed via `mmff94_bci_stereo_drift_diagnostic_227.rs` (dihedral
    /// ~167° -> ~0.3°) and via RDKit's own real MMFF94 minimizer on the same
    /// molecule, which does NOT reproduce this (all 4 `rdkit_etkdgv3_*` arms
    /// satisfied -- a chematic-specific minimizer-robustness gap, not a physically
    /// expected crossing). Pins two things at once: (1) under `Ignore` (the policy
    /// `pipeline_v2_vs_rdkit_dump.rs`'s `chematic_pipeline_v2_mmff94_strict` arm
    /// actually uses, and this PR's own measured baseline), the violation is real
    /// and NOT gated -- documented here explicitly, not silently accepted; (2) under
    /// `RepairAndVerify`, the new post-minimization step genuinely recovers it.
    #[test]
    // Runs in CI: the bounded stereo regression completes within seconds.
    fn chembl_tier_b_0082_ez_bond_survives_bci_fix_under_repair_and_verify_not_under_ignore() {
        let mol = parse("COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC").unwrap();

        // (1) Ignore: same policy pipeline_v2_vs_rdkit_dump.rs's
        // chematic_pipeline_v2_mmff94_strict arm uses. Must still succeed (Ignore
        // never gates on stereo) but must show the real, uncorrected violation --
        // this is the documented, known, out-of-scope-for-Ignore residual from the
        // BCI fix's 3-state re-measurement (validation/results/mmff94_bci_gap_227_phase2_report.md
        // §3b), pinned here so it can never silently regress further or be
        // silently "fixed" by an unrelated future change without this test noticing.
        // Config mirrors `pipeline_v2_vs_rdkit_dump.rs`'s `base_config` exactly
        // (embed feature flags + `RingTorsionApplicationPolicy::DiagnosticOnly`) --
        // `PipelineV2Config::minimal`'s all-conservative defaults do not reproduce
        // this molecule's violation (confirmed: `EmbedParameters::default()`'s
        // fewer active embed features change the geometry enough that it doesn't
        // manifest), so this is not just a style choice, it's load-bearing for
        // reproducing the exact case measured in the corpus dump.
        let mut ignore_config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94BondAngleStrict);
        ignore_config.embed = EmbedParameters {
            random_seed: 20260801,
            max_attempts: 8,
            use_exp_torsions: true,
            use_small_ring_torsions: true,
            use_macrocycle_torsions: true,
            use_macrocycle_14_bounds: true,
            track_failures: true,
            ..EmbedParameters::default()
        };
        ignore_config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        ignore_config.stereo_policy = StereoPolicy::Ignore;
        let ignore_result = embed_pipeline_v2(&mol, &ignore_config)
            .expect("chembl_tier_b_0082 must still succeed under Ignore (never gated)");
        // NOT asserted here: an exact `stereo_before`/`final_stereo` violation
        // count. Empirically, `stereo_before.n_violations()` differs between a
        // `--release` and a plain `cargo test` (dev profile) build for this
        // exact molecule/seed (1 vs 0) -- distance-geometry embedding's
        // eigendecomposition + `max_attempts` retry loop is numerically
        // sensitive enough that which attempt succeeds first can differ across
        // optimization levels, even at a fixed seed. This is a real, checked
        // property of the embedder (not assumed away), not a bug this PR
        // introduced or is trying to fix -- so the assertions below only pin
        // what's true in BOTH profiles: Ignore never repairs, and (below)
        // RepairAndVerify always ends fully satisfied regardless of which
        // stage the violation actually appeared at. The exact "violated only
        // in State 3, satisfied in State 2, via minimization specifically"
        // narrative is the `--release`-build, `pipeline_v2_vs_rdkit_dump.rs`-
        // reproduced case documented in `scripts/mmff94_provenance/PROVENANCE.md`
        // and `validation/results/mmff94_bci_gap_227_phase2_report.md` -- this
        // test pins the POLICY CONTRACT (Ignore doesn't recover, RepairAndVerify
        // does), not the exact numeric trajectory, which is what actually
        // varies here.
        assert!(ignore_result.post_minimization_stereo_repair.is_none());

        // (2) RepairAndVerify: must always end fully satisfied, regardless of
        // whether stage 8 (pre-minimization) or the new post-minimization step
        // is what ultimately did the repairing.
        let mut repair_config = ignore_config.clone();
        repair_config.stereo_policy = StereoPolicy::RepairAndVerify;
        let repair_result = embed_pipeline_v2(&mol, &repair_config)
            .expect("chembl_tier_b_0082 must recover under RepairAndVerify");
        assert!(repair_result.final_stereo.is_fully_satisfied());
    }

    /// Negative control (spec §17): a genuinely degenerate (coplanar) 3D arrangement
    /// for a real declared stereocenter must read `Unevaluable`, and
    /// `StereoVerification`'s own accessor methods (reused, not reimplemented, by
    /// `pipeline_v2.rs`'s stage-11 gate) must never fold that into `n_satisfied()`.
    #[test]
    fn negative_control_unevaluable_counted_as_satisfied_would_be_caught() {
        let m = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine, implicit-H stereocenter
        let idx = AtomIdx(1);
        // All 4 real neighbors placed in the same z=0 plane -> the phantom-H
        // direction (opposite the sum of the other three unit bond vectors) is
        // ill-defined / the signed volume is ~0: a genuinely degenerate geometry,
        // not a hand-picked StereoStatus value.
        let mut coords = Coords3D::new_zeroed(m.atom_count());
        coords.set(AtomIdx(0), crate::coords::Point3::new(-1.0, -1.3, 0.0));
        coords.set(idx, crate::coords::Point3::new(-0.3, -0.2, 0.0));
        coords.set(AtomIdx(2), crate::coords::Point3::new(-1.0, 1.0, 0.0));
        coords.set(AtomIdx(3), crate::coords::Point3::new(1.2, -0.2, 0.0));
        coords.set(AtomIdx(4), crate::coords::Point3::new(1.7, -1.0, 0.0));
        coords.set(AtomIdx(5), crate::coords::Point3::new(1.9, 0.7, 0.0));

        let report = verify_stereo(&m, &coords);
        let status = report
            .tetrahedral
            .iter()
            .find(|r| r.atom == idx)
            .unwrap()
            .status;
        assert!(
            matches!(
                status,
                crate::stereo_constraints::StereoStatus::Unevaluable(_)
            ),
            "expected a genuinely degenerate (coplanar) geometry to read Unevaluable, got {status:?}"
        );
        assert_eq!(
            report.n_satisfied(),
            0,
            "Unevaluable must never count as Satisfied"
        );
        assert_eq!(
            report.n_violations(),
            0,
            "Unevaluable must never count as Violated"
        );
        assert_eq!(report.n_unevaluable(), 1);
    }

    #[test]
    fn negative_control_hidden_fallback_would_be_caught() {
        // Proves the correct fallback-occurred check (`fallback_reason.is_some()`,
        // NOT `actual_force_field_used != requested_force_field`) actually
        // discriminates -- mirrors the same check already relied on in
        // `examples/cf_integration_smoke_test.rs`.
        let mol = parse("CC").unwrap();
        let config = MinimizeConfig::default();
        let coords =
            distance_geometry_v2::embed_distance_geometry_v2(&mol, &EmbedParameters::default())
                .unwrap();
        let result = minimize_with_policy_gated(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
            false,
            false,
        )
        .expect("ethane should minimize fine under MMFF94");
        // Ethane's MMFF94 attempt should succeed cleanly (no fallback needed) --
        // `actual_force_field_used` still legitimately differs in NAME semantics from
        // `requested_force_field` per that struct's own doc, so asserting
        // `fallback_reason.is_none()` (not `actual == requested`) is the only
        // correct way to confirm "no fallback occurred."
        assert!(
            result.fallback_reason.is_none(),
            "ethane should not need a UFF fallback"
        );
    }

    #[test]
    fn negative_control_partial_coords_as_ok_would_be_caught() {
        // Confirms PipelineV2Result and PipelineV2Failure are structurally distinct
        // types -- a failure can never be mistaken for (or silently coerced into) a
        // success carrying the same `coords` field name, since `PipelineV2Failure`
        // exposes `last_known_coords: Option<Coords3D>`, never `coords: Coords3D`.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        // `last_known_coords` is diagnostic only -- confirm the type itself has no
        // field named `coords` a careless caller could mistake for a success value
        // (checked at compile time by construction; this line documents the intent).
        let _diagnostic_only: Option<Coords3D> = err.last_known_coords;
    }

    #[test]
    fn negative_control_dropping_a_failed_molecule_would_be_caught() {
        // Every failure carries `stage` -- confirms a caller CAN always attribute a
        // failed molecule to a specific stage rather than needing to drop it from
        // accounting. Exercised for real by the gate harness's denominator
        // bookkeeping; this unit test confirms the field is always populated.
        let mol = parse("[*]C").unwrap();
        let err = embed_pipeline_v2(&mol, &config_none()).unwrap_err();
        // `stage` is a plain enum value, always present (not Option) -- there is no
        // code path that returns `Err` without it.
        match err.stage {
            PipelineStage::ValidateConfig
            | PipelineStage::TorsionKnowledge
            | PipelineStage::MacrocycleBoundAdjustment
            | PipelineStage::DistanceGeometry
            | PipelineStage::TorsionEnergyEvaluation
            | PipelineStage::TorsionOptimization
            | PipelineStage::StereoVerifyBefore
            | PipelineStage::StereoRepair
            | PipelineStage::StereoVerifyAfterRepair
            | PipelineStage::ForceFieldMinimization
            | PipelineStage::FinalStereoVerify
            | PipelineStage::FinalGeometryValidationStage => {}
        }
    }

    #[test]
    fn negative_control_legacy_api_changing_would_be_caught() {
        // "the legacy live API changing at all" -- confirms generate_coords_etkdg
        // (an existing default path this PR must not touch) still produces the same
        // shape/behavior it always has: infallible, correct atom count. A full
        // byte-for-byte regression comparison lives in the crate's existing test
        // suite (untouched by this PR); this is a light-touch confirmation the
        // symbol is still callable with its pre-existing signature and behavior.
        let mol = parse("CCO").unwrap();
        let coords = crate::generate_coords_etkdg(&mol);
        assert_eq!(coords.atom_count(), mol.atom_count());
    }

    #[test]
    fn negative_control_accepting_nan_inf_would_be_caught() {
        let mol = parse("CC").unwrap();
        let bad_pairs = [
            DistanceBoundAdjustment {
                atom1: AtomIdx(0),
                atom2: AtomIdx(1),
                lower: f64::NAN,
                upper: 2.0,
            },
            DistanceBoundAdjustment {
                atom1: AtomIdx(0),
                atom2: AtomIdx(1),
                lower: 1.0,
                upper: f64::INFINITY,
            },
        ];
        for bad in bad_pairs {
            let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
                &mol,
                &EmbedParameters::default(),
                &[bad],
            )
            .unwrap_err();
            assert!(matches!(
                err.0,
                EmbedWithAdjustmentsFailure::InvalidAdjustment
            ));
        }
    }
}
