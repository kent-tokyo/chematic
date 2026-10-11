use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn worst_bond(mol: &Molecule, coords: &Coords3D) -> f64 {
        mol.bonds()
            .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
            .fold(0.0_f64, f64::max)
    }

    /// FORMERLY a known limitation (`dg_fft::build_bound_matrix`'s angle-constraint
    /// loop treated every pair of a center atom's neighbors as a 1-3
    /// (through-center) relationship and unconditionally tightened their bound with
    /// the *generic* ideal angle (~109.5°/120°) -- but in a 3-membered ring, that
    /// "1-3" pair is *also* a direct 1-2 bonded pair, so the generic-angle bound
    /// overwrote the correct, much shorter bond-length bound with a contradictory
    /// one). Fixed by skipping the angle-derived bound for any neighbor pair that
    /// is itself directly bonded (`dg_fft::build_bond_angle_bounds`): a 3-membered
    /// ring's three bond-length constraints already fully determine its shape, so
    /// nothing is lost. See `cyclopropane_ring_closing_bond_uses_bond_length_bound_not_angle`
    /// for the exact numbers this replaced.
    #[test]
    fn three_membered_rings_embed_successfully() {
        for smiles in [
            "C1CC1",         // cyclopropane
            "C1CO1",         // epoxide
            "C1CN1",         // aziridine
            "C1CS1",         // thiirane
            "C1CC1c1ccccc1", // cyclopropylbenzene
            "C1CC1C(=O)O",   // cyclopropanecarboxylic acid
        ] {
            let mol = parse(smiles).unwrap();
            let params = EmbedParameters::default();
            let coords = embed_distance_geometry_v2(&mol, &params)
                .unwrap_or_else(|e| panic!("{smiles}: expected a successful embed, got {e:?}"));
            let worst = worst_bond(&mol, &coords);
            assert!(worst < 2.0, "{smiles}: worst bond length {worst}");
        }
        // Controls: 4- and 5-membered rings are unaffected (this is specific to
        // 3-membered rings, not "any small ring").
        for smiles in ["C1CCC1", "C1CCCC1"] {
            let mol = parse(smiles).unwrap();
            let params = EmbedParameters::default();
            assert!(
                embed_distance_geometry_v2(&mol, &params).is_ok(),
                "{smiles}: 4/5-membered rings should still embed fine"
            );
        }
    }

    #[test]
    fn cyclopropane_ring_closing_bond_uses_bond_length_bound_not_angle() {
        // The ring-closing C-C pair must keep the tight bond-length bound (upper
        // ≈ 1.59 Å) rather than being overwritten by the generic ~109.5°-angle-
        // derived bound (which previously forced lower ≈ 2.414 Å > upper, the
        // exact contradiction `three_membered_rings_embed_successfully` used to
        // fail closed on).
        let mol = parse("C1CC1").unwrap();
        let (lower, upper) = crate::dg_fft::build_bound_matrix(&mol);
        let mut found = false;
        for (_, bond) in mol.bonds() {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            assert!(
                lower[i][j] <= upper[i][j],
                "bonded pair ({i},{j}): lower={} > upper={}",
                lower[i][j],
                upper[i][j]
            );
            assert!((upper[i][j] - 1.590).abs() < 0.01, "upper={}", upper[i][j]);
            found = true;
        }
        assert!(found, "cyclopropane has no bonds?");
    }

    #[test]
    fn ethane_bond_length_reasonable() {
        let mol = parse("CC").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).expect("ethane should embed");
        let d = coords.get(AtomIdx(0)).distance(&coords.get(AtomIdx(1)));
        assert!((d - 1.54).abs() < 0.2, "C-C distance {d}");
    }

    #[test]
    fn benzene_ring_closes_and_bonds_reasonable() {
        let mol = parse("c1ccccc1").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).expect("benzene should embed");
        let worst = worst_bond(&mol, &coords);
        assert!(worst < 2.0, "benzene worst bond {worst}");
    }

    #[test]
    fn same_seed_reproducible() {
        let mol = parse("CCCCCCCCCC").unwrap(); // decane
        let params = EmbedParameters {
            random_seed: 12345,
            ..EmbedParameters::default()
        };
        let c1 = embed_distance_geometry_v2(&mol, &params).unwrap();
        let c2 = embed_distance_geometry_v2(&mol, &params).unwrap();
        for i in 0..mol.atom_count() {
            let p1 = c1.get(AtomIdx(i as u32));
            let p2 = c2.get(AtomIdx(i as u32));
            assert_eq!(
                p1, p2,
                "same seed must reproduce identical coords at atom {i}"
            );
        }
    }

    #[test]
    fn different_seed_gives_different_output() {
        let mol = parse("CCCCCCCCCC").unwrap(); // decane
        // Seeds 0/1 specifically -- under the old buggy `Prng::from_seed`
        // (`Self(seed | 1)`), 0 and 1 both mapped to the identical internal
        // state, so this exact pair is what would have caught that
        // regression. Seeds 1/2 (the original choice here) aren't an
        // aliased pair even under the old bug, so they couldn't have
        // detected it -- confirmed by independent verification during
        // review.
        let params_a = EmbedParameters {
            random_seed: 0,
            ..EmbedParameters::default()
        };
        let params_b = EmbedParameters {
            random_seed: 1,
            ..EmbedParameters::default()
        };
        let c1 = embed_distance_geometry_v2(&mol, &params_a).unwrap();
        let c2 = embed_distance_geometry_v2(&mol, &params_b).unwrap();
        let mut any_diff = false;
        for i in 0..mol.atom_count() {
            let p1 = c1.get(AtomIdx(i as u32));
            let p2 = c2.get(AtomIdx(i as u32));
            if p1 != p2 {
                any_diff = true;
            }
        }
        assert!(
            any_diff,
            "different seeds should not reproduce identical coords"
        );
    }

    #[test]
    fn track_failures_defaults_to_no_breakdown() {
        let mol = parse("CC").unwrap();
        let params = EmbedParameters::default(); // track_failures: false
        let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        assert!(stats.failure_counts.is_empty());
        assert_eq!(stats.attempts_used, 1);
    }

    #[test]
    fn track_failures_records_real_counts_on_forced_failure() {
        // O=C1CC[C@H]2CCC[C@H]12: both declared stereocenters are ring-fused with no
        // acyclic bridge substituent, so `repair_stereo` structurally cannot fix
        // either one -- empirically confirmed to fail all of seeds 0..8 with
        // max_attempts=1 (verified while writing this test, not assumed). Exercises
        // the genuine per-attempt failure path: every attempt is actually made and
        // recorded, unlike the old immediate-refusal behavior this replaced.
        let mol = parse("O=C1CC[C@H]2CCC[C@H]12").unwrap();
        let params = EmbedParameters {
            random_seed: 0,
            max_attempts: 3,
            enforce_chirality: true,
            track_failures: true,
            ..EmbedParameters::default()
        };
        let err = embed_distance_geometry_v2_detail(&mol, &params).unwrap_err();
        assert_eq!(err.0, EmbedFailureCause::StereoConstraintFailed);
        assert_eq!(
            err.1.attempts_used, 3,
            "every attempt must be genuinely made, not skipped"
        );
        assert_eq!(
            err.1
                .failure_counts
                .get(&EmbedFailureCause::StereoConstraintFailed),
            Some(&3),
            "each of the 3 attempts must be individually recorded"
        );
        assert!(err.1.stereo_repair_attempted);
    }

    #[test]
    fn enforce_chirality_repairs_and_succeeds_for_a_bridge_eligible_stereocenter() {
        // 2-butanol has one declared stereocenter with a plain acyclic substituent
        // (the ethyl group), so `repair_stereo` can always fix a wrong raw embedding
        // -- but `refine_coords`'s bounds-driven relaxation that runs after the
        // repair has no chirality awareness of its own and can occasionally flip a
        // correctly-repaired center back to `Violated` on a given attempt
        // (empirically confirmed for this exact molecule at random_seed=2,
        // max_attempts=1, while designing this: repair succeeds pre-refine, then
        // POST-REFINE is Violated again). That is exactly what the per-attempt
        // retry loop (`max_attempts`, distinct derived seed each try) exists to
        // absorb -- see the module doc's "`enforce_chirality`" section -- so this
        // asserts the retry loop delivers a satisfying geometry, not that any
        // single attempt does.
        let mol = parse("C[C@H](O)CC").unwrap();
        for seed in 0..20u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                ..EmbedParameters::default()
            };
            let (coords, _stats) =
                embed_distance_geometry_v2_detail(&mol, &params).unwrap_or_else(|e| {
                    panic!("seed {seed} must succeed within 8 attempts, got {e:?}")
                });
            assert!(
                verify_stereo(&mol, &coords).is_fully_satisfied(),
                "seed {seed}: returned geometry must satisfy declared stereo"
            );
        }
    }

    #[test]
    fn enforce_chirality_false_is_unaffected_by_declared_stereo() {
        // Same molecule/seeds as the ring-fused failing case above, but with
        // enforce_chirality left at its default (false) -- must always succeed
        // (no stereo checking at all), proving the new logic is fully gated behind
        // the flag and doesn't leak into the default path.
        let mol = parse("O=C1CC[C@H]2CCC[C@H]12").unwrap();
        for seed in 0..8u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 1,
                enforce_chirality: false,
                ..EmbedParameters::default()
            };
            let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params)
                .unwrap_or_else(|e| panic!("seed {seed} must succeed, got {e:?}"));
            assert!(!stats.stereo_repair_attempted);
        }
    }

    #[test]
    fn enforce_chirality_noop_without_declared_stereo() {
        let mol = parse("CCCC").unwrap(); // butane, no stereo at all
        let params = EmbedParameters {
            enforce_chirality: true,
            ..EmbedParameters::default()
        };
        let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        assert!(!stats.stereo_repair_attempted);
    }

    #[test]
    fn atom_limit_exceeded_is_typed_not_silent() {
        // Build a long acyclic chain past DG_MAX_ATOMS to hit the typed limit path
        // without needing to actually construct 500+ Bonds by hand -- reuse the real
        // limit constant so this test tracks it if it ever changes.
        let smiles: String = "C".repeat(DG_MAX_ATOMS + 1);
        let mol = parse(&smiles).unwrap();
        let params = EmbedParameters::default();
        let err = embed_distance_geometry_v2(&mol, &params).unwrap_err();
        assert_eq!(err, EmbedFailureCause::AtomLimitExceeded);
    }

    #[test]
    fn use_random_coords_path_still_produces_valid_geometry() {
        let mol = parse("c1ccccc1").unwrap();
        let params = EmbedParameters {
            use_random_coords: true,
            ..EmbedParameters::default()
        };
        let (coords, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        assert!(stats.used_random_coords);
        for i in 0..coords.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    #[test]
    fn wildcard_atom_is_invalid_topology() {
        let mol = parse("[*]C").unwrap();
        let params = EmbedParameters::default();
        let err = embed_distance_geometry_v2(&mol, &params).unwrap_err();
        assert_eq!(err, EmbedFailureCause::InvalidTopology);
    }

    #[test]
    fn empty_molecule_is_trivially_ok() {
        let mol = chematic_core::MoleculeBuilder::new().build();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        assert_eq!(coords.atom_count(), 0);
    }

    #[test]
    fn smoothing_invariants_hold_on_a_real_molecule() {
        // Direct check of the invariant helper against dg_fft's own bound matrix +
        // smoothing on a fused-ring molecule (naphthalene), matching the Coordinator's
        // requested "before vs after" check.
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let (lower0, upper0) = crate::dg_fft::build_bound_matrix(&mol);
        let mut lower = lower0.clone();
        let mut upper = upper0.clone();
        crate::dg_fft::smooth_bounds(&mut lower, &mut upper);
        assert!(smoothing_preserves_invariants(
            &lower0, &upper0, &lower, &upper
        ));
    }

    #[test]
    fn negative_eigenvalues_are_reported_not_silently_dropped() {
        // Not every molecule produces negative eigenvalues, but the field must always
        // be *present* and consistent (0 negative eigenvalues => magnitude 0.0).
        let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene: fused ring, good stress case
        let params = EmbedParameters::default();
        let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        if stats.negative_eigenvalues_beyond_embedding_dim == 0 {
            assert_eq!(stats.max_negative_eigenvalue_magnitude, 0.0);
        } else {
            assert!(stats.max_negative_eigenvalue_magnitude > 0.0);
        }
    }

    #[test]
    fn bounds_conformance_reports_nonzero_residual_on_a_real_molecule() {
        // refine_coords is a heuristic, not a guaranteed-convergent solver -- confirm
        // bounds_conformance actually measures real (nonzero) residual violation on a
        // stress case, not just always reporting a trivial 0/0.
        let mol = parse("CN(C)CCOC(c1ccccc1)c1ccccc1").unwrap(); // diphenhydramine
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        let bc = bounds_conformance(&mol, &coords);
        assert!(bc.n_pairs > 0);
        // Not asserting a specific violation count (that's this PR's honestly-reported
        // measured number, not a hardcoded expectation) -- just that the function
        // computes something real and internally consistent.
        assert!(bc.n_violations <= bc.n_pairs);
        if bc.n_violations == 0 {
            assert_eq!(bc.max_rel_violation, 0.0);
        } else {
            assert!(bc.max_rel_violation > 0.0);
        }
    }

    /// Stochastic gate (4-layer acceptance gate, layer 4): "uniqueness of derived
    /// per-attempt seed streams" -- if `max_attempts` allows multiple internal
    /// retries, each attempt must draw from a genuinely different seed, not
    /// accidentally repeat one (which would make a retry after a bad stochastic
    /// draw just redraw the same bad geometry). Checked directly against
    /// `derive_attempt_seed` rather than at the full-embedding level, since this is
    /// an internal-derivation invariant, not an observable-geometry one.
    #[test]
    fn derive_attempt_seed_is_distinct_across_attempts() {
        use std::collections::HashSet;
        for base in [0u64, 1, 42, 0xC0FF_EE42_D157_6E02, u64::MAX] {
            let seeds: HashSet<u64> = (0..32)
                .map(|attempt| derive_attempt_seed(base, attempt))
                .collect();
            assert_eq!(
                seeds.len(),
                32,
                "base seed {base:#x}: expected 32 distinct per-attempt seeds, got {}",
                seeds.len()
            );
        }
    }

    // -------------------------------------------------------------------------
    // `apply_declared_ez_bounds` (issue #285 but2ene_Z root-cause fix)
    // -------------------------------------------------------------------------

    /// Primary gate: the RAW (`max_attempts: 1`, no retry) embedding itself must land
    /// in the declared-E/Z-compatible region, for every SMILES notation of the same
    /// declared stereochemistry, across many seeds. Before this fix, but2ene_Z's raw
    /// embedding was deterministically wrong-signed for all of 10 tested seeds (see
    /// the PR body's root-cause diagnosis) -- `max_attempts` retry could not help,
    /// because the bound matrix itself excluded the correct geometry. A test that only
    /// checked "does `embed_distance_geometry_v2_detail` eventually succeed" (its
    /// default `max_attempts: 8`) would not catch a regression back to that state, so
    /// this asserts every individual raw attempt is already correct.
    #[test]
    fn ez_bounds_but2ene_all_notation_variants_raw_embed_all_seeds_satisfied() {
        let cases: &[(&str, &str)] = &[
            ("Z, C/C=C\\C", r"C/C=C\C"),
            ("Z, C\\C=C/C", r"C\C=C/C"),
            ("E, C/C=C/C", "C/C=C/C"),
            ("E, C\\C=C\\C", r"C\C=C\C"),
        ];
        for (label, smiles) in cases.iter().copied() {
            let mol = parse(smiles).unwrap();
            for seed in 0..10u64 {
                let params = EmbedParameters {
                    random_seed: seed,
                    max_attempts: 1,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                let (coords, _stats) = embed_distance_geometry_v2_detail(&mol, &params)
                    .unwrap_or_else(|e| {
                        panic!("but2ene ({label}) seed {seed}: raw embed must succeed, got {e:?}")
                    });
                assert!(
                    verify_stereo(&mol, &coords).is_fully_satisfied(),
                    "but2ene ({label}) seed {seed}: raw embed must already satisfy declared E/Z"
                );
            }
        }
    }

    /// Same raw-embed-all-seeds gate, generalized across the corpus's other small/
    /// rigid declared-E/Z alkenes (`scripts/etkdg_vs_rdkit_gap.py::CORPUS`'s
    /// `alkene_ez` group minus `cinnamic_acid_E`, see the next test) -- confirms the
    /// fix isn't a but2ene-specific coincidence.
    #[test]
    fn ez_bounds_rigid_alkene_corpus_raw_embed_all_seeds_satisfied() {
        let cases: &[(&str, &str)] = &[
            ("but2ene_E", "C/C=C/C"),
            ("but2ene_Z", r"C/C=C\C"),
            ("chloropropene_E", "C(/C=C/C)Cl"),
            ("chloropropene_Z", r"C(/C=C\C)Cl"),
            ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
            ("pent2ene_E", "CC/C=C/C"),
            ("pent2ene_Z", r"CC/C=C\C"),
        ];
        for (name, smiles) in cases.iter().copied() {
            let mol = parse(smiles).unwrap();
            for seed in 0..10u64 {
                let params = EmbedParameters {
                    random_seed: seed,
                    max_attempts: 1,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                let (coords, _stats) = embed_distance_geometry_v2_detail(&mol, &params)
                    .unwrap_or_else(|e| {
                        panic!("{name} seed {seed}: raw embed must succeed, got {e:?}")
                    });
                assert!(
                    verify_stereo(&mol, &coords).is_fully_satisfied(),
                    "{name} seed {seed}: raw embed must already satisfy declared E/Z"
                );
            }
        }
    }

    /// `cinnamic_acid_E` is the corpus's one alkene with genuine, pre-existing
    /// per-seed conformational variability (measured before this fix, see the module
    /// doc's original but2ene_Z diagnosis and the PR body: 3/10 raw sign-flip rate,
    /// larger/flexible molecule, NOT the small-rigid-bound-infeasibility mechanism
    /// this fix addresses) -- raw embedding is not expected to be all-seed-green for
    /// it, and this fix does not change that (its declared configuration is E, whose
    /// analytic 1-4 distance already sat inside the pre-existing VDW-derived bound).
    /// What must still hold: the ordinary retry loop (`max_attempts: 8`, the default)
    /// resolves it at the SAME rate as before this fix (independently measured on
    /// unmodified `main`, seeds 0..5: `[true, true, true, false, true]`, 4/5) -- a
    /// control against this PR accidentally making the flexible-molecule case worse
    /// while fixing the rigid one. Asserting the exact pre-fix count, not just
    /// "mostly succeeds", so a regression down to e.g. 3/5 or 2/5 is still caught.
    #[test]
    fn ez_bounds_cinnamic_acid_e_retry_loop_still_resolves_flexible_variance() {
        let mol = parse("OC(=O)/C=C/c1ccccc1").unwrap();
        let mut passes = 0usize;
        for seed in 0..5u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                ..EmbedParameters::default()
            };
            if embed_distance_geometry_v2_detail(&mol, &params).is_ok() {
                passes += 1;
            }
        }
        // Was 4/5 before issue #291's post-refinement repair safety net
        // (`try_embed_once`'s second `repair_stereo` pass when the pre-refinement
        // repair doesn't survive `refine_coords`) -- that net recovers the one
        // seed that used to fail here too, independent of ring-fused centers.
        assert_eq!(
            passes, 5,
            "cinnamic_acid_E: expected 5/5 seeds to resolve within 8 attempts \
             (post-refinement repair safety net, issue #291), got {passes}/5"
        );
    }

    /// Issue #285's two named fixtures (`chembl_tier_b_0126`/`chembl_tier_b_0168`,
    /// each with one declared E double bond alongside two declared tetrahedral
    /// centers, one ring-fused): confirms this fix generalizes to a realistic
    /// drug-like molecule, not just an isolated 4-atom toy alkene, and doesn't break
    /// the standard retry loop's ability to resolve the (unrelated, unaffected by
    /// this PR) tetrahedral-center variance those molecules still have. Exact counts
    /// independently measured on unmodified `main`, seeds 0..5, `max_attempts: 8`:
    /// `chembl_tier_b_0126` 5/5, `chembl_tier_b_0168` 4/5 (`[true, true, true, true,
    /// false]`) -- asserted exactly, not just "mostly succeeds", so a regression is
    /// still caught. Re-measured after issue #291's post-refinement repair safety
    /// net (`try_embed_once`'s second `repair_stereo` pass): `chembl_tier_b_0168`'s
    /// one previously-failing seed now recovers too, so both are 5/5.
    #[test]
    fn ez_bounds_chembl_tier_b_0126_0168_retry_loop_succeeds() {
        let cases: &[(&str, &str, usize)] = &[
            (
                "chembl_tier_b_0126",
                "CC(=O)/C=C/CC1C(=O)N2[C@@H](C(=O)O)C(C)(C)S(=O)(=O)[C@@H]12",
                5,
            ),
            (
                "chembl_tier_b_0168",
                "CC(=O)/C=C/CC1C(=O)N2[C@@H](C(=O)O)C(C)(C)S(=O)(=O)[C@H]12",
                5,
            ),
        ];
        for (name, smiles, expected_passes) in cases.iter().copied() {
            let mol = parse(smiles).unwrap();
            let mut passes = 0usize;
            for seed in 0..5u64 {
                let params = EmbedParameters {
                    random_seed: seed,
                    max_attempts: 8,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                if embed_distance_geometry_v2_detail(&mol, &params).is_ok() {
                    passes += 1;
                }
            }
            assert_eq!(
                passes, expected_passes,
                "{name}: expected {expected_passes}/5 seeds to resolve within 8 attempts \
                 (unchanged pre-existing tetrahedral-center variance), got {passes}/5"
            );
        }
    }

    /// Full re-measure of the 29 stereo-bearing molecules from
    /// `scripts/etkdg_vs_rdkit_gap.py::CORPUS` (same protocol as the module doc's
    /// original `enforce_chirality` measurement: 5 base seeds, `max_attempts: 8`),
    /// checked against the exact before-this-fix pass/fail set (independently
    /// re-measured on unmodified `main` while designing this fix, see the PR body) --
    /// not a re-derived guess. Two things must both hold:
    /// - `but2ene_Z` moves from failing (1/5 base seeds before) to fully passing
    ///   (5/5) -- the fix actually fires on the corpus, not just the isolated test
    ///   fixtures above.
    /// - Every molecule that fully passed (5/5) *before* this fix still fully passes
    ///   *after* -- zero newly-broken molecules. `testosterone`/`cholesterol` (ring-
    ///   fused tetrahedral centers) and `cinnamic_acid_E` (flexible-molecule variance)
    ///   are excluded from that "before" set on purpose: this PR is declared-scoped to
    ///   E/Z bounds only (see the module doc), so their pre-existing, unrelated
    ///   failure/partial-pass modes are expected to be unaffected, not fixed here.
    #[test]
    fn ez_bounds_29_corpus_regression_but2ene_z_fixed_nothing_else_broken() {
        const CORPUS_29: &[(&str, &str)] = &[
            ("l_alanine", "N[C@@H](C)C(=O)O"),
            ("d_alanine", "N[C@H](C)C(=O)O"),
            ("l_serine", "N[C@@H](CO)C(=O)O"),
            ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O"),
            ("2_butanol_R", "C[C@H](O)CC"),
            ("2_butanol_S", "C[C@@H](O)CC"),
            ("2_chlorobutane_R", "C[C@H](Cl)CC"),
            ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
            ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
            ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O"),
            ("chfclbr_R", "[C@H](F)(Cl)Br"),
            ("chfclbr_S", "[C@@H](F)(Cl)Br"),
            ("quaternary_1_R", "[C@](F)(Cl)(Br)I"),
            ("quaternary_1_S", "[C@@](F)(Cl)(Br)I"),
            ("quaternary_2_R", "[C@](C)(N)(O)F"),
            ("quaternary_2_S", "[C@@](C)(N)(O)F"),
            ("but2ene_E", "C/C=C/C"),
            ("but2ene_Z", r"C/C=C\C"),
            ("chloropropene_E", "C(/C=C/C)Cl"),
            ("chloropropene_Z", r"C(/C=C\C)Cl"),
            ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1"),
            ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
            ("pent2ene_E", "CC/C=C/C"),
            ("pent2ene_Z", r"CC/C=C\C"),
            (
                "penicillin_core",
                "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
            ),
            (
                "testosterone",
                "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
            ),
            (
                "cholesterol",
                "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
            ),
            (
                "atorvastatin_fragment",
                "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
            ),
            ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O"),
        ];

        // Independently re-measured on unmodified `main` (before this fix) with the
        // identical protocol below -- not a re-derived guess. Matches the module
        // doc's documented 25-27/29 range (this run: 25/29).
        const FULLY_PASSED_BEFORE_FIX: &[&str] = &[
            "2_butanol_R",
            "2_butanol_S",
            "2_chlorobutane_R",
            "atorvastatin_fragment",
            "but2ene_E",
            "chfclbr_R",
            "chfclbr_S",
            "chloropropene_E",
            "chloropropene_Z",
            "cinnamic_acid_Z",
            "d_alanine",
            "gly_ala_gly",
            "ibuprofen_S",
            "l_alanine",
            "l_serine",
            "l_threonine",
            "menthol",
            "naproxen_S",
            "penicillin_core",
            "pent2ene_E",
            "pent2ene_Z",
            "quaternary_1_R",
            "quaternary_1_S",
            "quaternary_2_R",
            "quaternary_2_S",
        ];

        let mut fully_passed_after: Vec<&str> = Vec::new();
        let mut but2ene_z_seeds_passed = 0usize;
        for (name, smiles) in CORPUS_29.iter().copied() {
            let mol = parse(smiles).unwrap();
            let mut passes = 0usize;
            for base_seed in 0..5u64 {
                let params = EmbedParameters {
                    random_seed: base_seed,
                    max_attempts: 8,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                if embed_distance_geometry_v2_detail(&mol, &params).is_ok() {
                    passes += 1;
                }
            }
            if name == "but2ene_Z" {
                but2ene_z_seeds_passed = passes;
            }
            if passes == 5 {
                fully_passed_after.push(name);
            }
        }

        assert_eq!(
            but2ene_z_seeds_passed, 5,
            "but2ene_Z must now pass all 5 base seeds (was 1/5 before this fix)"
        );
        for &name in FULLY_PASSED_BEFORE_FIX {
            assert!(
                fully_passed_after.contains(&name),
                "{name} fully passed before this fix and must still fully pass after -- \
                 regression introduced by the declared-E/Z bound change"
            );
        }
    }

    /// Issue #291 residual (testosterone, one of the two ring-fused molecules
    /// `repair_tetrahedral_center`'s bridge-eligibility check alone cannot fix, since
    /// every real neighbor of the affected centers is a ring atom): confirms
    /// `materialize_implicit_h_for_chirality` combined with the post-refinement
    /// repair safety net in `try_embed_once` above actually delivers *correct*
    /// geometry, not just an internally-reported "satisfied". Checked against
    /// `chematic_chem::assign_cip` (declared) vs `stereo3d::assign_stereo_from_3d`
    /// (perceived from the *returned* coordinates) -- a genuinely different code
    /// path than this module's own `verify_stereo`, so it can't be fooled the same
    /// way a bug in `verify_stereo`'s own input (e.g. the now-fixed `add_hydrogens`
    /// stereo-order loss) could fool `verify_stereo` alone.
    ///
    /// Without the post-refinement safety net, every one of these 5 seeds failed
    /// with `StereoConstraintFailed` (measured directly while diagnosing this) even
    /// though `repair_stereo` alone, called directly on a raw embed, reached 6/6 on
    /// every seed -- `refine_coords`'s chirality-blind bound correction was undoing
    /// the pre-refinement repair every time.
    #[test]
    fn materialize_implicit_h_for_chirality_fixes_testosterone_with_correct_geometry() {
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let declared = chematic_chem::assign_cip(&mol);
        assert_eq!(
            declared.assignments.len(),
            6,
            "sanity: 6 declared stereocenters"
        );

        let mut passes = 0usize;
        for seed in 0..5u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                materialize_implicit_h_for_chirality: true,
                ..EmbedParameters::default()
            };
            let coords = match embed_distance_geometry_v2_detail(&mol, &params) {
                Ok((coords, _)) => coords,
                Err(_) => continue,
            };
            // `assign_stereo_from_3d` only assigns centers with exactly 4 heavy-atom
            // neighbors (no implicit H) -- skip declared centers it can't independently
            // check (e.g. any with an implicit H), don't treat that as a mismatch.
            let perceived = crate::stereo3d::assign_stereo_from_3d(&mol, &coords);
            for &(idx, code) in &declared.assignments {
                if let Some(perceived_code) = perceived.get(idx) {
                    assert_eq!(
                        perceived_code, code,
                        "testosterone seed={seed}: atom {idx:?} declared {code:?} but \
                         3D-perceived {perceived_code:?} -- embed reported success with \
                         wrong chirality"
                    );
                }
            }
            passes += 1;
        }
        assert!(
            passes >= 4,
            "testosterone: expected at least 4/5 seeds to embed successfully with \
             materialize_implicit_h_for_chirality (post-refinement safety net), got {passes}/5"
        );
    }
}
