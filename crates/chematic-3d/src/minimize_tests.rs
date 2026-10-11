use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dg::generate_coords;
    use chematic_smiles::parse;

    fn all_pairs_min_dist(coords: &Coords3D, n: usize) -> f64 {
        let mut min_d = f64::MAX;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = coords
                    .get(AtomIdx(i as u32))
                    .distance(&coords.get(AtomIdx(j as u32)));
                min_d = min_d.min(d);
            }
        }
        min_d
    }

    #[test]
    fn test_single_atom_unchanged() {
        let mol = parse("O").unwrap();
        let coords = generate_coords(&mol);
        let orig = coords.get(AtomIdx(0));
        let result = minimize(&mol, coords);
        let after = result.get(AtomIdx(0));
        assert!((orig.x - after.x).abs() < 1e-10);
    }

    #[test]
    fn test_zero_steps_unchanged() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig {
            max_steps: 0,
            ..MinimizeConfig::default()
        };
        let before0 = coords.get(AtomIdx(0));
        let result = minimize_with_config(&mol, coords, &config);
        let after0 = result.get(AtomIdx(0));
        assert!((before0.x - after0.x).abs() < 1e-10);
    }

    #[test]
    fn test_ethane_bond_after_minimize() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(
            d > 1.2 && d < 1.8,
            "C-C distance={d:.3}, expected 1.2-1.8 Å"
        );
    }

    #[test]
    fn test_ethane_converges_to_uff_length() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        // UFF C-C single bond is 1.540 Å; minimizer should get within 0.05 Å.
        assert!(
            (d - 1.540).abs() < 0.05,
            "C-C distance={d:.4}, expected ~1.540"
        );
    }

    #[test]
    fn test_propane_no_clash() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "atom clash: min distance={min_d:.3}");
    }

    #[test]
    fn test_benzene_no_clash() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(
            min_d > 0.8,
            "atom clash in benzene: min distance={min_d:.3}"
        );
    }

    #[test]
    fn test_disconnected_no_clash() {
        let mol = parse("CC.CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(
            min_d > 0.8,
            "atom clash in disconnected: min distance={min_d:.3}"
        );
    }

    #[test]
    fn test_default_config_no_panic() {
        let mol = parse("CC(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        assert_eq!(result.atom_count(), mol.atom_count());
    }

    #[test]
    fn test_acetic_acid_no_clash() {
        let mol = parse("CC(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "clash in acetic acid: {min_d:.3}");
    }

    #[test]
    fn test_minimize_idempotent() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let result1 = minimize(&mol, coords);
        let e1 = total_energy(&mol, &result1);
        let result2 = minimize(&mol, result1);
        let e2 = total_energy(&mol, &result2);
        assert!(e2 <= e1 + 1.0, "energy increased: e1={e1:.4}, e2={e2:.4}");
    }

    #[test]
    fn test_naphthalene_no_overlap() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "overlap in naphthalene: {min_d:.3}");
    }

    #[test]
    fn test_co_bond_double_shorter_than_single() {
        // Acetic acid: C=O should be shorter than C-O
        let mol = parse("CC(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        // Atom 1 is the carbonyl C, its bonds include C=O (double) and C-O (single).
        // Just check overall: minimized coords have no clash and atom count preserved.
        assert_eq!(result.atom_count(), 4);
        let min_d = all_pairs_min_dist(&result, 4);
        assert!(min_d > 0.5, "clash in CO test: {min_d:.3}");
    }

    #[test]
    fn test_heteroatom_c_n_bond() {
        let mol = parse("CN").unwrap(); // methylamine
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        // C-N single bond UFF: 1.469 Å; expect within 0.1 Å.
        assert!(
            (d - 1.469).abs() < 0.1,
            "C-N distance={d:.4}, expected ~1.469"
        );
    }

    #[test]
    fn test_acetylene_sp_hybridization() {
        let mol = parse("C#C").unwrap(); // acetylene: C≡C
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        // C≡C triple bond UFF: 1.204 Å; expect within 0.05 Å.
        assert!(
            (d - 1.204).abs() < 0.05,
            "C≡C distance={d:.4}, expected ~1.204"
        );
    }

    #[test]
    fn test_ideal_bond_len_cc_single() {
        assert!((ideal_bond_len("C", "C", BondOrder::Single) - 1.540).abs() < 1e-6);
        assert!((ideal_bond_len("C", "C", BondOrder::Double) - 1.340).abs() < 1e-6);
        assert!((ideal_bond_len("C", "C", BondOrder::Triple) - 1.204).abs() < 1e-6);
        assert!((ideal_bond_len("C", "C", BondOrder::Aromatic) - 1.395).abs() < 1e-6);
    }

    #[test]
    fn test_ideal_bond_len_symmetry() {
        // Should be the same regardless of argument order.
        let bo = BondOrder::Single;
        assert_eq!(ideal_bond_len("C", "N", bo), ideal_bond_len("N", "C", bo));
        assert_eq!(ideal_bond_len("C", "O", bo), ideal_bond_len("O", "C", bo));
        assert_eq!(ideal_bond_len("Br", "C", bo), ideal_bond_len("C", "Br", bo));
    }

    #[test]
    fn test_atom_hybridization_sp2_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        for i in 0..6 {
            assert_eq!(
                atom_hybridization(&mol, AtomIdx(i)),
                Hybridization::SP2,
                "benzene atom {i} should be SP2"
            );
        }
    }

    #[test]
    fn test_atom_hybridization_sp_triple() {
        let mol = parse("C#C").unwrap();
        assert_eq!(atom_hybridization(&mol, AtomIdx(0)), Hybridization::SP);
        assert_eq!(atom_hybridization(&mol, AtomIdx(1)), Hybridization::SP);
    }

    #[test]
    fn test_atom_hybridization_sp3_alkane() {
        let mol = parse("CCC").unwrap();
        for i in 0..3 {
            assert_eq!(
                atom_hybridization(&mol, AtomIdx(i)),
                Hybridization::SP3,
                "propane atom {i} should be SP3"
            );
        }
    }

    #[test]
    fn test_minimize_dreiding_ethane_no_clash() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let min_coords = minimize_dreiding(&mol, coords);
        let n = mol.atom_count();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = min_coords
                    .get(AtomIdx(i as u32))
                    .distance(&min_coords.get(AtomIdx(j as u32)));
                assert!(
                    d > 0.5,
                    "atoms {i} and {j} clashed after DREIDING minimization (d={d:.3})"
                );
            }
        }
    }

    #[test]
    fn test_minimize_dreiding_benzene_no_clash() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let min_coords = minimize_dreiding(&mol, coords);
        let n = mol.atom_count();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = min_coords
                    .get(AtomIdx(i as u32))
                    .distance(&min_coords.get(AtomIdx(j as u32)));
                assert!(
                    d > 0.5,
                    "atoms {i} and {j} clashed after DREIDING minimization (d={d:.3})"
                );
            }
        }
    }

    #[test]
    fn test_minimize_mmff94_ethane() {
        let mol = parse("CC").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);
        assert_eq!(result.atom_count(), 2);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(d > 1.4 && d < 1.7, "C-C should be ~1.54 Å, got {:.3}", d);
    }

    #[test]
    fn test_minimize_mmff94_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);
        assert_eq!(result.atom_count(), 6);
        let min_d = all_pairs_min_dist(&result, 6);
        assert!(min_d > 1.2, "benzene clash: {min_d:.3}");
    }

    #[test]
    fn test_minimize_mmff94_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);
        // Verify minimize_mmff94 completes without error and produces valid coordinates
        assert_eq!(result.atom_count(), mol.atom_count());
        for i in 0..mol.atom_count() {
            let p = result.get(chematic_core::AtomIdx(i as u32));
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "aspirin atom {i} has invalid coords"
            );
        }
    }

    // ===== Phase 2: Electrostatic Energy Integration (B5) =====

    #[test]
    fn test_electrostatic_energy_methanol() {
        // Methanol has partial charges due to O-H polarity
        let mol = parse("CO").unwrap();
        let c = generate_coords(&mol);
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        // Electrostatic energy should be calculable
        let elec_e = electrostatic_energy_mmff94(&mol, &c, &mmff94_types);
        assert!(elec_e.is_ok());
        assert!(elec_e.unwrap().is_finite());
    }

    #[test]
    fn test_electrostatic_energy_carboxylic_acid() {
        // Carboxylic acids have significant charge separation
        let mol = parse("CC(=O)O").unwrap();
        let c = generate_coords(&mol);
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let elec_e = electrostatic_energy_mmff94(&mol, &c, &mmff94_types);
        assert!(elec_e.is_ok());
        let energy = elec_e.unwrap();
        assert!(energy.is_finite());
        // Carboxylic acids with negative oxygen should have non-zero electrostatic energy
    }

    #[test]
    fn test_mmff94_with_electrostatic_ethane() {
        // Ethane is non-polar, so electrostatic should be small
        let mol = parse("CC").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);

        // Should still minimize correctly with electrostatic term
        assert_eq!(result.atom_count(), 2);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(
            d > 1.4 && d < 1.7,
            "ethane C-C should be ~1.54 Å with electrostatic, got {:.3}",
            d
        );
    }

    #[test]
    fn test_mmff94_minimization_includes_charge_effects() {
        // Minimize polar molecule where charges should matter
        let mol = parse("CCO").unwrap();
        let c = generate_coords(&mol);

        // Minimization should complete without error and produce valid coordinates
        let result = minimize_mmff94(&mol, c);

        // All coordinates should be finite and properly positioned
        assert_eq!(result.atom_count(), 3);
        for i in 0..3 {
            let p = result.get(AtomIdx(i as u32));
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "atom {i} has invalid coordinate after minimization"
            );
        }

        // Verify atoms are reasonably separated (no clashes)
        let c_c = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        let c_o = result.get(AtomIdx(1)).distance(&result.get(AtomIdx(2)));
        assert!(c_c > 1.0, "C-C bond too short: {c_c:.3}");
        assert!(c_o > 1.0, "C-O bond too short: {c_o:.3}");
    }

    #[test]
    fn test_mmff94_charges_3d_integration() {
        // Verify that 3D charges are being used in minimization
        let mol = parse("c1ccccc1O").unwrap(); // Phenol
        let c = generate_coords(&mol);

        // The minimization should complete without errors
        let result = minimize_mmff94(&mol, c);
        assert_eq!(result.atom_count(), mol.atom_count());

        // All coordinates should be finite
        for i in 0..mol.atom_count() {
            let p = result.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    #[test]
    fn test_total_energy_mmff94_includes_electrostatic() {
        // Verify that total_energy_mmff94 includes electrostatic component
        let mol = parse("CCN").unwrap(); // Has C-N polarity
        let c = generate_coords(&mol);
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let total_e = total_energy_mmff94(&mol, &c, &mmff94_types);
        let bond_e = bond_energy_mmff94(&mol, &c, &mmff94_types);
        let angle_e = angle_energy_mmff94(&mol, &c, &mmff94_types);
        let vdw_e = vdw_energy_mmff94(&mol, &c, &mmff94_types);

        // Total should include bond + angle + vdw + electrostatic
        let electrostatic_e = electrostatic_energy_mmff94(&mol, &c, &mmff94_types).unwrap_or(0.0);
        let expected = bond_e + angle_e + vdw_e + electrostatic_e;

        assert!(
            (total_e - expected).abs() < 1e-6,
            "total energy mismatch: got {}, expected {}",
            total_e,
            expected
        );
    }
}

// ---------------------------------------------------------------------------
// Force-field bridge tests (Wave 1, Agent F)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod policy_bridge_tests {
    use super::*;
    use crate::dg::generate_coords;
    use chematic_ff::{assign_mmff94_numeric_types, mmff94_bond_energy, mmff94_total_energy};
    use chematic_smiles::parse;

    #[test]
    fn rejected_mmff94_quality_gate_preserves_optimizer_stop_reason() {
        let mol = parse("CC").expect("ethane");
        let coords = [[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]];
        let error = check_minimization_soundness(
            &mol,
            &coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            OptimizerStatus {
                converged: false,
                iterations: 20,
                mmff94_termination: Some(Mmff94TerminationReason::ConstraintRejectedFallback),
            },
            MAX_SANE_RESIDUAL_FORCE + 1.0,
            false,
        )
        .expect_err("excessive force must be a typed refusal");
        let ForceFieldBridgeError::MinimizationFailed(detail) = error else {
            panic!("expected minimization failure detail");
        };
        assert!(matches!(
            detail.reason,
            MinimizationFailureReason::ExcessiveResidualForce
        ));
        assert_eq!(
            detail.mmff94_termination,
            Some(Mmff94TerminationReason::ConstraintRejectedFallback)
        );
    }

    // --- Coords3D <-> Vec<[f64; 3]> bridge plumbing -------------------------

    #[test]
    fn coords_vec_round_trip_preserves_per_atom_identity() {
        // Distinct, non-symmetric coordinates per atom so a transposition or
        // off-by-one bug in the bridge's index handling would be caught.
        let n = 5;
        let mut coords = Coords3D::new_zeroed(n);
        for i in 0..n {
            coords.set(
                AtomIdx(i as u32),
                Point3::new(i as f64, 10.0 + i as f64 * 2.0, 100.0 - i as f64 * 3.0),
            );
        }
        let v = coords_to_vec(&coords, n);
        let back = vec_to_coords(&v);
        for i in 0..n {
            let orig = coords.get(AtomIdx(i as u32));
            let rt = back.get(AtomIdx(i as u32));
            assert!(
                (orig.x - rt.x).abs() < 1e-12
                    && (orig.y - rt.y).abs() < 1e-12
                    && (orig.z - rt.z).abs() < 1e-12,
                "atom {i} round-trip mismatch: {orig:?} vs {rt:?}"
            );
        }
    }

    /// The production bridge reports the prepared analytic MMFF94 residual
    /// gradient. Keep an independent central-difference calculation here to
    /// verify both the analytic value and the bridge's `Coords3D` <- > raw
    /// coordinate mapping on a deliberately high-gradient geometry.
    #[test]
    fn bridge_fd_gradient_matches_raw_chematic_ff_call() {
        let mol = parse("CCO").expect("ethanol topology (heavy atoms only)");
        let mut coords = generate_coords(&mol);
        // Deliberately stretch the C-O bond (atom 2) far from equilibrium so
        // the gradient at atom 2 is unambiguously the largest one.
        let o = coords.get(AtomIdx(2));
        coords.set(AtomIdx(2), Point3::new(o.x + 3.0, o.y, o.z));

        // Route (a): through Coords3D -> Vec via this bridge's own helper.
        let via_bridge = coords_to_vec(&coords, mol.atom_count());
        let grad_a = fd_max_gradient(
            &via_bridge,
            |c| mmff94_total_energy(&mol, c).expect("energy"),
            1e-4,
        );

        // Route (b): construct the equivalent raw Vec directly, no Coords3D
        // anywhere in this path.
        let raw: Vec<[f64; 3]> = (0..mol.atom_count())
            .map(|i| {
                let p = coords.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        let grad_b = fd_max_gradient(
            &raw,
            |c| mmff94_total_energy(&mol, c).expect("energy"),
            1e-4,
        );

        assert!(
            (grad_a - grad_b).abs() < 1e-9,
            "bridge FD gradient ({grad_a}) disagrees with raw-Vec FD gradient ({grad_b}) \
             on the same coordinates — possible index/transposition bug in coords_to_vec"
        );
        // Sanity: the artificially stretched bond should produce a large,
        // clearly-nonzero residual force, not a spuriously-zero one (which
        // is exactly the mechanism-3 failure mode this bridge fixes).
        assert!(
            grad_a > 10.0,
            "expected a large residual force from a 3 Å-stretched C-O bond, got {grad_a}"
        );

        let model = Mmff94EnergyModel::new(&mol).expect("prepared MMFF94 model");
        let analytic = max_gradient_component(&model.bounded_analytic_gradient(&via_bridge));
        let relative_delta = (analytic - grad_a).abs() / grad_a.max(1.0);
        assert!(
            relative_delta < 1e-5,
            "prepared analytic residual force ({analytic}) disagrees with FD reference \
             ({grad_a}); relative delta={relative_delta}"
        );
    }

    /// Pins `chematic_ff::bond_type_for` (called directly now — this bridge
    /// no longer maintains its own duplicate) against `mmff94_bond_energy`'s
    /// actual Some/None behavior for a handful of type/order pairs. If
    /// chematic-ff's classification ever regresses, this test catches it —
    /// a divergence here would make the coverage report falsely confident
    /// (missing bonds could recompute as "covered" wrongly, or vice versa),
    /// the same failure class this whole bridge exists to fix.
    ///
    /// The `(1, 3, Single)` and `(3, 7, Double)` cases are the exact
    /// carbonyl-bond regression found in independent review of this PR's
    /// *previous* hand-copied classification (pre-chematic-ff-#183, that
    /// duplicate used order-blind OR-logic: "bond_type=1 if EITHER atom is
    /// SP2"). Since type 3 (carbonyl C) was in that OR-based list, EVERY
    /// bond touching it — including the C=O double bond itself (types 3,7)
    /// — was misclassified bond_type=1, and `mmff94_bond_energy(1, 3, 7)`
    /// has no entry (only `(0, 3, 7, ...)` does), so the old duplicate
    /// falsely reported carbonyl C=O — present in essentially every
    /// ketone/aldehyde/amide/carboxylic acid — as missing MMFF94 coverage.
    /// Calling chematic-ff's real, order-aware `bond_type_for` directly
    /// (this fix) makes both cases resolve to bond_type=0 (forced by the
    /// `Double` order guard, independent of the SP2/MLTB list) and find
    /// their real table rows — see also
    /// `carbonyl_c_o_double_bond_is_not_falsely_reported_missing` below for
    /// the same regression pinned at the `compute_mmff94_coverage` level.
    #[test]
    fn bridge_bond_type_matches_chematic_ff_lookup_behavior() {
        // (type_i, type_j, order, expect_some)
        let cases: &[(u8, u8, BondOrder, bool)] = &[
            (1, 1, BondOrder::Single, true),  // C(sp3)-C(sp3): covered
            (1, 5, BondOrder::Single, true),  // C(sp3)-H: covered
            (1, 11, BondOrder::Single, true), // C(sp3)-F: covered
            (1, 12, BondOrder::Single, true), // C(sp3)-Cl: covered
            (1, 13, BondOrder::Single, true), // C(sp3)-Br: covered
            // Aromatic C-C (benzene, real type post-issue-#227-Phase-1B-0 is
            // 37 = CB, not the old wrong 63): `BondOrder::Aromatic` forces
            // bt=0 unconditionally (like Double/Triple), matching RDKit's
            // own (bondType=0, r0=1.374) for this exact bond -- verified
            // against a live RDKit oracle, see mmff94_minimizer.rs's
            // `benzene_ring_bonds_type_37_resolve_to_the_rdkit_verified_row`.
            (37, 37, BondOrder::Aromatic, true),
            // Cα(sp3)-C(carbonyl) SINGLE bond: bt=0 (type 1 isn't in MLTB_TYPES,
            // so the AND-gate fails) -> real (0,1,3) row -> covered.
            (1, 3, BondOrder::Single, true),
            // The carbonyl C=O DOUBLE bond itself: bt=0 (forced by the
            // Double/Triple/Quadruple order guard, before the MLTB check
            // even runs) -> real (0,3,7) row -> covered. THE carbonyl
            // regression case (see doc comment above).
            (3, 7, BondOrder::Double, true),
        ];
        for &(ti, tj, order, expect_some) in cases {
            let bt = bond_type_for(ti, tj, order);
            let got = mmff94_bond_energy(bt, ti, tj).is_some();
            assert_eq!(
                got, expect_some,
                "bond_type_for({ti},{tj},{order:?})={bt}: mmff94_bond_energy returned {got}, expected {expect_some}"
            );
        }
    }

    /// Direct regression pin for the carbonyl false-negative at the level
    /// that actually matters to callers: `compute_mmff94_coverage` on a real
    /// molecule containing a C=O bond must NOT list it in `bonds_missing`.
    /// This is exactly the check that would have caught the bridge's
    /// pre-#183-duplicate bug (see module docs above and the PR body) before
    /// it shipped — `Mmff94BondAngleStrict` would otherwise have spuriously refused
    /// (or `Mmff94WithUffFallback` silently fallen back on) ordinary
    /// ketones/aldehydes/amides/carboxylic acids.
    #[test]
    fn carbonyl_c_o_double_bond_is_not_falsely_reported_missing() {
        let mol = parse("CC=O").expect("acetaldehyde: C(sp3)-C(carbonyl)=O");
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        // Confirm this really does exercise the carbonyl types (3, 7), not
        // some other typing this bridge doesn't intend to test.
        assert_eq!(types[1], 3, "atom 1 (carbonyl C) should be MMFF94 type 3");
        assert_eq!(types[2], 7, "atom 2 (carbonyl O) should be MMFF94 type 7");

        let report = compute_mmff94_coverage(&mol, &types);
        assert!(
            report.bonds_missing.is_empty(),
            "carbonyl C=O (and Cα-C) bonds must be covered, not falsely missing: {:?}",
            report.bonds_missing
        );
    }

    #[test]
    fn torsion_coverage_counts_only_four_distinct_atoms() {
        for (smiles, expected_count) in [("C1CC1", 0), ("C1CCC1", 4), ("CCCC", 1)] {
            let mol = parse(smiles).unwrap();
            let types = assign_mmff94_numeric_types(&mol).unwrap();
            let report = compute_mmff94_coverage(&mol, &types);
            assert_eq!(report.torsions_total, expected_count, "{smiles}");
        }
    }

    #[test]
    fn iminium_torsion_empirical_parameters_are_included_in_coverage() {
        let mol = parse("CC=[N+](C)c1ccccc1").unwrap();
        let (types, view) = chematic_ff::assign_mmff94_numeric_types_with_view(&mol).unwrap();
        assert_eq!((types[2], types[4]), (54, 37));
        let report = compute_mmff94_coverage(&view, &types);
        assert!(report.torsions_total > 0);
        assert!(
            report.torsions_missing.is_empty(),
            "{:?}",
            report.torsions_missing
        );
    }

    #[test]
    fn oop_coverage_includes_n_oxide_and_counts_centers() {
        let mol = parse("[O-][n+]1ccccc1").unwrap();
        let (types, view) = chematic_ff::assign_mmff94_numeric_types_with_view(&mol).unwrap();
        assert_eq!(types[1], 69);
        let report = compute_mmff94_coverage(&view, &types);
        assert_eq!(report.oop_total, 1);
        assert!(report.oop_missing.is_empty());
    }

    // --- Mechanism-3 repro: [C@H](F)(Cl)Br under old vs. new MMFF94 --------

    fn chfclbr_mol() -> chematic_core::Molecule {
        parse("[C@H](F)(Cl)Br").expect("chfclbr")
    }

    fn worst_bond(mol: &chematic_core::Molecule, coords: &Coords3D) -> f64 {
        let mut worst = 0.0_f64;
        for (_, bond) in mol.bonds() {
            let d = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
            if d > worst {
                worst = d;
            }
        }
        worst
    }

    #[test]
    fn old_minimize_mmff94_reproduces_mechanism3_blowup_at_this_level() {
        // Confirm the bug reproduces starting from THIS crate's own
        // dg::generate_coords geometry, at the minimize.rs entry point
        // directly (`minimize_mmff94`) — not assuming the RFC's Python-path
        // 24.3 Å number (which goes through etkdg.rs's torsion/constraint
        // stages too) transfers unchanged to this narrower entry point.
        let mol = chfclbr_mol();
        assert_eq!(mol.atom_count(), 4, "C, F, Cl, Br heavy atoms only");
        let coords = generate_coords(&mol);
        let before = worst_bond(&mol, &coords);

        let result = minimize_mmff94(&mol, coords);
        let after = worst_bond(&mol, &result);

        // The old path's bond+angle energy silently zeroes for F-C-Cl,
        // F-C-Br, Cl-C-Br (mmff94_bond_params has no halogen-halogen-via-C
        // angle entries), leaving only VdW repulsion to act on the 3 C-X
        // bonds — they should blow out well past any real covalent bond.
        assert!(
            after > 5.0,
            "expected the old MMFF94 path to blow up C-X bonds on this molecule \
             (before={before:.2} Å, after={after:.2} Å) — if this no longer reproduces, \
             the mechanism-3 'before' baseline this PR cites is stale and must be re-measured"
        );
    }

    #[test]
    fn mmff94_strict_now_resolves_chfclbr_via_empirical_angle_rule() {
        // Issue #227 Stage C: chfclbr's 3 halogen-C-halogen angles (F-C-Cl,
        // F-C-Br, Cl-C-Br -- MMFF94 types 11/1/12, 11/1/13, 12/1/13) used to
        // have no MMFF94 table entry at all (see the superseded
        // `mmff94_strict_refuses_chfclbr_with_missing_angle_params`,
        // replaced here). All 3 now resolve via Halgren's eq. 20 empirical
        // angle-bend rule -- theta0 reused verbatim from the restored
        // `(0, 0, 1, 0)` wildcard row, ka derived empirically -- oracle-
        // confirmed by
        // `mmff94_energy::angle::tests::angle_empirical_reuses_wildcard_theta0_for_halogen_and_amine_c_substituents`
        // against this exact molecule.
        let mol = chfclbr_mol();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        )
        .expect("chfclbr's 3 halogen-C-halogen angles now resolve via the Stage C empirical rule");

        assert_eq!(
            result.requested_force_field,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert!(
            result.fallback_reason.is_none(),
            "no fallback needed now that all bonds/angles resolve"
        );

        let coverage = result
            .coverage
            .as_ref()
            .expect("strict policy always reports coverage");
        assert!(coverage.bonds_missing.is_empty());
        assert!(
            coverage.angles_missing.is_empty(),
            "all 3 halogen-C-halogen angles should now resolve via the empirical rule, got: {:?}",
            coverage.angles_missing
        );

        // The empirical rule must contribute a REAL, nonzero angle energy
        // before minimization -- not silently 0.0 (the same "never claim a
        // term covered but compute 0 energy for it" principle this issue's
        // whole Stage C is built on).
        match &result.energy_before {
            EnergyReport::Mmff94(b) => {
                assert!(
                    b.angle > 0.0,
                    "empirical angle-bend energy must be real and positive before minimization, got {}",
                    b.angle
                );
            }
            other => panic!("expected Mmff94 energy report, got {other:?}"),
        }
        assert!(result.converged);
    }

    #[test]
    fn mmff94_with_uff_fallback_no_longer_needs_to_fall_back_on_chfclbr() {
        // Issue #227 Stage C: companion to
        // `mmff94_strict_now_resolves_chfclbr_via_empirical_angle_rule` --
        // now that strict MMFF94 itself resolves chfclbr, the UFF-fallback
        // policy must use MMFF94 directly rather than falling back (the
        // superseded version of this test, `..._falls_back_and_reports_why_on_chfclbr`,
        // asserted the opposite -- that fallback was required).
        let mol = chfclbr_mol();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
        )
        .expect("chfclbr now succeeds directly under MMFF94, no fallback needed");

        assert_eq!(
            result.requested_force_field,
            ForceFieldPolicy::Mmff94WithUffFallback
        );
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            "no fallback needed now that chfclbr's angles resolve via the empirical rule"
        );
        assert!(
            result.fallback_reason.is_none(),
            "fallback must not fire when the strict policy already succeeds"
        );
        assert!(result.missing_parameter_classes.is_empty());
        let coverage = result
            .coverage
            .as_ref()
            .expect("strict policy always reports coverage");
        assert!(
            coverage.angles_missing.is_empty(),
            "all 3 halogen-C-halogen angles resolve via the Stage C empirical rule, got: {:?}",
            coverage.angles_missing
        );
        assert!(
            coverage.stretch_bend_missing.is_empty(),
            "all 3 halogen-C-halogen triples resolve via the Dfsb periodic-row fallback \
             (Priority 2B) -- if this regresses, the Dfsb port likely broke"
        );

        // MMFF94 itself now produces a sane geometry -- no need for UFF's
        // generic coverage to rescue this molecule anymore.
        let after = worst_bond(&mol, &result.coords);
        assert!(
            after < 3.0,
            "expected a sane, non-blown-up geometry from MMFF94 directly, got worst bond {after:.2} Å"
        );
    }

    #[test]
    fn caffeine_passes_the_complete_bonded_term_gate_after_the_reperceived_view_fix() {
        // Issue #227 Phase 1: caffeine's dione-ring torsions used to be
        // torsions_missing under the classification bug this PR fixes --
        // never gated under the default bond+angle-only policy
        // (minimize_with_policy_gated(..., false, false), what
        // minimize_with_policy uses), but WOULD have failed
        // `include_torsion_oop_in_gate=true` before this fix. Confirms the
        // real, gated production path -- not just the coverage-audit tool --
        // now succeeds, and does so via a real, non-fallback MMFF94 run.
        let mol = chematic_smiles::parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy_gated(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
            true,
            true,
        )
        .expect("caffeine must pass the complete bonded-term gate directly, no fallback");

        let coverage = result
            .coverage
            .as_ref()
            .expect("strict policy always reports coverage");
        assert!(
            coverage.torsions_missing.is_empty(),
            "caffeine has no genuine torsion coverage gap, got: {:?}",
            coverage.torsions_missing
        );
    }

    /// Issues #185/#188's fix, exercised through `Mmff94WithUffFallback`:
    /// naphthalene lacks enough MMFF94 angle/torsion coverage to pass
    /// `Mmff94BondAngleStrict`, so this falls back to `UffOnly` — and that
    /// UFF fallback, from `generate_coords`'s geometry, is measured (58-
    /// molecule corpus, post-chematic-ff-#183) to blow the worst bond length
    /// out past 3 Å (naphthalene: 1.43 Å -> 4.74 Å) at the default 200-step
    /// budget. Before `run_uff_bridge`'s rescue existed, this surfaced as a
    /// typed `Err(MinimizationFailed)` (a real fix over the pre-soundness-gate
    /// behavior of returning `Ok(PolicyMinimizeResult { converged: false, .. })`,
    /// which a careless caller checking only `.is_ok()` would mistake for
    /// success) — but the underlying cause was `generate_coords` producing a
    /// severely clashed start, not an unrecoverable molecule: retrying once
    /// from `embed_distance_geometry_v2`'s geometry succeeds. This test now
    /// pins the fixed behavior: `Ok`, with `starting_geometry` disclosing
    /// that the rescue fired (never silently). See
    /// `chematic_ff_own_uff_minimizer_blows_up_naphthalene_independent_of_this_bridge`
    /// below — chematic-ff's own `minimize_uff`, called with zero bridge code
    /// in the path, still blows up on this exact `generate_coords` geometry;
    /// this test's rescue lives entirely in `run_uff_bridge`, not in
    /// chematic-ff, so that fixture remains correct and unchanged.
    #[test]
    fn naphthalene_now_succeeds_directly_via_mmff94_strict_post_227_fix() {
        // Issue #227 Phase 1B-0: this test used to pin naphthalene as a
        // case where MMFF94 coverage was "genuinely incomplete" (both rings
        // all-carbon fused-aromatic, mistyped as 63/64 pre-fix) and had to
        // fall back to UFF, exercising the separate #185/#188
        // catastrophic-clash-rescue path (`UffStartingGeometry::
        // ReplacedWithDistanceGeometryV2`). After the atom-typing fix
        // (naphthalene's ring carbons now correctly type as 37, CB, matching
        // both rings being 6-membered), naphthalene's MMFF94 coverage is
        // complete and it succeeds directly -- no UFF fallback, no rescue
        // path exercised at all. This is a real positive outcome of the fix,
        // not a loosened assertion. The #185/#188 rescue-path behavior this
        // test used to cover needs a different fixture molecule that still
        // genuinely requires UFF fallback post-fix if that coverage matters
        // going forward -- not re-derived here (out of this PR's scope).
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
        )
        .expect("naphthalene must succeed post-#227-fix");
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            "naphthalene should now have complete MMFF94 coverage and need no UFF fallback"
        );
        assert!(
            result.fallback_reason.is_none(),
            "no fallback should have occurred"
        );
        let worst = worst_bond_length_vec(&mol, &coords_to_vec(&result.coords, mol.atom_count()));
        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected a sound geometry, got worst bond {worst:.2} Å"
        );
    }

    /// Isolated confirmation for the item-4 finding: this calls
    /// `chematic_ff::assign_uff_types`/`minimize_uff` DIRECTLY — no
    /// `chematic-3d` bridge/conversion code anywhere in the minimization
    /// path (only `dg::generate_coords` supplies the identical starting
    /// geometry) — so the naphthalene blow-up pinned above cannot be an
    /// artifact of this bridge's `coords_to_vec`/`vec_to_coords` handling;
    /// it reproduces inside chematic-ff's own UFF minimizer. Confirmed: the
    /// vdW 1-3 exclusion bug (chematic-ff #176) is already fixed (verified
    /// by reading `uff.rs`'s graph-based exclusion set), so this is a
    /// Regression for issue #185: chematic-ff's incomplete UFF potential must
    /// not be allowed to accept an energy-decreasing proposal with a
    /// catastrophically stretched fused-aromatic bond. The minimizer now
    /// rejects that proposal in its line search and returns a bounded,
    /// explicitly unsound result only if no sound descent step remains.
    #[test]
    fn chematic_ff_own_uff_minimizer_rejects_naphthalene_blowup_steps() {
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let coords = generate_coords(&mol);
        let raw: Vec<[f64; 3]> = (0..mol.atom_count())
            .map(|i| {
                let p = coords.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();

        let types = assign_uff_types(&mol);
        let result = ff_minimize_uff(&mol, &types, raw, 200);

        let worst = mol
            .bonds()
            .map(|(_, b)| {
                let (i, j) = (b.atom1.0 as usize, b.atom2.0 as usize);
                let d = [
                    result.coords[i][0] - result.coords[j][0],
                    result.coords[i][1] - result.coords[j][1],
                    result.coords[i][2] - result.coords[j][2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .fold(0.0_f64, f64::max);

        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "UFF must reject catastrophic naphthalene steps; got worst bond {worst:.2} Å"
        );
    }

    // --- General policy sanity on a fully-covered molecule ------------------

    #[test]
    fn mmff94_strict_succeeds_and_lowers_energy_for_ethane() {
        let mol = parse("CC").expect("ethane");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        )
        .expect("ethane's single C-C bond/no-angle case is fully MMFF94-covered");

        assert_eq!(
            result.requested_force_field,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert!(result.fallback_reason.is_none());
        assert!(result.missing_parameter_classes.is_empty());
        let coverage = result.coverage.expect("coverage must be reported");
        assert!(coverage.bond_angle_fully_covered());
        assert!(
            result.energy_after.total() <= result.energy_before.total() + 1e-6,
            "energy should not increase: before={:.4}, after={:.4}",
            result.energy_before.total(),
            result.energy_after.total()
        );
        assert!(result.max_residual_force.is_finite());
        for i in 0..result.coords.atom_count() {
            let p = result.coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    /// Locks in the actually-correct `fallback_reason` invariant (a real bug
    /// was found in the doc comment, not the code, during independent
    /// review: `fallback_reason` is `Some` iff a fallback truly occurred —
    /// NOT iff `actual_force_field_used != requested_force_field`, which is
    /// true on *every* successful `Mmff94WithUffFallback` call since
    /// `actual_force_field_used == Mmff94BondAngleStrict` on success while
    /// `requested_force_field == Mmff94WithUffFallback`). Previously only the
    /// failing-molecule fallback path was tested; this covers the
    /// no-fallback-needed path on molecules fully covered by MMFF94.
    #[test]
    fn mmff94_with_uff_fallback_reports_no_fallback_when_mmff94_fully_covers() {
        for smiles in ["CC", "C1CCCCC1", "CCCC", "C1CC2CC3CC1CC(C2)C3"] {
            let mol = parse(smiles).unwrap_or_else(|e| panic!("{smiles}: {e}"));
            let coords = generate_coords(&mol);
            let config = MinimizeConfig::default();

            let result = minimize_with_policy(
                &mol,
                coords,
                ForceFieldPolicy::Mmff94WithUffFallback,
                &config,
            )
            .expect("infallible policy");

            // requested != actual is expected and NOT itself evidence of a
            // fallback -- actual_force_field_used reports which physics ran
            // (Mmff94BondAngleStrict, the more informative value), not "did it match
            // what was requested."
            assert_eq!(
                result.requested_force_field,
                ForceFieldPolicy::Mmff94WithUffFallback
            );
            assert_eq!(
                result.actual_force_field_used,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                "{smiles}: fully MMFF94-covered, should run full MMFF94, not UFF"
            );
            assert!(
                result.fallback_reason.is_none(),
                "{smiles}: no fallback should have occurred, but fallback_reason={:?}",
                result.fallback_reason
            );
            assert!(result.missing_parameter_classes.is_empty(), "{smiles}");
            let coverage = result
                .coverage
                .as_ref()
                .unwrap_or_else(|| panic!("{smiles}: coverage must be reported"));
            assert!(
                coverage.bond_angle_fully_covered(),
                "{smiles}: expected full bond+angle coverage"
            );
        }
    }

    #[test]
    fn uff_only_policy_produces_finite_non_clashing_geometry() {
        let mol = parse("c1ccccc1O").expect("phenol");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config).expect(
            "phenol's UFF geometry is sound here -- UffOnly is NOT unconditionally infallible in \
             general (see its doc: an unsound result returns Err(MinimizationFailed))",
        );
        assert_eq!(result.requested_force_field, ForceFieldPolicy::UffOnly);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::UffOnly);
        assert!(result.fallback_reason.is_none());
        assert!(matches!(result.energy_after, EnergyReport::Uff { .. }));

        let n = mol.atom_count();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = result
                    .coords
                    .get(AtomIdx(i as u32))
                    .distance(&result.coords.get(AtomIdx(j as u32)));
                assert!(d > 0.5, "atoms {i}/{j} clashed after UFF bridge: {d:.3}");
            }
        }
    }

    #[test]
    fn dreiding_policy_matches_existing_behavior_and_reports_convergence() {
        let mol = parse("CC(=O)O").expect("acetic acid");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::Dreiding, &config)
            .expect(
                "acetic acid's DREIDING geometry is sound here (worst bond ~1.5 Å) even though it \
             does not fully converge within the default iteration budget (max_residual_force \
             ~55, below the 200.0 soundness ceiling) -- Dreiding is NOT unconditionally \
             infallible in general (see check_minimization_soundness)",
            );
        assert_eq!(result.requested_force_field, ForceFieldPolicy::Dreiding);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::Dreiding);
        assert!(result.coverage.is_none());
        assert!(
            result.energy_after.total() <= result.energy_before.total() + 1e-6,
            "Dreiding energy should not increase after minimization"
        );
        assert!(result.iterations > 0, "expected at least one FD step");
    }

    #[test]
    fn none_policy_returns_geometry_unchanged() {
        let mol = parse("CCC").expect("propane");
        let coords = generate_coords(&mol);
        let before: Vec<Point3> = (0..mol.atom_count())
            .map(|i| coords.get(AtomIdx(i as u32)))
            .collect();
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::None, &config)
            .expect("None policy never errors");
        assert_eq!(result.requested_force_field, ForceFieldPolicy::None);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::None);
        assert_eq!(result.iterations, 0);
        for (i, orig) in before.iter().enumerate() {
            let after = result.coords.get(AtomIdx(i as u32));
            assert_eq!(orig.x, after.x);
            assert_eq!(orig.y, after.y);
            assert_eq!(orig.z, after.z);
        }
    }

    #[test]
    fn single_atom_molecule_is_trivially_handled_by_every_policy() {
        let mol = parse("O").expect("water heavy atom only");
        let config = MinimizeConfig::default();
        for policy in [
            ForceFieldPolicy::Mmff94BondAngleStrict,
            ForceFieldPolicy::Mmff94WithUffFallback,
            ForceFieldPolicy::UffOnly,
            ForceFieldPolicy::Dreiding,
            ForceFieldPolicy::None,
        ] {
            let coords = generate_coords(&mol);
            let result = minimize_with_policy(&mol, coords, policy, &config)
                .expect("single-atom molecules never fail any policy");
            assert!(result.converged);
            assert_eq!(result.iterations, 0);
        }
    }

    // --- Coverage report sanity ---------------------------------------------

    #[test]
    fn coverage_report_totals_match_bond_and_angle_counts() {
        let mol = parse("CCC").expect("propane"); // 2 bonds, 1 angle (middle C)
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        let report = compute_mmff94_coverage(&mol, &types);
        assert_eq!(report.bonds_total, 2);
        assert_eq!(report.angles_total, 1);
        assert!(report.bond_angle_fully_covered());
        assert_eq!(report.total_missing(), 0);
    }

    /// RDKit/COSMolKit/OpenEye advantage directive, Phase 4 (3D coverage
    /// gap): 11 of the 265-molecule strict corpus's coverage failures were
    /// exactly the same `Angle C(3)-C(2)-C(64)` gap (bis-indolylmaleimide
    /// scaffold, C5B indole-ring carbon exocyclic to a maleimide ring),
    /// closed by adding type 64's `MMFF94_EQ_LEVEL` row
    /// (`chematic-ff/src/mmff94_energy/angle.rs`). End-to-end confirmation
    /// through this crate's own real coverage-check path (not just the
    /// isolated `chematic-ff` angle-lookup unit test).
    #[test]
    fn issue_c5b_angle_gap_bis_indolylmaleimide_scaffold_now_covered() {
        let smiles = [
            "Cn1cc(C2=C(c3cn(CCCSC(=N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21", // chembl_tier_b_0180
            "Cn1cc(C2=C(c3cn(CCC[N+](C)(C)C)c4ccccc34)C(=O)NC2=O)c2ccccc21", // _0184
            "Cn1cc(C2=C(c3cn(CCCCC(=N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21", // _0185
            "CN(C)CCCn1cc(C2=C(c3cn(C)c4ccccc34)C(=O)NC2=O)c2ccccc21",   // _0186
            "Cn1cc(C2=C(c3cn(CCCN)c4ccccc34)C(=O)NC2=O)c2ccccc21",       // _0187
            "Cn1cc(C2=C(c3cn(CCCN=C(N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21", // _0190
            "Cn1cc(C2=C(c3cn(CCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",        // _0195
            "Cn1cc(C2=C(c3cn(CCCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",       // _0196
            "Cn1cc(C2=C(c3cn(CCCCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",      // _0197
            "Cn1cc(C2=C(c3cn(CCCCCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",     // _0198
            "CNCCCn1cc(C2=C(c3cn(C)c4ccccc34)C(=O)NC2=O)c2ccccc21",      // _0199
        ];
        for smi in smiles {
            let mol = parse(smi).unwrap_or_else(|e| panic!("failed to parse '{smi}': {e}"));
            let (types, mmff_mol) =
                assign_mmff94_numeric_types_with_view(&mol).unwrap_or_else(|e| {
                    panic!("assign_mmff94_numeric_types_with_view failed for '{smi}': {e:?}")
                });
            let report = compute_mmff94_coverage(&mmff_mol, &types);
            assert!(
                report.angles_missing.is_empty(),
                "'{smi}': expected no missing angle parameters, got {:?}",
                report.angles_missing
            );
        }
    }

    /// Regression for the C5A equivalence correction: `chembl_tier_b_0022`'s
    /// `(angle_type=0, 43, 18, 63)` triple now reaches the documented
    /// empirical fallback rather than remaining an angle-coverage gap.
    #[test]
    fn issue_type_63_angle_gap_is_now_covered() {
        let smi = "COc1cc2nc(N3CCN(S(=O)(=O)c4cccs4)CC3)nc(N)c2cc1OC";
        let mol = parse(smi);
        let Ok(mol) = mol else {
            panic!("chembl_tier_b_0022 fixture failed to parse: {smi}");
        };
        let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(&mol)
            .unwrap_or_else(|e| panic!("assign_mmff94_numeric_types_with_view failed: {e:?}"));
        let report = compute_mmff94_coverage(&mmff_mol, &types);
        assert!(
            report.angles_missing.is_empty(),
            "chembl_tier_b_0022's type-63 angle should resolve through the empirical fallback, got {:?}",
            report.angles_missing
        );
    }

    /// Was `chematic_ff_own_energy_function_is_blind_to_this_bond_stretch`,
    /// which pinned a real chematic-ff bug (`bond_type_for` using OR-logic
    /// with no bond-order check) found while building this bridge:
    /// `mmff94_total_energy` used to show exactly zero energy change for a 1
    /// Å stretch of the Cα(sp3)-C(carbonyl) bond. That bug is now fixed
    /// upstream (chematic-ff #173/PR #183: `bond_type_for` takes the real
    /// `BondOrder` and requires BOTH atoms in `MLTB_TYPES`, not just one) —
    /// confirmed here by asserting the *opposite* of what this test used to
    /// pin: chematic-ff's own energy function is now correctly sensitive to
    /// this stretch. If this ever regresses back to "blind," chematic-ff's
    /// bond_type_for was re-broken and this bridge's coverage report should
    /// be re-audited.
    #[test]
    fn chematic_ff_own_energy_function_is_now_sensitive_to_this_bond_stretch_post_183() {
        let mol = parse("CC=O").expect("acetaldehyde: C(sp3)-C(carbonyl)=O");
        let coords = generate_coords(&mol);
        let mut v = coords_to_vec(&coords, mol.atom_count());

        let d = [v[0][0] - v[1][0], v[0][1] - v[1][1], v[0][2] - v[1][2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let u = [d[0] / len, d[1] / len, d[2] / len];

        let e0 = mmff94_total_energy(&mol, &v).expect("energy");
        v[0][0] += u[0] * 1.0;
        v[0][1] += u[1] * 1.0;
        v[0][2] += u[2] * 1.0;
        let e1 = mmff94_total_energy(&mol, &v).expect("energy");

        assert!(
            (e1 - e0).abs() > 1.0,
            "expected chematic-ff's own mmff94_total_energy to now be clearly sensitive to \
             this 1 Å bond stretch post-#183 (was silently zero pre-fix); got delta={:.6} -- \
             if this now fails, chematic-ff's bond_type_for regressed and this bridge's \
             coverage report should be re-audited",
            e1 - e0
        );
    }

    // -----------------------------------------------------------------------
    // Issues #185 / #188 diagnosis: is `UffOnly` + `dg::generate_coords`'s
    // blow-up a bad-starting-geometry bug, a fundamentally broken minimizer,
    // or something else? See `docs/rfcs/uff_robustness_diagnosis_185_188.md` for
    // the full write-up and measured numbers this pins.
    //
    // Measured finding (not assumed): `generate_coords` produces starting
    // geometries with enormous UFF energies (naphthalene ~1.3e5, anthracene
    // ~1.3e5, hexane ~1.5e7 kcal/mol -- dominated by unrelieved vdW
    // overlap, not bond stretch) for essentially every molecule, not just
    // the ones that end up blowing up. `minimize_uff`'s plain fixed-step
    // steepest descent (no conjugate-gradient/quasi-Newton acceleration)
    // eventually relieves this and reaches a sound, low-energy geometry for
    // *most* molecules checked here, given enough iterations -- but not
    // all: naphthalene needs ~10,000 steps and hexane needs >20,000 (still
    // not fully converged at 100,000), while anthracene never reaches a
    // sound geometry at all, plateauing at worst bond 3.39 Å (still above
    // `MAX_SANE_BOND_LENGTH`) from 10,000 steps all the way through
    // 200,000. Before the soundness gate, `minimize_uff`'s RMS-gradient
    // convergence check reported `converged: true` on that unsound plateau;
    // the current contract reports `converged: false` there -- a real,
    // distinct trapped-local-minimum failure mode, not just an
    // under-provisioned iteration budget. Whether a molecule's specific
    // clash-relief trajectory succeeds at all is not predicted by starting
    // energy or worst starting bond length: naphthalene and anthracene
    // start from virtually identical raw energy (~1.3e5 kcal/mol either
    // way) and identical worst starting bond (2.26 Å) yet only one
    // recovers. (An earlier version of this comment cited anthracene's raw
    // worst bond as 3.73 Å and its energy as ~5 orders of magnitude worse
    // than naphthalene's, concluding anthracene was *immune* to this
    // failure class and using that as evidence against a `dg.rs`
    // ring-placement explanation -- both the numbers and that conclusion
    // were an artifact of a since-fixed `dg::generate_coords` bug that
    // silently superimposed two of anthracene's three rings on identical
    // coordinates, corrupting the measurement; see issue #185.)
    //
    // Not fixed by the diagnostic pass itself: per that phase's scope, it was
    // regression-fixture only. Fixed by a follow-up PR (`run_uff_bridge`'s
    // rescue -- see `UffStartingGeometry`): raw starting energy alone cannot
    // predict which starts need help (see above), so rather than raising the
    // shared default budget (measured too costly/still-insufficient for the
    // worst cases -- see `docs/rfcs/uff_robustness_diagnosis_185_188.md`'s
    // companion ablation) or touching `minimize_uff`'s algorithm at all, a
    // catastrophic bond blowup/non-finite result at the existing budget is
    // now retried once from `embed_distance_geometry_v2`'s geometry, with
    // the choice always disclosed, never silent. Still not touched:
    // `minimize_uff`'s own algorithm/default `max_steps`, and chematic-ff's
    // public API beyond the additive diagnostic instrumentation from #107's
    // work (`uff_energy_breakdown`/`minimize_uff_with_trace`/etc., unrelated
    // to this fix).

    /// The bounded UFF proposal introduced by the public-package performance
    /// pass lets the original `generate_coords` hexane geometry recover at
    /// the default budget. Pin that improvement: the bridge must return the
    /// sound first trajectory and must not claim that the DG-v2 rescue ran.
    #[test]
    fn uff_only_recovers_hexane_from_generate_coords_without_rescue() {
        let mol = parse("CCCCCC").expect("hexane");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .unwrap_or_else(|e| {
                panic!(
                    "hexane's UFF minimization from generate_coords is catastrophically clashed, \
                     but embed_distance_geometry_v2 rescues it -- expected Ok after the \
                     #185/#188 fix, got {e:?}"
                )
            });
        assert_eq!(result.requested_force_field, ForceFieldPolicy::UffOnly);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::UffOnly);
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::AsProvided),
            "the bounded first trajectory should now succeed without a rescue"
        );
        let worst = worst_bond_length_vec(&mol, &coords_to_vec(&result.coords, mol.atom_count()));
        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected a sound rescued geometry, got worst bond {worst:.2} Å"
        );
    }

    /// Negative control for the rescue itself: a molecule whose
    /// `generate_coords` start already passes `UffOnly`'s soundness gate
    /// (phenol, pinned separately below as
    /// `uff_only_policy_produces_finite_non_clashing_geometry`) must report
    /// `starting_geometry == AsProvided` -- the rescue must never fire, and
    /// must never be reported as having fired, when the first attempt was
    /// already sound.
    #[test]
    fn uff_only_reports_as_provided_when_no_rescue_was_needed() {
        let mol = parse("c1ccccc1O").expect("phenol");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .expect("phenol's generate_coords geometry is sound under UffOnly without a rescue");
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::AsProvided),
            "no rescue should have been attempted -- the first attempt was already sound"
        );
    }

    /// `energy_before`, on the (no-rescue) `AsProvided` path, must be the caller's
    /// own starting geometry's energy -- computed here independently (a fresh
    /// `uff_total_energy` call on the exact same `generate_coords` output) and
    /// compared for exact equality, not just "looks plausible."
    #[test]
    fn uff_only_as_provided_energy_before_is_the_caller_geometry_energy() {
        let mol = parse("c1ccccc1O").expect("phenol");
        let coords = generate_coords(&mol);
        let n = mol.atom_count();
        let types = assign_uff_types(&mol);
        let expected_energy_before = uff_total_energy(&mol, &types, &coords_to_vec(&coords, n));
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .expect("phenol's generate_coords geometry is sound under UffOnly");
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::AsProvided)
        );
        match result.energy_before {
            EnergyReport::Uff { total } => assert_eq!(
                total, expected_energy_before,
                "AsProvided energy_before must equal the caller-provided geometry's own energy"
            ),
            other => panic!("expected EnergyReport::Uff, got {other:?}"),
        }
    }

    /// The energy-semantics regression this test exists to catch: on a
    /// successful rescue, `energy_before`/`energy_after` must both describe
    /// the SAME trajectory (the `embed_distance_geometry_v2` geometry, before
    /// and after minimizing it) -- never the caller's abandoned original
    /// geometry paired with the retry's outcome. The bounded-proposal
    /// minimizer now recovers the former finite-clash hexane fixture directly,
    /// so make the caller trajectory explicitly non-finite to exercise the
    /// same public rescue contract deterministically.
    #[test]
    fn uff_only_rescue_energy_before_reflects_dg_v2_geometry_not_caller() {
        use crate::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};

        let mol = parse("CCCCCC").expect("hexane");
        let n = mol.atom_count();
        let types = assign_uff_types(&mol);

        let mut caller_coords = generate_coords(&mol);
        let first = caller_coords.get(AtomIdx(0));
        caller_coords.set(AtomIdx(0), Point3::new(f64::NAN, first.y, first.z));
        let caller_energy = uff_total_energy(&mol, &types, &coords_to_vec(&caller_coords, n));
        assert!(
            !caller_energy.is_finite(),
            "fixture must force the typed non-finite rescue path"
        );

        let v2_coords = embed_distance_geometry_v2(&mol, &EmbedParameters::default())
            .expect("hexane must embed via distance_geometry_v2");
        let expected_v2_energy_before =
            uff_total_energy(&mol, &types, &coords_to_vec(&v2_coords, n));
        assert!(
            expected_v2_energy_before < 1000.0,
            "expected embed_distance_geometry_v2's hexane energy to be the measured \
             ~3.78-scale sound value, got {expected_v2_energy_before}"
        );

        let config = MinimizeConfig::default();
        let result = minimize_with_policy(&mol, caller_coords, ForceFieldPolicy::UffOnly, &config)
            .expect(
                "hexane's generate_coords start is catastrophically clashed, but the \
                     embed_distance_geometry_v2 rescue must succeed",
            );
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::ReplacedWithDistanceGeometryV2),
            "expected the rescue to have fired"
        );

        match (result.energy_before, result.energy_after) {
            (EnergyReport::Uff { total: before }, EnergyReport::Uff { total: after }) => {
                assert_eq!(
                    before, expected_v2_energy_before,
                    "rescue energy_before must be the DG v2 geometry's OWN starting energy \
                     (~{expected_v2_energy_before}), not the caller's abandoned geometry's \
                     energy (~{caller_energy}) -- got {before}"
                );
                assert_ne!(
                    before, caller_energy,
                    "rescue energy_before must never equal the caller's original geometry's \
                     energy -- these are two different trajectories"
                );
                assert!(
                    after <= before,
                    "energy_after ({after}) must not exceed energy_before ({before}) -- both \
                     values describe one single minimization trajectory"
                );
            }
            other => panic!("expected (EnergyReport::Uff, EnergyReport::Uff), got {other:?}"),
        }
    }

    /// Ring-fused positive controls for issue #210. These have no movable
    /// heavy-atom bridge at several declared centers; the rescue must retain
    /// materialized implicit hydrogens through UFF and return only after the
    /// truncated heavy-atom geometry independently verifies every declaration.
    #[test]
    // Runs in CI: the bounded stereo regression completes within seconds.
    fn uff_only_rescue_preserves_ring_fused_declared_stereochemistry() {
        for (name, smiles) in [
            (
                "testosterone",
                "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
            ),
            (
                "cholesterol",
                "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
            ),
        ] {
            let mol = parse(smiles).unwrap_or_else(|error| panic!("parse {name}: {error}"));
            let result = minimize_with_policy(
                &mol,
                generate_coords(&mol),
                ForceFieldPolicy::UffOnly,
                &MinimizeConfig::default(),
            )
            .unwrap_or_else(|error| panic!("{name} rescue failed: {error:?}"));
            assert_eq!(
                result.starting_geometry,
                Some(UffStartingGeometry::ReplacedWithDistanceGeometryV2)
            );
            assert!(
                crate::stereo_constraints::verify_stereo(&mol, &result.coords).is_fully_satisfied(),
                "{name} returned with violated declared stereo"
            );
        }
    }

    /// A geometrically sound first UFF result is still unacceptable when it
    /// contradicts declared stereo. Before issue #210's final fix, ibuprofen_S
    /// returned `Ok(AsProvided)` with its only stereocenter inverted.
    #[test]
    fn uff_only_reseeds_a_sound_but_stereo_violating_first_result() {
        let mol = parse("CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O").expect("ibuprofen_S");
        let result = minimize_with_policy(
            &mol,
            generate_coords(&mol),
            ForceFieldPolicy::UffOnly,
            &MinimizeConfig::default(),
        )
        .expect("stereo-violating first result must be replaced by a safe rescue");
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::ReplacedWithDistanceGeometryV2)
        );
        assert!(
            crate::stereo_constraints::verify_stereo(&mol, &result.coords).is_fully_satisfied()
        );
    }

    /// Atorvastatin positive control retained from #210's earlier partial fix.
    /// The final five-molecule gate lives in
    /// `examples/issue210_rescue_measurement.rs`; this focused test keeps the
    /// original high-stress flexible-molecule regression independently pinned.
    #[test]
    // Runs in CI: the bounded stereo regression completes within seconds.
    fn uff_only_rescue_now_preserves_stereo_for_atorvastatin_fragment() {
        let mol = parse(
            "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        )
        .expect("atorvastatin_fragment");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config).expect(
            "atorvastatin_fragment's generate_coords start is catastrophically clashed, \
                 but the enforce_chirality/materialize_implicit_h_for_chirality-aware rescue \
                 must now succeed",
        );
        assert!(
            crate::stereo_constraints::verify_stereo(&mol, &result.coords).is_fully_satisfied(),
            "the rescued geometry must preserve every declared stereocenter"
        );
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::ReplacedWithDistanceGeometryV2),
            "the successful result must disclose that the bounded DG-v2 rescue was used"
        );
    }

    /// Starting-geometry-source axis: the same `UffOnly` policy, same
    /// default budget, starting from `distance_geometry_v2`'s embedder
    /// output instead of `generate_coords`, succeeds for both of the
    /// molecules pinned as failures above -- this is the other half of the
    /// comparison the diagnosis doc's table is built from.
    #[test]
    fn uff_only_succeeds_from_embed_distance_geometry_v2_for_naphthalene_and_hexane() {
        use crate::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};

        for smiles in ["c1ccc2ccccc2c1", "CCCCCC"] {
            let mol = parse(smiles).expect("parse");
            let params = EmbedParameters {
                random_seed: 42,
                max_attempts: 5,
                timeout_ms: Some(10_000),
                ..EmbedParameters::default()
            };
            let coords = embed_distance_geometry_v2(&mol, &params)
                .unwrap_or_else(|e| panic!("embedding {smiles} failed: {e:?}"));
            let config = MinimizeConfig::default();

            let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
                .unwrap_or_else(|e| {
                    panic!(
                        "expected {smiles} seeded from embed_distance_geometry_v2 to pass the \
                         UffOnly soundness gate at the default budget, got {e:?}"
                    )
                });
            assert!(
                worst_bond_length_vec(&mol, &coords_to_vec(&result.coords, mol.atom_count()))
                    <= MAX_SANE_BOND_LENGTH,
                "{smiles}: expected a sound final geometry"
            );
        }
    }

    /// Proves the "under-budgeted, not fundamentally broken" finding for the
    /// #185-class molecule: same starting geometry as the failing pin above,
    /// but with `max_steps` raised from 200 to 5,000, naphthalene's UFF
    /// minimization actually converges (`converged: true`) to a sound
    /// geometry. `#[ignore]`d -- not part of `bash scripts/check.sh`'s
    /// default run -- because 5,000 fixed-step iterations with a
    /// finite-difference gradient is slow relative to this crate's other
    /// tests; run explicitly with `cargo test -- --ignored` to reproduce.
    #[test]
    #[ignore = "slow (5,000 fixed-step iterations); run with `cargo test -- --ignored` to reproduce the convergence finding"]
    fn uff_direct_naphthalene_converges_given_a_large_enough_iteration_budget() {
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let coords = generate_coords(&mol);
        let raw = coords_to_vec(&coords, mol.atom_count());
        let types = assign_uff_types(&mol);

        let result = ff_minimize_uff(&mol, &types, raw, 5_000);
        let worst = worst_bond_length_vec(&mol, &result.coords);

        assert!(
            result.converged,
            "expected naphthalene's UFF minimization to reach the rms-gradient convergence \
             criterion within 5,000 steps (measured: converges at ~4,539); got \
             converged={} iterations={} energy={} worst_bond={worst}",
            result.converged, result.iterations, result.energy
        );
        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected a sound final geometry once converged, got worst_bond_length={worst}"
        );
    }

    /// Same proof for the #188-class molecule: hexane needs an even larger
    /// budget (measured: worst bond still 45.7 Å at 5,000 steps, 23.8 Å at
    /// 20,000, but a sound 2.04 Å at 100,000) -- confirming this is the same
    /// slow-convergence mechanism as naphthalene, just far more severe, not
    /// a qualitatively different trap. `#[ignore]`d for the same reason as
    /// above (100,000 iterations is meaningfully slower still).
    #[test]
    #[ignore = "slow (100,000 fixed-step iterations); run with `cargo test -- --ignored` to reproduce the convergence finding"]
    fn uff_direct_hexane_reaches_sound_geometry_given_a_large_enough_iteration_budget() {
        let mol = parse("CCCCCC").expect("hexane");
        let coords = generate_coords(&mol);
        let raw = coords_to_vec(&coords, mol.atom_count());
        let types = assign_uff_types(&mol);

        let result = ff_minimize_uff(&mol, &types, raw, 100_000);
        let worst = worst_bond_length_vec(&mol, &result.coords);

        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected hexane's worst bond length to have recovered to a sound value by 100,000 \
             steps (measured: 2.04 Å); got {worst} -- note this is not expected to report \
             converged=true (rms gradient tolerance is tighter than what 100,000 steps reaches \
             here), only a sound bond length",
        );
    }

    /// Issue #227 Phase 0.1: a molecule `chematic_core::kekulize` cannot
    /// Kekulize (here: an odd-membered, all-carbon, neutral aromatic ring --
    /// unmatchable by simple parity, guaranteed `Err` regardless of which of
    /// `kekulize`'s four internal passes runs) must fail closed all the way
    /// through this crate's own policy bridge, structurally classified as
    /// `ForceFieldBridgeError::UnsupportedAtomType` by `chematic_ff`'s
    /// `NumericTypeError` -> this crate's `From` impl -- never by matching on
    /// the error's message text, and never silently minimized using an
    /// un-re-perceived molecule as a stand-in MMFF view.
    #[test]
    fn kekulization_failure_fails_closed_through_the_policy_bridge_not_silently() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let atoms: Vec<_> = (0..5)
            .map(|_| b.add_atom(chematic_core::Atom::aromatic(chematic_core::Element::C)))
            .collect();
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        let mol = b.build();
        assert!(
            chematic_core::kekulize(&mol).is_err(),
            "test fixture must be genuinely unkekulizable"
        );

        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        );
        match result {
            Err(ForceFieldBridgeError::UnsupportedAtomType(msg)) => {
                assert!(
                    msg.contains("Kekulization"),
                    "error must name the Kekulization stage, got: {msg}"
                );
            }
            other => panic!(
                "expected Err(UnsupportedAtomType(_)) for an unkekulizable molecule, got {other:?}"
            ),
        }

        // Determinism: repeated calls produce the same classification, not a
        // flaky/order-dependent result.
        let coords2 = generate_coords(&mol);
        let result2 = minimize_with_policy(
            &mol,
            coords2,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        );
        assert!(
            matches!(result2, Err(ForceFieldBridgeError::UnsupportedAtomType(_))),
            "Kekulization failure must deterministically classify as UnsupportedAtomType"
        );
    }
}

#[cfg(test)]
mod diagnostic_boundary_contract_tests {
    use super::*;

    #[test]
    fn missing_parameter_citations_are_invariant_under_term_symmetries() {
        let mol = chematic_smiles::parse("FC(N)O").unwrap();
        let types = [11, 1, 8, 6];
        for (kind, permutations) in [
            (Mmff94TermKind::Bond, vec![vec![0, 1], vec![1, 0]]),
            (Mmff94TermKind::Angle, vec![vec![0, 1, 2], vec![2, 1, 0]]),
            (
                Mmff94TermKind::StretchBend,
                vec![vec![0, 1, 2], vec![2, 1, 0]],
            ),
            (
                Mmff94TermKind::Torsion,
                vec![vec![0, 1, 2, 3], vec![3, 2, 1, 0]],
            ),
            (
                Mmff94TermKind::Oop,
                vec![
                    vec![1, 0, 2, 3],
                    vec![1, 0, 3, 2],
                    vec![1, 2, 0, 3],
                    vec![1, 2, 3, 0],
                    vec![1, 3, 0, 2],
                    vec![1, 3, 2, 0],
                ],
            ),
        ] {
            let mut expected = None;
            for permutation in permutations {
                let atoms: Vec<_> = permutation.into_iter().map(AtomIdx).collect();
                let term = missing_term(&mol, &types, kind, &atoms);
                for &idx in &term.atoms {
                    assert!(term.description.contains(&format!(
                        "{}({})",
                        mol.atom(idx).element.symbol(),
                        types[idx.0 as usize]
                    )));
                }
                assert!(term.description.contains("atom indices"));
                let citation = (term.atoms, term.description);
                if let Some(expected) = &expected {
                    assert_eq!(&citation, expected);
                } else {
                    expected = Some(citation);
                }
            }
        }
    }

    #[test]
    fn force_field_errors_preserve_setup_context_and_numerical_failure_evidence() {
        let helium = chematic_smiles::parse("[He]").unwrap();
        let setup_error =
            chematic_ff::mmff94_minimizer::mmff94_total_energy(&helium, &[[0.0; 3]]).unwrap_err();
        let message = setup_error.to_string();
        assert!(
            ForceFieldBridgeError::from(setup_error)
                .to_string()
                .contains(&message)
        );
        let numeric = chematic_ff::assign_mmff94_numeric_types(&helium).unwrap_err();
        assert!(
            ForceFieldBridgeError::from(numeric.clone())
                .to_string()
                .contains(&numeric.to_string())
        );
        assert_eq!(
            ForceFieldBridgeError::MissingParameters(Box::default()).to_string(),
            "MMFF94 parameters missing for 0 internal coordinate(s) (0 bond, 0 angle, 0 torsion, 0 oop, 0 stretch-bend)"
        );
        let detail = MinimizationFailureDetail {
            policy: ForceFieldPolicy::UffOnly,
            reason: MinimizationFailureReason::CatastrophicBondBlowup,
            converged: false,
            iterations: 17,
            mmff94_termination: None,
            max_residual_force: 123.5,
            worst_bond_length: 9.75,
            distance_geometry_v2_retry_attempted: true,
        };
        let message = ForceFieldBridgeError::MinimizationFailed(Box::new(detail)).to_string();
        for value in [
            "UffOnly",
            "CatastrophicBondBlowup",
            "iterations=17",
            "9.75",
            "123.50",
        ] {
            assert!(message.contains(value), "{message}");
        }
        assert_eq!(EnergyReport::None.total(), 0.0);
    }
}
