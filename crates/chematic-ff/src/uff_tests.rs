use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn assign_types_ethanol() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        assert_eq!(types.len(), 3);
        // C sp3 → C_3, O sp3 → O_3
        let type_map: std::collections::HashMap<_, _> = types.into_iter().collect();
        for (_, atom) in mol.atoms() {
            let idx = mol
                .atoms()
                .find(|(_, a)| a.element == atom.element)
                .map(|(i, _)| i);
            if atom.element.atomic_number() == 6 {
                assert!(matches!(
                    type_map[&idx.unwrap()],
                    UffType::C_3 | UffType::C_2
                ));
            }
        }
    }

    #[test]
    fn assign_types_benzene_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        let types = assign_uff_types(&mol);
        // All C aromatic → C_R
        for (_, t) in &types {
            assert_eq!(*t, UffType::C_R);
        }
    }

    #[test]
    fn energy_finite() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [2.5, 1.2, 0.0]];
        let e = uff_total_energy(&mol, &types, &coords);
        assert!(e.is_finite(), "energy should be finite: {e}");
    }

    #[test]
    fn fused_aromatic_uff_includes_planarity_terms() {
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let types = assign_uff_types(&mol);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        assert!(
            !prepared.inversions.is_empty(),
            "fused aromatic trigonal centers must have UFF inversion terms"
        );

        let mut folded = vec![[0.0, 0.0, 0.0]; mol.atom_count()];
        for (idx, point) in folded.iter_mut().enumerate() {
            point[0] = idx as f64;
            point[1] = (idx as f64 * 0.7).dsin();
            point[2] = if idx % 2 == 0 { 0.8 } else { -0.6 };
        }
        let mut planar = folded.clone();
        for point in &mut planar {
            point[2] = 0.0;
        }
        assert!(
            prepared.energy(&planar) < prepared.energy(&folded),
            "the aromatic inversion term must prefer a planar fused-ring geometry"
        );
    }

    #[test]
    fn butane_uff_includes_periodic_torsion_and_matches_energy_difference() {
        let mol = parse("CCCC").expect("butane");
        let types = assign_uff_types(&mol);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        assert_eq!(prepared.torsions.len(), 1);
        assert_eq!(prepared.torsions[0].periodicity, 3);
        let anti = [
            [0.0, 0.0, 0.0],
            [1.54, 0.0, 0.0],
            [2.95, 1.0, 0.0],
            [4.35, 1.0, 0.0],
        ];
        let gauche = [
            [0.0, 0.0, 0.0],
            [1.54, 0.0, 0.0],
            [2.95, 1.0, 0.0],
            [3.80, 1.85, 0.85],
        ];
        assert_ne!(
            PreparedUffEnergy::torsion_energy(prepared.torsions[0], &anti),
            PreparedUffEnergy::torsion_energy(prepared.torsions[0], &gauche)
        );
        let gradient = prepared.analytic_gradient(&gauche);
        let delta = 1e-5;
        for atom in 0..gauche.len() {
            for axis in 0..3 {
                let mut plus = gauche;
                let mut minus = gauche;
                plus[atom][axis] += delta;
                minus[atom][axis] -= delta;
                let expected = (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                assert!(
                    (gradient[atom][axis] - expected).abs() < 3e-3,
                    "torsion gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                    gradient[atom][axis]
                );
            }
        }
    }

    #[test]
    fn uff_torsion_selects_sp2_and_mixed_central_bond_rules() {
        let butadiene = parse("C=CC=C").expect("butadiene");
        let butadiene_types = assign_uff_types(&butadiene);
        let prepared = PreparedUffEnergy::new(&butadiene, &butadiene_types);
        assert_eq!(prepared.torsions.len(), 1);
        assert_eq!(prepared.torsions[0].periodicity, 2);
        assert!(prepared.torsions[0].barrier.is_finite());

        let allyl = parse("C=CCC").expect("allyl chain");
        let allyl_types = assign_uff_types(&allyl);
        let prepared = PreparedUffEnergy::new(&allyl, &allyl_types);
        assert_eq!(prepared.torsions.len(), 1);
        assert_eq!(prepared.torsions[0].periodicity, 6);
        assert_eq!(prepared.torsions[0].barrier, 1.0);
    }

    #[test]
    fn uff_torsion_does_not_depend_on_bond_endpoint_storage_order() {
        use chematic_core::{Atom, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(Atom::organic(Element::C)))
            .collect();
        for pair in [(1, 0), (2, 1), (3, 2)] {
            builder
                .add_bond(atoms[pair.0], atoms[pair.1], BondOrder::Single)
                .unwrap();
        }
        let mol = builder.build();
        let types = assign_uff_types(&mol);
        let prepared = PreparedUffEnergy::new(&mol, &types);

        assert_eq!(
            prepared.torsions.len(),
            1,
            "a C-C-C-C torsion must be prepared even when every bond is stored in reverse endpoint order"
        );
    }

    #[test]
    fn minimize_reduces_energy() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0], // stretched bond
            [3.5, 1.2, 0.0],
        ];
        let e0 = uff_total_energy(&mol, &types, &coords);
        let result = minimize_uff(&mol, &types, coords, 200);
        assert!(
            result.energy < e0,
            "minimisation should reduce energy: {e0} → {}",
            result.energy
        );
    }

    #[test]
    fn constrained_minimizer_rejects_constraint_violating_steps() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [3.5, 1.2, 0.0]];
        let initial = coords.clone();
        let result = minimize_uff_with_constraint(&mol, &types, coords, 20, |candidate| {
            candidate == initial.as_slice()
        });
        assert_eq!(result.coords, initial);
        assert_eq!(result.energy, uff_total_energy(&mol, &types, &initial));
        assert!(!result.converged);
        assert!(
            !result.rejected_unsound_step,
            "a caller constraint rejection must not be reported as an unsound UFF step"
        );
    }

    #[test]
    fn minimize_uff_reports_sound_on_ordinary_ethanol() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [2.5, 1.2, 0.0]];
        let result = minimize_uff(&mol, &types, coords, 200);
        assert!(
            result.sound,
            "an ordinary small molecule minimizing normally should report sound"
        );
        assert!(
            result.worst_bond_length <= MAX_SANE_UFF_BOND_LENGTH,
            "sound result must expose a bond length within the soundness limit"
        );
    }

    #[test]
    fn minimize_uff_reports_unsound_on_a_blown_up_bond() {
        // max_iter=0 returns the initial coords untouched (the `for iter in
        // 0..max_iter` loop never runs), so this deterministically exercises
        // `sound`'s bond-length check against a deliberately-stretched
        // C-C bond (5.0 Å, well past MAX_SANE_UFF_BOND_LENGTH) without
        // depending on steepest descent actually getting stuck there.
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [6.0, 1.2, 0.0]];
        let result = minimize_uff(&mol, &types, coords, 0);
        assert!(!result.sound, "a 5.0 Å C-C bond must be reported unsound");
        assert!(
            !result.converged,
            "an unsound result must never be reported as converged"
        );
        assert!(
            result.worst_bond_length > MAX_SANE_UFF_BOND_LENGTH,
            "unsound result must expose the stretched bond measurement"
        );
    }

    #[test]
    fn minimize_uff_reports_unsound_on_non_finite_coordinates() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0], [2.5, 1.2, 0.0]];
        let result = minimize_uff(&mol, &types, coords, 0);
        assert!(!result.sound, "non-finite coordinates must be unsound");
        assert!(
            result.worst_bond_length.is_infinite(),
            "non-finite geometry must fail closed in the exposed bond metric"
        );
    }

    /// Propane skeleton (C0-C1-C2, heavy atoms only — implicit H fills
    /// valence) placed so the C0-C1-C2 angle is exactly `theta_deg`, both
    /// C-C bonds at UFF's own C_3-C_3 r0 (~1.514 Å) so the bond-stretch term
    /// stays near zero and any energy blow-up as the angle closes is
    /// attributable to the angle and vdW terms only — the same isolation
    /// issue #176's own propane repro used.
    fn propane_at_angle(theta_deg: f64) -> (Molecule, Vec<(AtomIdx, UffType)>, Vec<[f64; 3]>) {
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let types = assign_uff_types(&mol);

        let r = 1.514_f64; // UFF C_3-C_3 bond length
        let half = theta_deg.to_radians() / 2.0;
        let coords = vec![
            [r * half.dcos(), r * half.dsin(), 0.0],
            [0.0, 0.0, 0.0],
            [r * half.dcos(), -r * half.dsin(), 0.0],
        ];
        (mol, types, coords)
    }

    #[test]
    fn uff_vdw_excludes_1_3_pair_propane_no_runaway() {
        // Issue #176's own measured pre-fix blow-up: 109.5°→16, 90°→110,
        // 70°→1326, 60°→6800 kcal/mol, driven entirely by the (wrongly
        // included) C0···C2 1-3 pair's LJ repulsion as it gets squeezed by
        // the closing angle. After excluding true 1-3 pairs, energy should
        // stay bounded (dominated by the smooth angle-bending term) even at
        // very closed angles.
        let mut energies = Vec::new();
        for &theta in &[109.5, 90.0, 70.0, 60.0, 45.0] {
            let (mol, types, coords) = propane_at_angle(theta);
            let e = uff_total_energy(&mol, &types, &coords);
            assert!(e.is_finite(), "energy at {theta}° should be finite: {e}");
            energies.push((theta, e));
        }
        for &(theta, e) in &energies {
            assert!(
                e < 200.0,
                "1-3 exclusion should prevent runaway vdW: at {theta}° energy={e} (issue #176 measured 6800 at 60° pre-fix)"
            );
        }
    }

    #[test]
    fn uff_vdw_still_repels_genuine_1_4_pair_butane() {
        // Positive control mirroring `mmff94_minimizer.rs`'s own
        // `vdw_more_repulsive_at_short_range`: a *real* 1-4 pair (butane's
        // terminal carbons) must still get full vdW repulsion — the 1-3 fix
        // must not over-exclude non-angle pairs.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        let mol = b.build();
        let types = assign_uff_types(&mol);

        let coords_close = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.6, 0.0, 0.0], // C3 forced very close to C0 (genuine 1-4)
        ];
        let coords_far = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [8.0, 0.0, 0.0],
        ];
        let e_close = uff_total_energy(&mol, &types, &coords_close);
        let e_far = uff_total_energy(&mol, &types, &coords_far);
        assert!(e_close.is_finite() && e_far.is_finite());
        assert!(
            e_close > e_far,
            "genuine 1-4 pair must still repel at short range: close={e_close} far={e_far}"
        );
    }

    #[test]
    fn uff_minimizer_uses_local_proposals_for_extreme_finite_vdw_gradients() {
        // Regression for the public 265-molecule 3D benchmark: three inputs
        // entered the UFF minimizer with a finite 1-4 contact just above the
        // 0.5 Å soft floor. The old raw `step * gradient` proposal jumped far
        // outside molecular scale, then exhausted the line-search floor
        // without ever trying a local move. This compact butane geometry
        // reproduces that gradient class without depending on the external
        // benchmark manifest.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(Atom::new(Element::C));
        let c1 = builder.add_atom(Atom::new(Element::C));
        let c2 = builder.add_atom(Atom::new(Element::C));
        let c3 = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(c0, c1, BondOrder::Single).unwrap();
        builder.add_bond(c1, c2, BondOrder::Single).unwrap();
        builder.add_bond(c2, c3, BondOrder::Single).unwrap();
        let mol = builder.build();
        let types = assign_uff_types(&mol);
        let initial = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.51, 0.02, 0.0],
        ];
        let initial_energy = uff_total_energy(&mol, &types, &initial);
        assert!(
            initial_energy > 1.0e8,
            "fixture must exercise the extreme-gradient basin"
        );

        let result = minimize_uff(&mol, &types, initial, 200);
        assert!(
            result.sound,
            "bounded proposals must retain a sound covalent geometry"
        );
        assert!(result.energy.is_finite());
        assert!(
            result.energy < initial_energy * 1.0e-4,
            "the local trust radius must escape the finite vdW clash: {} -> {}",
            initial_energy,
            result.energy
        );
    }

    #[test]
    fn uff_gradient_is_descent_direction_for_closed_angle_propane() {
        // At a strained, closed-angle propane geometry, stepping a small
        // distance along the analytic -gradient must lower the energy.
        let (mol, types, coords) = propane_at_angle(60.0);
        let e0 = uff_total_energy(&mol, &types, &coords);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        let grad = uff_gradient(&prepared, &coords);
        let grad_norm: f64 = grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(
            grad_norm > 1e-8,
            "gradient should be nonzero at strained geometry"
        );

        let step = 1e-4 / grad_norm;
        let stepped: Vec<[f64; 3]> = coords
            .iter()
            .zip(&grad)
            .map(|(c, g)| [c[0] - step * g[0], c[1] - step * g[1], c[2] - step * g[2]])
            .collect();
        let e1 = uff_total_energy(&mol, &types, &stepped);
        assert!(
            e1 < e0,
            "step along -gradient should decrease energy: {e0} → {e1}"
        );
    }

    #[test]
    fn analytic_uff_gradient_matches_finite_difference_reference() {
        let (mol, types, coords) = propane_at_angle(75.0);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        let actual = prepared.analytic_gradient(&coords);
        let delta = 1e-5;
        for atom in 0..coords.len() {
            for axis in 0..3 {
                let mut plus = coords.clone();
                let mut minus = coords.clone();
                plus[atom][axis] += delta;
                minus[atom][axis] -= delta;
                let expected = (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                assert!(
                    (actual[atom][axis] - expected).abs() < 2e-3,
                    "gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                    actual[atom][axis]
                );
            }
        }
    }

    #[test]
    fn analytic_uff_gradient_matches_vdw_14_finite_difference_reference() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        let mut builder = MoleculeBuilder::new();
        let atoms = [
            builder.add_atom(Atom::new(Element::C)),
            builder.add_atom(Atom::new(Element::C)),
            builder.add_atom(Atom::new(Element::C)),
            builder.add_atom(Atom::new(Element::C)),
        ];
        for pair in atoms.windows(2) {
            builder
                .add_bond(pair[0], pair[1], BondOrder::Single)
                .expect("butane probe bond should be valid");
        }
        let mol = builder.build();
        let types = assign_uff_types(&mol);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        let coords = [
            [0.0, 0.0, 0.0],
            [1.52, 0.1, 0.0],
            [3.01, 0.45, 0.2],
            [4.30, 0.85, 0.55],
        ];
        let actual = prepared.analytic_gradient(&coords);
        let delta = 1e-5;
        for atom in 0..coords.len() {
            for axis in 0..3 {
                let mut plus = coords;
                let mut minus = coords;
                plus[atom][axis] += delta;
                minus[atom][axis] -= delta;
                let expected = (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                assert!(
                    (actual[atom][axis] - expected).abs() < 2e-3,
                    "1-4 vdW gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                    actual[atom][axis]
                );
            }
        }
    }

    #[test]
    fn uff_handles_zinc_complex() {
        // Zinc as a metal centre — UFF should assign Zn type
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let zn = b.add_atom(Atom::new(Element::ZN));
        let n1 = b.add_atom(Atom::new(Element::N));
        let n2 = b.add_atom(Atom::new(Element::N));
        b.add_bond(zn, n1, BondOrder::Single).unwrap();
        b.add_bond(zn, n2, BondOrder::Single).unwrap();
        let mol = b.build();
        let types = assign_uff_types(&mol);
        let zn_type = types.iter().find(|(_, t)| *t == UffType::Zn);
        assert!(zn_type.is_some(), "Zn should get UffType::Zn");
    }

    #[test]
    fn uff_supported_metal_and_halogen_types_have_finite_prepared_energy() {
        // This is a local soundness probe for the explicit element-to-UFF
        // dispatch table. It does not claim reference-parameter parity or
        // chemical validity for arbitrary metal coordination geometries.
        use chematic_core::{Atom, Element, MoleculeBuilder};

        let elements = [
            Element::LI,
            Element::NA,
            Element::K,
            Element::MG,
            Element::CA,
            Element::AL,
            Element::SI,
            Element::FE,
            Element::CO,
            Element::NI,
            Element::CU,
            Element::ZN,
            Element::AG,
            Element::AU,
            Element::HG,
            Element::CL,
            Element::BR,
            Element::I,
        ];
        for element in elements {
            let mut builder = MoleculeBuilder::new();
            let first = builder.add_atom(Atom::new(element));
            let second = builder.add_atom(Atom::new(Element::C));
            builder
                .add_bond(first, second, BondOrder::Single)
                .expect("two-atom probe bond should be valid");
            let mol = builder.build();
            let types = assign_uff_types(&mol);
            let prepared = PreparedUffEnergy::new(&mol, &types);
            let coords = [[0.0, 0.0, 0.0], [1.8, 0.2, 0.1]];
            let energy = prepared.energy(&coords);
            let gradient = prepared.analytic_gradient(&coords);
            assert!(energy.is_finite(), "{element:?} energy must be finite");
            assert!(
                gradient.iter().flatten().all(|value| value.is_finite()),
                "{element:?} gradient must be finite: {gradient:?}"
            );
            let delta = 1e-5;
            for atom in 0..coords.len() {
                for axis in 0..3 {
                    let mut plus = coords;
                    let mut minus = coords;
                    plus[atom][axis] += delta;
                    minus[atom][axis] -= delta;
                    let expected =
                        (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                    assert!(
                        (gradient[atom][axis] - expected).abs() < 2e-3,
                        "{element:?} gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                        gradient[atom][axis]
                    );
                }
            }
        }
    }

    #[test]
    fn every_uff_type_has_finite_prepared_energy_and_gradient() {
        // Dispatch coverage above exercises elements that are commonly
        // reached through the public assigner. This table additionally keeps
        // every declared UffType parameter bundle executable, including
        // hybridization-specific and fallback-only variants.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

        let all_types = [
            UffType::C_3,
            UffType::C_2,
            UffType::C_1,
            UffType::C_R,
            UffType::N_3,
            UffType::N_2,
            UffType::N_1,
            UffType::N_R,
            UffType::O_3,
            UffType::O_2,
            UffType::O_1,
            UffType::O_R,
            UffType::S_3,
            UffType::S_2,
            UffType::S_R,
            UffType::P_3,
            UffType::P_R,
            UffType::H_,
            UffType::F_,
            UffType::Cl,
            UffType::Br,
            UffType::I_,
            UffType::Li,
            UffType::Na,
            UffType::K,
            UffType::Ca,
            UffType::Mg,
            UffType::Fe,
            UffType::Co,
            UffType::Ni,
            UffType::Cu,
            UffType::Zn,
            UffType::Mn,
            UffType::Cr,
            UffType::V_,
            UffType::Mo,
            UffType::W_,
            UffType::Pd,
            UffType::Pt,
            UffType::Au,
            UffType::Ag,
            UffType::Hg,
            UffType::Al,
            UffType::Si,
            UffType::Unknown,
        ];
        assert_eq!(all_types.len(), 45);

        let mut builder = MoleculeBuilder::new();
        let first = builder.add_atom(Atom::new(Element::C));
        let second = builder.add_atom(Atom::new(Element::C));
        builder
            .add_bond(first, second, BondOrder::Single)
            .expect("two-atom probe bond should be valid");
        let mol = builder.build();
        let coords = [[0.0, 0.0, 0.0], [1.8, 0.2, 0.1]];
        let delta = 1e-5;

        for uff_type in all_types {
            let types = vec![(first, uff_type), (second, UffType::C_3)];
            let prepared = PreparedUffEnergy::new(&mol, &types);
            let energy = prepared.energy(&coords);
            let gradient = prepared.analytic_gradient(&coords);
            assert!(energy.is_finite(), "{uff_type:?} energy must be finite");
            assert!(
                gradient.iter().flatten().all(|value| value.is_finite()),
                "{uff_type:?} gradient must be finite: {gradient:?}"
            );
            for atom in 0..coords.len() {
                for axis in 0..3 {
                    let mut plus = coords;
                    let mut minus = coords;
                    plus[atom][axis] += delta;
                    minus[atom][axis] -= delta;
                    let expected =
                        (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                    assert!(
                        (gradient[atom][axis] - expected).abs() < 2e-3,
                        "{uff_type:?} gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                        gradient[atom][axis]
                    );
                }
            }
        }
    }
}
