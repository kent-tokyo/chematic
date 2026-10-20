//! Original analytic controls for Open Babel #2651; no upstream molecule corpus.
use chematic_3d::{
    Coords3D, GeometryCheckConfig, GeometryCheckError as Error, Point3, check_geometry,
};
use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

fn coords(rows: &[[f64; 3]]) -> Coords3D {
    Coords3D {
        points: rows.iter().map(|p| Point3::new(p[0], p[1], p[2])).collect(),
    }
}

fn screen(smiles: &str, rows: &[[f64; 3]]) -> chematic_3d::GeometryCheckReport {
    let mol = chematic_smiles::parse(smiles).unwrap();
    check_geometry(&mol, &coords(rows), &GeometryCheckConfig::default()).unwrap()
}

#[test]
fn finite_but_collapsed_or_stretched_bonds_are_reported() {
    for length in [0., 0.1, 10.] {
        let r = screen("CC", &[[0.; 3], [length, 0., 0.]]);
        assert_eq!(r.bonds_checked, 1);
        assert_eq!(r.bond_violations.len(), 1);
        assert!(r.nonbonded_clashes.is_empty());
        let v = &r.bond_violations[0];
        assert_eq!(v.atoms, [AtomIdx(0), AtomIdx(1)]);
        assert_eq!(v.distance, length);
        assert!((v.reference_length - 1.54).abs() < 1e-6);
        assert_eq!(Some(v.relative_deviation), r.max_relative_bond_deviation);
    }
    assert!(
        screen("CC", &[[0.; 3], [1.54, 0., 0.]])
            .bond_violations
            .is_empty()
    );
}

#[test]
fn tiny_finite_distances_do_not_underflow_into_false_coincidences() {
    let r = screen("CC", &[[0.; 3], [1e-200, 0., 0.]]);
    assert_eq!(r.bond_violations[0].distance, 1e-200);
    let mol = chematic_smiles::parse("C.C").unwrap();
    let config = GeometryCheckConfig {
        nonbonded_min_distance: 1e-250,
        ..GeometryCheckConfig::default()
    };
    let r = check_geometry(&mol, &coords(&[[0.; 3], [1e-200, 0., 0.]]), &config).unwrap();
    assert!(r.nonbonded_clashes.is_empty());
}

#[test]
fn covalent_order_scales_are_not_flat_single_bond_defaults() {
    for (s, expected) in [
        ("CC", 1.54),
        ("C=C", 1.54 * 0.87),
        ("C#C", 1.54 * 0.78),
        ("c1ccccc1", 1.54 * 0.93),
    ] {
        let mol = chematic_smiles::parse(s).unwrap();
        let r = check_geometry(
            &mol,
            &Coords3D::new_zeroed(mol.atom_count()),
            &GeometryCheckConfig::default(),
        )
        .unwrap();
        assert_eq!(r.bonds_checked, mol.bond_count());
        assert_eq!(r.bond_violations.len(), mol.bond_count());
        for v in &r.bond_violations {
            assert!((v.reference_length - expected).abs() < 1e-6, "{s}: {v:?}");
        }
    }
}

#[test]
fn clashes_include_one_three_pairs_and_fragments_in_original_index_order() {
    let r = screen("CCC.C", &[[0.; 3], [1.54, 0., 0.], [0.; 3], [0.; 3]]);
    let pairs: Vec<_> = r
        .nonbonded_clashes
        .iter()
        .map(|v| v.atoms.map(|a| a.0))
        .collect();
    assert_eq!(pairs, [[0, 2], [0, 3], [2, 3]]);
    assert!(r.bond_violations.is_empty());
    assert_eq!(r.nonbonded_clashes[0].distance, 0.);
}

#[test]
fn threshold_equalities_and_configuration_are_explicit() {
    let mol = chematic_smiles::parse("CC.C").unwrap();
    let c = coords(&[[0.; 3], [1.54, 0., 0.], [0., 1.2, 0.]]);
    let r = check_geometry(&mol, &c, &GeometryCheckConfig::default()).unwrap();
    assert!(r.nonbonded_clashes.is_empty());
    let config = GeometryCheckConfig {
        bond_relative_tolerance: 0.,
        nonbonded_min_distance: 2.,
    };
    let strict = check_geometry(&mol, &c, &config).unwrap();
    assert_eq!(strict.config, config);
    assert_eq!(strict.nonbonded_clashes.len(), 2);
    // Use the actual f32-backed radius reference to avoid decimal-rounding ties.
    let reference = strict.bond_violations[0].reference_length;
    let tie = coords(&[[0.; 3], [reference * 1.5, 0., 0.], [20., 0., 0.]]);
    let tie = check_geometry(&mol, &tie, &GeometryCheckConfig::default()).unwrap();
    assert!(tie.bond_violations.is_empty());
    assert_eq!(tie.max_relative_bond_deviation, Some(0.5));
}

#[test]
fn unsupported_bonds_and_wildcards_are_not_successful_length_checks() {
    for order in [
        BondOrder::Dative,
        BondOrder::Zero,
        BondOrder::Quadruple,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(Atom::new(Element::C));
        let c = b.add_atom(Atom::new(Element::N));
        b.add_bond(a, c, order).unwrap();
        let r = check_geometry(
            &b.build(),
            &coords(&[[0.; 3], [0.; 3]]),
            &GeometryCheckConfig::default(),
        )
        .unwrap();
        assert_eq!(r.bonds_checked, 0);
        assert_eq!(r.skipped_bonds.len(), 1);
        assert_eq!(r.max_relative_bond_deviation, None);
        assert!(r.nonbonded_clashes.is_empty());
    }
    let r = screen("*C", &[[0.; 3], [0.; 3]]);
    assert_eq!(r.skipped_bonds.len(), 1);
    assert_eq!(r.bonds_checked, 0);
}

#[test]
fn mirror_and_planar_tetrahedral_outcomes_remain_distinct() {
    let s = "F[C@](Cl)(Br)I";
    let rows = [
        [1., 1., 1.],
        [0.; 3],
        [1., -1., -1.],
        [-1., 1., -1.],
        [-1., -1., 1.],
    ];
    let a = screen(s, &rows);
    let mirror = rows.map(|p| [-p[0], p[1], p[2]]);
    let b = screen(s, &mirror);
    assert_eq!(a.stereo.n_declared(), 1);
    assert_eq!(a.stereo.n_unevaluable(), 0);
    assert_eq!(a.stereo.n_violations() + b.stereo.n_violations(), 1);
    assert_eq!(a.stereo.n_satisfied() + b.stereo.n_satisfied(), 1);
    let planar = screen(s, &rows.map(|p| [p[0], p[1], 0.]));
    assert_eq!(planar.stereo.n_unevaluable(), 1);
    assert_eq!(planar.stereo.n_violations(), 0);
    let unspecified = screen("FC(Cl)(Br)I", &rows);
    assert_eq!(unspecified.stereo.n_declared(), 0);
}

#[test]
fn ez_and_unsupported_coordination_keep_per_element_results() {
    let rows = [[-1., 1., 0.], [0.; 3], [1.34, 0., 0.], [2.34, -1., 0.]];
    let e = screen("C/C=C/C", &rows);
    let z = screen("C/C=C\\C", &rows);
    assert_eq!(e.stereo.n_declared(), 1);
    assert_eq!(e.stereo.n_violations() + z.stereo.n_violations(), 1);
    let sp = screen("Cl[Pt@SP1](Cl)(N)N", &[[0.; 3]; 5]);
    assert_eq!(sp.stereo.n_unevaluable(), 1);
}

#[test]
fn malformed_inputs_return_errors_before_any_indexed_work() {
    let mol = chematic_smiles::parse("CC").unwrap();
    let config = GeometryCheckConfig::default();
    for count in [0, 1, 3] {
        assert!(matches!(
            check_geometry(&mol, &Coords3D::new_zeroed(count), &config),
            Err(Error::CoordinateCount { .. })
        ));
    }
    assert!(matches!(
        check_geometry(&MoleculeBuilder::new().build(), &coords(&[]), &config),
        Err(Error::EmptyMolecule)
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            check_geometry(&mol, &coords(&[[0.; 3], [value, 0., 0.]]), &config),
            Err(Error::NonfiniteCoordinate { atom: AtomIdx(1) })
        ));
    }
    for value in [1_000_001., -1_000_001., f64::MAX] {
        assert!(matches!(
            check_geometry(&mol, &coords(&[[0.; 3], [value, 0., 0.]]), &config),
            Err(Error::CoordinateRange { atom: AtomIdx(1) })
        ));
    }
    for tolerance in [-0.1, 1., f64::NAN, f64::INFINITY] {
        let c = GeometryCheckConfig {
            bond_relative_tolerance: tolerance,
            ..config
        };
        assert!(matches!(
            check_geometry(&mol, &Coords3D::new_zeroed(2), &c),
            Err(Error::BondTolerance)
        ));
    }
    for cutoff in [-1., 0., f64::NAN, f64::INFINITY] {
        let c = GeometryCheckConfig {
            nonbonded_min_distance: cutoff,
            ..config
        };
        assert!(matches!(
            check_geometry(&mol, &Coords3D::new_zeroed(2), &c),
            Err(Error::ClashCutoff)
        ));
    }
}

#[test]
fn atom_limit_and_single_atom_semantics() {
    let mut b = MoleculeBuilder::new();
    for _ in 0..2049 {
        b.add_atom(Atom::new(Element::C));
    }
    assert!(matches!(
        check_geometry(
            &b.build(),
            &Coords3D::new_zeroed(2049),
            &GeometryCheckConfig::default()
        ),
        Err(Error::AtomLimit)
    ));
    let r = screen("[He]", &[[0.; 3]]);
    assert_eq!(r.bonds_checked, 0);
    assert_eq!(r.max_relative_bond_deviation, None);
    assert!(r.nonbonded_clashes.is_empty());
}

#[test]
fn maximum_atom_count_and_coordinate_boundary_are_supported() {
    let mut b = MoleculeBuilder::new();
    for _ in 0..2048 {
        b.add_atom(Atom::new(Element::C));
    }
    let c = Coords3D {
        points: (0..2048)
            .map(|i| Point3::new(i as f64 * 2., 1_000_000., -1_000_000.))
            .collect(),
    };
    let r = check_geometry(&b.build(), &c, &GeometryCheckConfig::default()).unwrap();
    assert_eq!(r.atom_count, 2048);
    assert!(r.nonbonded_clashes.is_empty());
}

#[test]
fn enhanced_group_semantics_are_reported_as_unchecked() {
    let mut mol = chematic_smiles::parse("F[C@](Cl)(Br)I").unwrap();
    mol.add_stereo_group(chematic_core::StereoGroup::new(
        chematic_core::StereoGroupKind::Or(1),
        vec![AtomIdx(1)],
    ));
    let r = check_geometry(
        &mol,
        &Coords3D::new_zeroed(5),
        &GeometryCheckConfig::default(),
    )
    .unwrap();
    assert_eq!(r.stereo_groups_not_evaluated, 1);
    assert_eq!(mol.stereo_groups().len(), 1);
}

fn snapshot(mol: &Molecule) -> String {
    format!(
        "{:?} {:?} {:?} {:?}",
        mol.atoms().collect::<Vec<_>>(),
        mol.bonds().collect::<Vec<_>>(),
        mol.stereo_groups(),
        (0..mol.atom_count())
            .map(|i| mol.stereo_neighbor_order(AtomIdx(i as u32)))
            .collect::<Vec<_>>()
    )
}

#[test]
fn read_only_and_rigid_transform_invariance() {
    let mol = chematic_smiles::parse("F[C@](Cl)(Br)I").unwrap();
    let c = coords(&[
        [1., 1., 1.],
        [0.; 3],
        [1., -1., -1.],
        [-1., 1., -1.],
        [-1., -1., 1.],
    ]);
    let before_mol = snapshot(&mol);
    let before_coords = c.clone();
    let a = check_geometry(&mol, &c, &GeometryCheckConfig::default()).unwrap();
    let moved = Coords3D {
        points: c
            .points
            .iter()
            .map(|p| Point3::new(-p.y + 20., p.x - 10., p.z + 5.))
            .collect(),
    };
    let b = check_geometry(&mol, &moved, &GeometryCheckConfig::default()).unwrap();
    assert_eq!(a.bond_violations, b.bond_violations);
    assert_eq!(a.nonbonded_clashes, b.nonbonded_clashes);
    assert_eq!(a.stereo.tetrahedral, b.stereo.tetrahedral);
    assert_eq!(snapshot(&mol), before_mol);
    assert_eq!(c, before_coords);
}
