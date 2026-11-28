use crate::geometry_check::{
    GeometryCheckConfig, GeometryCheckError, MAX_GEOMETRY_CHECK_ATOMS,
    MAX_GEOMETRY_CHECK_COORDINATE, check_geometry, compare_nonbonded_clashes,
};

#[test]
fn clash_comparison_detects_replaced_pairs_when_counts_are_equal() {
    let molecule = chematic_smiles::parse("C.C.C").unwrap();
    let before = Coords3D {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.4, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
    };
    let after = Coords3D {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.4, 0.0, 0.0),
        ],
    };
    let comparison = compare_nonbonded_clashes(&molecule, &before, &after, 0.8).unwrap();
    assert_eq!((comparison.before_count, comparison.after_count), (1, 1));
    assert_eq!(comparison.new_clashes[0].atoms, [AtomIdx(0), AtomIdx(2)]);
    assert_eq!(
        comparison.resolved_clashes[0].atoms,
        [AtomIdx(0), AtomIdx(1)]
    );
    assert!(
        compare_nonbonded_clashes(&molecule, &before, &before, 0.8)
            .unwrap()
            .new_clashes
            .is_empty()
    );
    let boundary = Coords3D {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.8, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
    };
    assert_eq!(
        compare_nonbonded_clashes(&molecule, &boundary, &boundary, 0.8)
            .unwrap()
            .after_count,
        0
    );
    for cutoff in [0.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            compare_nonbonded_clashes(&molecule, &before, &after, cutoff).unwrap_err(),
            GeometryCheckError::ClashCutoff
        );
    }
}
use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Coords3D, Element, MoleculeBuilder, Point3,
};

#[test]
fn geometry_check_rejects_invalid_shapes_limits_and_numbers() {
    let variants = [
        GeometryCheckError::EmptyMolecule,
        GeometryCheckError::AtomLimit,
        GeometryCheckError::CoordinateCount {
            expected: 2,
            actual: 1,
        },
        GeometryCheckError::NonfiniteCoordinate { atom: AtomIdx(3) },
        GeometryCheckError::CoordinateRange { atom: AtomIdx(4) },
        GeometryCheckError::BondTolerance,
        GeometryCheckError::ClashCutoff,
    ];
    for error in variants {
        assert!(!error.to_string().is_empty());
        assert!(std::error::Error::source(&error).is_none());
    }

    let empty = MoleculeBuilder::new().build();
    assert_eq!(
        check_geometry(&empty, &Coords3D::default(), &Default::default()).unwrap_err(),
        GeometryCheckError::EmptyMolecule
    );

    let mut builder = MoleculeBuilder::new();
    for _ in 0..=MAX_GEOMETRY_CHECK_ATOMS {
        builder.add_atom(Atom::new(Element::C));
    }
    let too_large = builder.build();
    assert_eq!(
        check_geometry(
            &too_large,
            &Coords3D::new_zeroed(too_large.atom_count()),
            &Default::default(),
        )
        .unwrap_err(),
        GeometryCheckError::AtomLimit
    );

    let methane = chematic_smiles::parse("C").unwrap();
    assert!(matches!(
        check_geometry(&methane, &Coords3D::default(), &Default::default()),
        Err(GeometryCheckError::CoordinateCount { .. })
    ));
    for (point, expected) in [
        (
            Point3::new(f64::NAN, 0.0, 0.0),
            GeometryCheckError::NonfiniteCoordinate { atom: AtomIdx(0) },
        ),
        (
            Point3::new(MAX_GEOMETRY_CHECK_COORDINATE + 1.0, 0.0, 0.0),
            GeometryCheckError::CoordinateRange { atom: AtomIdx(0) },
        ),
    ] {
        let coords = Coords3D {
            points: vec![point],
        };
        assert_eq!(
            check_geometry(&methane, &coords, &Default::default()).unwrap_err(),
            expected
        );
    }
    let coords = Coords3D::new_zeroed(1);
    for config in [
        GeometryCheckConfig {
            bond_relative_tolerance: 1.0,
            ..Default::default()
        },
        GeometryCheckConfig {
            nonbonded_min_distance: 0.0,
            ..Default::default()
        },
    ] {
        assert!(check_geometry(&methane, &coords, &config).is_err());
    }
}

#[test]
fn geometry_check_reports_skips_stretched_bonds_and_clashes() {
    let mut molecule = chematic_smiles::parse("CCC.C").unwrap();
    molecule.set_bond_order(BondIdx(1), BondOrder::QueryAny);
    let coords = Coords3D {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.1, 0.0, 0.0),
            Point3::new(0.1, 0.0, 0.0),
        ],
    };
    let report = check_geometry(&molecule, &coords, &Default::default()).unwrap();
    assert_eq!(report.atom_count, 4);
    assert_eq!(report.bonds_checked, 1);
    assert_eq!(report.skipped_bonds, vec![BondIdx(1)]);
    assert_eq!(report.bond_violations.len(), 1);
    assert!(report.max_relative_bond_deviation.unwrap() > 1.0);
    assert!(
        report
            .nonbonded_clashes
            .iter()
            .any(|clash| clash.atoms == [AtomIdx(0), AtomIdx(3)])
    );
}
