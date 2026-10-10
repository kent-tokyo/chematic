use crate::{AngleConstraint, BondConstraint, ConstraintSet, Coords3D, Point3};
use chematic_core::AtomIdx;

#[test]
fn constraint_diagnostics_and_projection_preserve_pair_centroids() {
    let mol = chematic_smiles::parse("CC").unwrap();
    let constraint = BondConstraint::new(AtomIdx(0), AtomIdx(1), 1.5);
    let constraints = ConstraintSet {
        bonds: vec![constraint.clone()],
        angles: vec![],
    };
    for (distance, violation) in [(0.5, 0.95), (1.5, 0.0), (2.5, 0.95)] {
        let coords = Coords3D {
            points: vec![Point3::zero(), Point3::new(distance, 0.0, 0.0)],
        };
        assert!((constraint.violation(&coords) - violation).abs() < 1e-12);
        assert!((constraints.max_violation(&coords) - violation).abs() < 1e-12);
        let unchanged = crate::satisfy_constraints(&coords, &mol, &constraints, 0);
        assert_eq!(unchanged, coords);
        let projected = crate::satisfy_constraints(&coords, &mol, &constraints, 8);
        assert!(constraints.max_violation(&projected) < 1e-12);
        assert!((projected.points[0].x + projected.points[1].x - distance).abs() < 1e-12);
    }
    let coincident = Coords3D::new_zeroed(2);
    assert_eq!(
        crate::satisfy_constraints(&coincident, &mol, &constraints, 2),
        coincident
    );
    let angle = AngleConstraint::new(
        AtomIdx(0),
        AtomIdx(1),
        AtomIdx(2),
        std::f64::consts::FRAC_PI_2,
    );
    for (degrees, violation) in [(30.0_f64, 55.0_f64), (90.0, 0.0), (150.0, 55.0)] {
        let a = degrees.to_radians();
        let coords = Coords3D {
            points: vec![
                Point3::new(a.cos(), a.sin(), 0.0),
                Point3::zero(),
                Point3::new(1.0, 0.0, 0.0),
            ],
        };
        assert!((angle.violation(&coords) - violation.to_radians()).abs() < 1e-12);
        assert!(
            (crate::mol_transforms::get_bond_angle_deg(
                &coords,
                AtomIdx(0),
                AtomIdx(1),
                AtomIdx(2)
            ) - degrees)
                .abs()
                < 1e-10
        );
    }
}

#[test]
fn bounded_geometry_minimization_is_finite_and_rigid_motion_equivariant() {
    use crate::minimize::{
        ForceField, MinimizeConfig, minimize_dreiding_with_config, minimize_with_config,
    };
    // Diverse bond orders and non-bonded pairs exercise the legacy public
    // geometry paths. A translation must never change the predicted motion.
    for source in [
        "CC",
        "C=C",
        "C#N",
        "CCCO",
        "CN=O",
        "CSSC",
        "CP",
        "C[SiH3]",
        "CC(F)(Cl)Br",
        "c1ccncc1",
        "[H][H]",
        "NN",
        "N=N",
        "N#N",
        "OO",
        "O=O",
        "SS",
        "P#P",
        "[He].[Ne]",
    ] {
        let mol = chematic_smiles::parse(source).unwrap();
        let coords = Coords3D {
            points: (0..mol.atom_count())
                .map(|i| {
                    Point3::new(
                        i as f64 * 1.8,
                        (i as f64).sin() * 0.8,
                        (i as f64 * 0.7).cos() * 0.6,
                    )
                })
                .collect(),
        };
        let shifted = Coords3D {
            points: coords
                .points
                .iter()
                .map(|p| Point3::new(p.x + 5.0, p.y - 2.0, p.z + 1.0))
                .collect(),
        };
        let config = MinimizeConfig {
            max_steps: 4,
            force_field: ForceField::UFF,
            ..Default::default()
        };
        for dreiding in [false, true] {
            let run = |c| {
                if dreiding {
                    minimize_dreiding_with_config(&mol, c, &config)
                } else {
                    minimize_with_config(&mol, c, &config)
                }
            };
            let result = run(coords.clone());
            let moved = run(shifted.clone());
            assert_eq!(result.atom_count(), mol.atom_count());
            assert!(result.is_finite(), "{source}, DREIDING={dreiding}");
            for (p, q) in result.points.iter().zip(&moved.points) {
                assert!(
                    (p.x + 5.0 - q.x).abs() < 1e-6
                        && (p.y - 2.0 - q.y).abs() < 1e-6
                        && (p.z + 1.0 - q.z).abs() < 1e-6,
                    "{source}, DREIDING={dreiding}"
                );
            }
        }
    }
}

#[test]
fn two_equal_mass_atoms_have_analytic_principal_moments() {
    let mol = chematic_smiles::parse("CC").unwrap();
    let coords = Coords3D {
        points: vec![Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
    };
    assert!(crate::shape_descriptors::pmi1(&mol, &coords).abs() < 1e-12);
    for value in [
        crate::shape_descriptors::pmi2(&mol, &coords),
        crate::shape_descriptors::pmi3(&mol, &coords),
    ] {
        assert!((value - 24.022).abs() < 1e-12);
    }
}

#[test]
fn double_bond_diagnostics_agree_with_the_verifier_for_each_geometric_side() {
    use crate::stereo_constraints::*;
    use chematic_core::BondOrder;
    for source in ["F/C=C/F", "F/C=C\\F"] {
        let mol = chematic_smiles::parse(source).unwrap();
        let (bond, _) = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .unwrap();
        for y in [1., -1.] {
            let coords = Coords3D {
                points: vec![
                    Point3::new(0., 1., 0.),
                    Point3::zero(),
                    Point3::new(1., 0., 0.),
                    Point3::new(1., y, 0.),
                ],
            };
            let info = debug_double_bond(&mol, &coords, bond).unwrap();
            assert_eq!(info.bond, bond);
            assert_eq!(
                info.status,
                verify_stereo(&mol, &coords)
                    .double_bond
                    .iter()
                    .find(|r| r.bond == bond)
                    .unwrap()
                    .status
            );
            assert!((info.actual_angle_deg - info.actual_angle_rad.to_degrees()).abs() < 1e-12);
            assert!((info.margin_from_boundary_deg - 90.).abs() < 1e-12);
            assert_eq!(debug_all_double_bonds(&mol, &coords).len(), 1);
            assert!(debug_double_bond(&mol, &coords, chematic_core::BondIdx(0)).is_none());
        }
        assert!(debug_double_bond(&mol, &Coords3D::new_zeroed(4), bond).is_none());
    }
}

#[test]
fn shape_descriptor_errors_keep_coordinate_counts_and_bad_atom_indices() {
    use crate::rdkit_shape_descriptors::*;
    let mol = chematic_smiles::parse("CCCC").unwrap();
    for function in [
        rdkit_pmi1,
        rdkit_pmi2,
        rdkit_pmi3,
        rdkit_npr1,
        rdkit_npr2,
        rdkit_radius_of_gyration,
        rdkit_inertial_shape_factor,
        rdkit_eccentricity,
        rdkit_asphericity,
        rdkit_spherocity_index,
        rdkit_pbf,
    ] {
        let error = function(&mol, &Coords3D::new_zeroed(3)).unwrap_err();
        assert_eq!(
            error,
            RdkitDescriptorError::AtomCoordCountMismatch {
                atoms: 4,
                coords: 3
            }
        );
        assert!(error.to_string().contains("(4)") && error.to_string().contains("(3)"));
        let mut coords = Coords3D::new_zeroed(4);
        coords.points[2].z = f64::NAN;
        let error = function(&mol, &coords).unwrap_err();
        assert_eq!(
            error,
            RdkitDescriptorError::NonFiniteCoordinate { atom_index: 2 }
        );
        assert!(error.to_string().contains("index 2"));
        let empty = chematic_core::MoleculeBuilder::new().build();
        let error = function(&empty, &Coords3D::new_zeroed(0)).unwrap_err();
        assert_eq!(error, RdkitDescriptorError::ZeroAtoms);
        assert!(error.to_string().contains("zero atoms"));
    }
    let one = chematic_smiles::parse("C").unwrap();
    let e = rdkit_pbf(&one, &Coords3D::new_zeroed(1)).unwrap_err();
    assert_eq!(
        e,
        RdkitDescriptorError::TooFewAtoms {
            required: 4,
            found: 1
        }
    );
    assert!(e.to_string().contains(">= 4") && e.to_string().contains("found 1"));
    let e = rdkit_npr1(&mol, &Coords3D::new_zeroed(4)).unwrap_err();
    assert!(matches!(e, RdkitDescriptorError::DegenerateGeometry { .. }));
    assert!(e.to_string().contains("degenerate"));
}

#[test]
fn legacy_etkdg_amide_torsions_are_near_a_planar_angle_after_generation() {
    for source in [
        "CC(=O)NC",
        "CC(=O)N(C)C",
        "CCC(=O)NCC",
        "O=C(NCC)C",
        "CC(=O)N(C)CC",
        "CC(=O)N(C)CCC",
    ] {
        let mol = chematic_smiles::parse(source).unwrap();
        for noise in [0., 30.] {
            let coords = crate::etkdg::generate_coords_etkdg_with_noise(&mol, noise);
            assert_eq!(coords.atom_count(), mol.atom_count());
            assert!(coords.is_finite());
            let (n, _) = mol
                .atoms()
                .find(|(_, a)| a.element.atomic_number() == 7)
                .unwrap();
            let c = mol
                .neighbors(n)
                .find(|(a, _)| {
                    mol.neighbors(*a).any(|(o, b)| {
                        mol.atom(o).element.atomic_number() == 8
                            && mol.bond(b).order == chematic_core::BondOrder::Double
                    })
                })
                .unwrap()
                .0;
            let a = mol.neighbors(n).find(|(a, _)| *a != c).unwrap().0;
            let d = mol.neighbors(c).find(|(a, _)| *a != n).unwrap().0;
            let omega = crate::mol_transforms::get_dihedral(&coords, a, n, c, d)
                .unwrap()
                .to_degrees()
                .abs();
            assert!(
                omega.min((omega - 180.).abs()) <= 30. + 1e-8,
                "{source} noise {noise}: {omega}"
            );
        }
    }
}

#[test]
fn strict_pdb_errors_name_the_bad_field_and_its_physical_line() {
    use crate::pdb::*;
    let base = "ATOM     17  C   LIG A  23       0.000   0.000   0.000  1.00  0.00           C  ";
    assert_eq!(
        parse_pdb_atoms_strict(base, &Default::default())
            .unwrap()
            .len(),
        1
    );
    for (start, end, field) in [
        (6, 11, "serial"),
        (22, 26, "residue sequence"),
        (30, 38, "x"),
        (38, 46, "y"),
        (46, 54, "z"),
    ] {
        for bad in ["bad", "", "NaN", "inf"] {
            let text = format!(
                "REMARK test\n{}{:>width$}{}",
                &base[..start],
                bad,
                &base[end..],
                width = end - start
            );
            let e = parse_pdb_atoms_strict(&text, &Default::default()).unwrap_err();
            let PdbStrictError::Parse(detail) = &e else {
                panic!("{e}")
            };
            assert_eq!(detail.line, 2);
            assert_eq!(detail.field, field);
            assert!(e.to_string().contains(field) && e.to_string().contains('2'));
        }
    }
    for which in 0..4 {
        let mut l = PdbParseLimits::default();
        match which {
            0 => l.max_input_bytes = 1,
            1 => l.max_line_bytes = 1,
            2 => l.max_atoms = 0,
            _ => l.max_models = 0,
        }
        let e = parse_pdb_atoms_strict(&format!("MODEL        1\n{base}\nENDMDL"), &l).unwrap_err();
        let PdbStrictError::Resource(r) = &e else {
            panic!("{e}")
        };
        assert!(r.actual > r.limit && e.to_string().contains(r.resource));
    }
}

#[test]
fn three_dimensional_xyz_rejections_report_headers_symbols_lines_and_limits() {
    use crate::xyz::*;
    for (text, fragment) in [
        ("bad\ncomment\n", "atom count"),
        ("1\ncomment\nC 0 0", "line 3"),
        ("1\ncomment\nUnknownX 0 0 0", "UnknownX"),
        ("2\ncomment\nC 0 0 0", "line 4"),
    ] {
        let e = parse_xyz(text).err().unwrap();
        assert!(e.to_string().contains(fragment), "{e}");
    }
    let text = "1\ncomment\nC 0 0 0";
    for which in 0..3 {
        let mut l = XyzParseLimits::default();
        match which {
            0 => l.max_input_bytes = 1,
            1 => l.max_atoms = 0,
            _ => l.max_line_bytes = 1,
        }
        let e = parse_xyz_with_limits(text, l).err().unwrap();
        assert!(e.to_string().contains("limit"));
    }
}
