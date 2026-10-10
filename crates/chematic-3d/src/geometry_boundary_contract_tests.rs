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
