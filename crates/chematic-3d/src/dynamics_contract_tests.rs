use crate::md::{MDConfig, Thermostat, run_md};
use crate::{Coords3D, Point3};
use chematic_core::MoleculeBuilder;
use chematic_smiles::parse;

#[test]
fn zero_temperature_dynamics_preserve_translation_and_frame_schedule() {
    // Short, deterministic trajectories exercise force computation for different
    // bond orders and elements without making a long-time stability claim.
    for smiles in [
        "CCO", "CCN", "CCS", "CP(C)C", "C#CC", "C=CC", "c1ccccc1", "CNO", "NN", "N=N", "N#N", "OO",
        "O=O", "CF", "CCl", "CBr", "CI", "[SiH4]",
    ] {
        let mol = parse(smiles).unwrap();
        let coords = Coords3D {
            points: (0..mol.atom_count())
                .map(|i| {
                    let t = i as f64 * 1.7;
                    Point3::new(2. * t.cos(), 2. * t.sin(), i as f64 * 0.9)
                })
                .collect(),
        };
        let translated = Coords3D {
            points: coords
                .points
                .iter()
                .map(|p| Point3::new(p.x + 3., p.y - 2., p.z + 1.))
                .collect(),
        };
        for thermostat in [Thermostat::None, Thermostat::Berendsen { tau_fs: 100. }] {
            let config = MDConfig {
                timestep_fs: 0.01,
                steps: 3,
                save_every: 1,
                temperature_k: 0.,
                thermostat,
                coulomb: true,
            };
            let first = run_md(&mol, coords.clone(), &config);
            let second = run_md(&mol, translated.clone(), &config);
            assert_eq!(first.frames.len(), 3);
            for (index, (a, b)) in first.frames.iter().zip(&second.frames).enumerate() {
                assert_eq!(a.step, index + 1);
                assert!(a.coords.is_finite());
                assert!(a.potential_energy.is_finite());
                assert!(a.kinetic_energy.is_finite() && a.kinetic_energy >= 0.);
                assert!(
                    (a.potential_energy - b.potential_energy).abs()
                        < 1e-5 * (1. + a.potential_energy.abs()),
                    "{smiles}"
                );
                for (p, q) in a.coords.points.iter().zip(&b.coords.points) {
                    assert!((q.x - p.x - 3.).abs() < 1e-6, "{smiles}");
                    assert!((q.y - p.y + 2.).abs() < 1e-6, "{smiles}");
                    assert!((q.z - p.z - 1.).abs() < 1e-6, "{smiles}");
                }
            }
        }
    }
}

#[test]
fn empty_and_zero_step_simulations_have_no_frames() {
    assert!(
        run_md(
            &MoleculeBuilder::new().build(),
            Coords3D::default(),
            &MDConfig::default()
        )
        .frames
        .is_empty()
    );
    let mol = parse("C").unwrap();
    let config = MDConfig {
        steps: 0,
        ..Default::default()
    };
    assert!(
        run_md(&mol, Coords3D::new_zeroed(1), &config)
            .frames
            .is_empty()
    );
}
