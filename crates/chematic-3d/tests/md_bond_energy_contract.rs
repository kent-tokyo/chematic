use chematic_3d::md::{MDConfig, Thermostat, run_md};
use chematic_3d::{Coords3D, Point3};
use chematic_core::AtomIdx;
use chematic_smiles::parse;

#[test]
fn isolated_bonds_have_minimum_md_energy_at_the_model_distance() {
    let config = MDConfig {
        timestep_fs: 0.0,
        steps: 1,
        temperature_k: 0.0,
        thermostat: Thermostat::None,
        save_every: 1,
        coulomb: false,
    };
    for (smiles, equilibrium) in [
        ("[H]C", 1.090),
        ("C=N", 1.279),
        ("C#N", 1.158),
        ("C=O", 1.217),
        ("[H]N", 1.010),
        ("[H]O", 0.960),
        ("N=O", 1.240),
    ] {
        let molecule = parse(smiles).unwrap();
        for (distance, expected) in [(equilibrium, 0.0), (equilibrium + 0.1, 3.5)] {
            let mut coords = Coords3D::new_zeroed(2);
            coords.set(AtomIdx(1), Point3::new(distance, 0.0, 0.0));
            let trajectory = run_md(&molecule, coords, &config);
            assert_eq!(trajectory.frames.len(), 1);
            assert!(
                (trajectory.frames[0].potential_energy - expected).abs() < 1e-9,
                "{smiles} at {distance}"
            );
        }
    }
}
