use super::*;

#[test]
fn interruptible_uff_matches_ordinary_and_preserves_interrupted_coordinates() {
    use crate::rdkit_mmff::tests::{COORDS, paracetamol};
    let field = RdkitUffField::new(&paracetamol(), &COORDS, 10.0, true);
    for limit in [0, 1, 200] {
        let mut expected = COORDS.as_flattened().to_vec();
        let mut actual = expected.clone();
        let (status, energy) = field.optimize(&mut expected, limit);
        let result = field.optimize_interruptible(&mut actual, limit, || false);
        let RdkitUffInterruptibleOptimization::Finished {
            status: actual_status,
            energy: actual_energy,
            ..
        } = result
        else {
            panic!("unexpected cancellation")
        };
        assert_eq!((actual_status, actual_energy), (status, energy));
        assert_eq!(actual, expected);
    }
    let mut coords = COORDS.as_flattened().to_vec();
    let original = coords.clone();
    assert_eq!(
        field.optimize_interruptible(&mut coords, 200, || true),
        RdkitUffInterruptibleOptimization::Interrupted {
            completed_iterations: 0,
            energy_at_position: None
        }
    );
    assert_eq!(coords, original);
    let mut polls = 0;
    let result = field.optimize_interruptible(&mut coords, 200, || {
        polls += 1;
        polls == 20
    });
    let RdkitUffInterruptibleOptimization::Interrupted {
        completed_iterations,
        energy_at_position: Some(energy),
    } = result
    else {
        panic!("expected interruption with known energy: {result:?}")
    };
    assert!(completed_iterations > 0);
    assert!((energy - field.energy(&coords)).abs() < 1e-10);
}

/// RDKit 2026.03.1 `UFFGetMoleculeForceField(vdwThresh=10)`:
/// `CalcEnergy()` and `Minimize(maxIts=200)` on paracetamol.
#[test]
fn uff_matches_rdkit_with_platform_tolerance() {
    use crate::rdkit_mmff::tests::{COORDS, paracetamol};
    let mol = paracetamol();
    assert!(rdkit_uff_has_all_params(&mol));
    let ff = RdkitUffField::new(&mol, &COORDS, 10.0, true);
    let mut pos = COORDS.as_flattened().to_vec();
    assert!((ff.energy(&pos) - 40.395092089572906).abs() <= 1e-12);
    let (status, e) = ff.optimize(&mut pos, 200);
    assert_eq!(status, 0);
    assert!((e - 19.492409722746174).abs() <= 1e-9);
    assert!((pos[0] - 3.764669178348734).abs() <= 1e-6);
}

#[test]
fn int_pow_matches_rdkit_recursion() {
    let x = 1.2345678901234567_f64;
    assert_eq!(int_pow(x, 5), (x * x) * (x * x) * x);
    assert_eq!(int_pow(x, 6), ((x * x) * x) * ((x * x) * x));
    let six = ((x * x) * x) * ((x * x) * x);
    assert_eq!(int_pow(x, 13), six * six * x);
}

#[test]
fn param_table_is_sorted_and_has_common_types() {
    assert!(UFF_PARAMS.windows(2).all(|w| w[0].0 < w[1].0));
    for l in [
        "C_3", "C_R", "C_2", "C_1", "N_3", "N_R", "O_3", "O_R", "H_", "S_3+2", "Cl",
    ] {
        assert!(params_for(l).is_some(), "{l}");
    }
}
