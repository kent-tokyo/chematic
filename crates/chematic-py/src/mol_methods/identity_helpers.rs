use super::*;

use chematic_3d::rdkit_embed::{RdkitEmbedError, RdkitEmbedOptions};

pub(super) type PyRdkitMmffConformer = (i32, f64, Vec<[f64; 3]>);
pub(super) type PyRdkitMmffConformerDiagnostics = (i32, f64, u32, Option<f64>, Vec<[f64; 3]>);
pub(super) type PyBondLengthConstraint = (usize, usize, f64, f64);
pub(super) type PyBondAngleConstraint = (usize, usize, usize, f64, f64);
pub(super) type PyTorsionConstraint = (usize, usize, usize, usize, f64, f64);

pub(super) fn internal_coordinate_constraints(
    bond_lengths: Option<Vec<PyBondLengthConstraint>>,
    bond_angles: Option<Vec<PyBondAngleConstraint>>,
    torsions: Option<Vec<PyTorsionConstraint>>,
) -> chematic_3d::rdkit_embed::InternalCoordinateConstraints {
    use chematic_3d::rdkit_embed::{
        BondAngleConstraint, BondLengthConstraint, InternalCoordinateConstraints,
        TorsionAngleConstraint,
    };

    InternalCoordinateConstraints {
        bond_lengths: bond_lengths
            .unwrap_or_default()
            .into_iter()
            .map(|(atom1, atom2, target, tolerance)| {
                BondLengthConstraint::around(atom1, atom2, target, tolerance)
            })
            .collect(),
        bond_angles: bond_angles
            .unwrap_or_default()
            .into_iter()
            .map(
                |(atom1, center, atom2, target_degrees, tolerance_degrees)| {
                    BondAngleConstraint::around_degrees(
                        atom1,
                        center,
                        atom2,
                        target_degrees,
                        tolerance_degrees,
                    )
                },
            )
            .collect(),
        torsions: torsions
            .unwrap_or_default()
            .into_iter()
            .map(
                |(atom1, atom2, atom3, atom4, target_degrees, tolerance_degrees)| {
                    TorsionAngleConstraint::around_degrees(
                        atom1,
                        atom2,
                        atom3,
                        atom4,
                        target_degrees,
                        tolerance_degrees,
                    )
                },
            )
            .collect(),
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn rdkit_embed_options(
    random_seed: i32,
    max_iterations: u32,
    use_exp_torsion_angle_prefs: bool,
    use_basic_knowledge: bool,
    et_version: u32,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_macrocycle_14_config: bool,
    enforce_chirality: bool,
    ignore_smoothing_failures: bool,
) -> PyResult<RdkitEmbedOptions> {
    if !(1..=2).contains(&et_version) {
        return Err(PyValueError::new_err(
            "Only version 1 and 2 of the experimental torsion-angle preferences (ETversion) supported",
        ));
    }
    Ok(RdkitEmbedOptions {
        random_seed,
        max_iterations,
        use_exp_torsion_angle_prefs,
        use_basic_knowledge,
        et_version,
        use_small_ring_torsions,
        use_macrocycle_torsions,
        use_macrocycle14config: use_macrocycle_14_config,
        enforce_chirality,
        force_trans_amides: true,
        ignore_smoothing_failures,
    })
}

pub(super) fn map_rdkit_embed_error(error: RdkitEmbedError) -> PyErr {
    match error {
        RdkitEmbedError::Failed => pyo3::exceptions::PyRuntimeError::new_err(error.to_string()),
        _ => PyValueError::new_err(error.to_string()),
    }
}
