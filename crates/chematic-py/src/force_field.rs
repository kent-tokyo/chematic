//! Persistent force-field bindings.
//!
//! The older `Mol.rdkit_*_energy()` and `Mol.rdkit_*_optimize()` methods build
//! their force field for every call.  This wrapper keeps the prepared term
//! lists and coordinates together so iterative modelling workflows can reuse
//! the expensive preparation step without changing the existing APIs.

use crate::Mol;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;

#[derive(Clone, Debug)]
enum PreparedField {
    Mmff {
        field: chematic_ff::rdkit_mmff::RdkitMmffField,
        variant: &'static str,
    },
    Uff(chematic_ff::rdkit_uff::RdkitUffField),
}

impl PreparedField {
    fn num_atoms(&self) -> usize {
        match self {
            Self::Mmff { field, .. } => field.num_atoms(),
            Self::Uff(field) => field.num_atoms(),
        }
    }

    fn kind(&self) -> &'static str {
        match self {
            Self::Mmff { variant, .. } => variant,
            Self::Uff(_) => "UFF",
        }
    }

    fn energy(&self, positions: &[f64]) -> f64 {
        match self {
            Self::Mmff { field, .. } => field.energy(positions),
            Self::Uff(field) => field.energy(positions),
        }
    }

    fn gradient(&self, positions: &[f64]) -> Vec<f64> {
        let mut gradient = vec![0.0; positions.len()];
        match self {
            Self::Mmff { field, .. } => field.gradient(positions, &mut gradient),
            Self::Uff(field) => field.gradient(positions, &mut gradient),
        }
        gradient
    }

    fn optimize(&self, positions: &mut [f64], max_iters: u32) -> (i32, f64) {
        match self {
            Self::Mmff { field, .. } => field.optimize(positions, max_iters),
            Self::Uff(field) => field.optimize(positions, max_iters),
        }
    }
}

fn flatten_positions(coords: Vec<[f64; 3]>, atom_count: usize) -> PyResult<Vec<f64>> {
    if coords.len() != atom_count {
        return Err(PyValueError::new_err(format!(
            "expected {atom_count} coordinates, got {}",
            coords.len()
        )));
    }
    if coords.iter().flatten().any(|value| !value.is_finite()) {
        return Err(PyValueError::new_err("coordinates must be finite"));
    }
    Ok(coords.into_iter().flatten().collect())
}

fn nested_positions(flat: &[f64]) -> Vec<[f64; 3]> {
    (0..flat.len() / 3)
        .map(|index| {
            let offset = index * 3;
            [flat[offset], flat[offset + 1], flat[offset + 2]]
        })
        .collect()
}

fn mmff_variant(variant: &str) -> PyResult<chematic_ff::rdkit_mmff::MmffVariant> {
    match variant {
        "MMFF94" => Ok(chematic_ff::rdkit_mmff::MmffVariant::Mmff94),
        "MMFF94s" | "MMFF94S" => Ok(chematic_ff::rdkit_mmff::MmffVariant::Mmff94s),
        _ => Err(PyValueError::new_err(format!(
            "unsupported MMFF variant {variant:?} (expected \"MMFF94\" or \"MMFF94s\")"
        ))),
    }
}

/// A prepared MMFF94 or UFF force field with mutable coordinates.
///
/// Construct with :meth:`Mol.mmff_force_field` or :meth:`Mol.uff_force_field`.
/// Reusing the object avoids repeating atom typing and force-field term
/// construction for each energy, gradient, or minimization call.
#[pyclass(name = "MolecularForceField", skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMolecularForceField {
    field: PreparedField,
    positions: Vec<f64>,
}

#[pymethods]
impl PyMolecularForceField {
    /// Constant relative dielectric for MMFF94/MMFF94s (default 1.0).
    /// UFF has no electrostatic term and raises ``ValueError``.
    #[getter]
    fn dielectric_constant(&self) -> PyResult<f64> {
        match &self.field {
            PreparedField::Mmff { field, .. } => Ok(field.dielectric_constant()),
            PreparedField::Uff(_) => Err(PyValueError::new_err(
                "UFF has no dielectric electrostatic term",
            )),
        }
    }

    /// Change constant screening immediately, preserving coordinates and terms.
    /// Finite positive values are required; overflowing coefficients are refused.
    fn set_dielectric_constant_<'py>(
        mut slf: PyRefMut<'py, Self>,
        value: f64,
    ) -> PyResult<PyRefMut<'py, Self>> {
        match &mut slf.field {
            PreparedField::Mmff { field, .. } => field
                .set_dielectric_constant(value)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            PreparedField::Uff(_) => {
                return Err(PyValueError::new_err(
                    "UFF has no dielectric electrostatic term",
                ));
            }
        }
        Ok(slf)
    }

    /// Force-field family (``"MMFF94"`` or ``"UFF"``).
    #[getter]
    fn kind(&self) -> &'static str {
        self.field.kind()
    }

    /// Number of atoms represented by the field.
    #[getter]
    fn num_atoms(&self) -> usize {
        self.field.num_atoms()
    }

    /// Current coordinates in molecule atom order.
    #[getter]
    fn positions(&self) -> Vec<[f64; 3]> {
        nested_positions(&self.positions)
    }

    /// Replace all coordinates and return ``self``.
    fn set_positions_<'py>(
        mut slf: PyRefMut<'py, Self>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<PyRefMut<'py, Self>> {
        slf.positions = flatten_positions(coords, slf.field.num_atoms())?;
        Ok(slf)
    }

    /// Return one atom position.
    fn position(&self, atom_index: usize) -> PyResult<[f64; 3]> {
        if atom_index >= self.field.num_atoms() {
            return Err(PyIndexError::new_err(format!(
                "atom index {atom_index} out of range for {} atoms",
                self.field.num_atoms()
            )));
        }
        let offset = atom_index * 3;
        Ok([
            self.positions[offset],
            self.positions[offset + 1],
            self.positions[offset + 2],
        ])
    }

    /// Replace one atom position and return ``self``.
    fn set_position_<'py>(
        mut slf: PyRefMut<'py, Self>,
        atom_index: usize,
        position: [f64; 3],
    ) -> PyResult<PyRefMut<'py, Self>> {
        if atom_index >= slf.field.num_atoms() {
            return Err(PyIndexError::new_err(format!(
                "atom index {atom_index} out of range for {} atoms",
                slf.field.num_atoms()
            )));
        }
        if position.iter().any(|value| !value.is_finite()) {
            return Err(PyValueError::new_err("coordinates must be finite"));
        }
        let offset = atom_index * 3;
        slf.positions[offset..offset + 3].copy_from_slice(&position);
        Ok(slf)
    }

    /// Energy at the current coordinates in kcal/mol.
    fn energy(&self) -> f64 {
        self.field.energy(&self.positions)
    }

    /// Unscaled Cartesian gradient at the current coordinates.
    fn gradient(&self) -> Vec<[f64; 3]> {
        nested_positions(&self.field.gradient(&self.positions))
    }

    /// Return ``(energy, gradient)`` without rebuilding the field.
    fn energy_gradient(&self) -> (f64, Vec<[f64; 3]>) {
        (self.energy(), self.gradient())
    }

    /// Minimize in place and return ``(status, energy)``.
    ///
    /// Status follows RDKit: 0 means converged and 1 means more iterations
    /// are needed. A failed line search raises ``ValueError``.
    #[pyo3(signature = (max_iters = 200))]
    fn minimize_(&mut self, max_iters: u32) -> PyResult<(i32, f64)> {
        let (status, energy) = self.field.optimize(&mut self.positions, max_iters);
        if status < 0 {
            return Err(PyValueError::new_err("bad direction in linearSearch"));
        }
        Ok((status, energy))
    }
}

#[pymethods]
impl Mol {
    /// Build a reusable RDKit-compatible MMFF94 force field.
    ///
    /// ``coords`` must contain one position per explicit atom. Non-bonded
    /// pairs are selected once at construction, matching RDKit's force-field
    /// lifecycle; coordinate updates reuse that prepared term set.
    #[pyo3(signature = (coords, variant = "MMFF94", non_bonded_thresh = 100.0, ignore_interfrag_interactions = true))]
    fn mmff_force_field(
        &self,
        coords: Vec<[f64; 3]>,
        variant: &str,
        non_bonded_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<PyMolecularForceField> {
        if !non_bonded_thresh.is_finite() || non_bonded_thresh < 0.0 {
            return Err(PyValueError::new_err(
                "non_bonded_thresh must be finite and non-negative",
            ));
        }
        let positions = flatten_positions(coords, self.inner.atom_count())?;
        let nested = nested_positions(&positions);
        let parsed_variant = mmff_variant(variant)?;
        let variant_name = match parsed_variant {
            chematic_ff::rdkit_mmff::MmffVariant::Mmff94 => "MMFF94",
            chematic_ff::rdkit_mmff::MmffVariant::Mmff94s => "MMFF94s",
        };
        let field = chematic_ff::rdkit_mmff::RdkitMmffField::with_terms(
            &self.inner,
            &nested,
            parsed_variant,
            non_bonded_thresh,
            ignore_interfrag_interactions,
            [true; 7],
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(PyMolecularForceField {
            field: PreparedField::Mmff {
                field,
                variant: variant_name,
            },
            positions,
        })
    }

    /// Build a reusable RDKit-compatible UFF force field.
    #[pyo3(signature = (coords, vdw_thresh = 10.0, ignore_interfrag_interactions = true))]
    fn uff_force_field(
        &self,
        coords: Vec<[f64; 3]>,
        vdw_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<PyMolecularForceField> {
        if !vdw_thresh.is_finite() || vdw_thresh < 0.0 {
            return Err(PyValueError::new_err(
                "vdw_thresh must be finite and non-negative",
            ));
        }
        let positions = flatten_positions(coords, self.inner.atom_count())?;
        let nested = nested_positions(&positions);
        let field = chematic_ff::rdkit_uff::RdkitUffField::new(
            &self.inner,
            &nested,
            vdw_thresh,
            ignore_interfrag_interactions,
        );
        Ok(PyMolecularForceField {
            field: PreparedField::Uff(field),
            positions,
        })
    }
}
