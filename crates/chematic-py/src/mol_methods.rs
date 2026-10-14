//! The `#[pymethods] impl Mol` block — all instance methods/properties on `Mol`.
//!
//! Kept as a single `impl` block (not further split) because PyO3 0.29 requires
//! the `multiple-pymethods` Cargo feature (untested in this crate) to spread
//! `#[pymethods]` for one `#[pyclass]` across more than one block.

use crate::EcfpBitInfo;
use crate::Mol;
use crate::RdkitBestAlignment;
use crate::RdkitMorganDetail;
use crate::formats::{bitvec2048_to_bytes, flat_to_coords3d};
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::sync::Arc;

type PyUffTerm = (String, Vec<usize>, Vec<f64>);
type PyExplicitHydrogenMmffResult = (Vec<[f64; 3]>, Vec<[f64; 3]>, f64, bool, usize, f64);

mod analysis;
mod descriptors;
mod extended;
mod fingerprints;
mod identity;
mod properties;
mod three_d;
mod transformations;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Maps the plain `radius`/`nbits` integers `rdkit_ecfp_config`/`rdkit_ecfp_config_detail`
/// accept from Python into chematic-fp's closed `RdkitMorganConfig` enums. An
/// unsupported value raises `ValueError` explicitly rather than being coerced to the
/// nearest supported one -- there is no guessed conversion for an option this API
/// doesn't (yet) support.
fn python_rdkit_morgan_config(
    radius: u32,
    nbits: usize,
    include_chirality: bool,
) -> PyResult<chematic_fp::RdkitMorganConfig> {
    let radius = match radius {
        0 => chematic_fp::RdkitMorganRadius::R0,
        1 => chematic_fp::RdkitMorganRadius::R1,
        2 => chematic_fp::RdkitMorganRadius::R2,
        3 => chematic_fp::RdkitMorganRadius::R3,
        other => {
            return Err(PyValueError::new_err(format!(
                "unsupported radius {other} -- rdkit_ecfp_config only supports 0, 1, 2, or 3 \
                 (each independently verified against a live RDKit oracle)"
            )));
        }
    };
    let fp_size = match nbits {
        128 => chematic_fp::RdkitMorganFpSize::B128,
        256 => chematic_fp::RdkitMorganFpSize::B256,
        512 => chematic_fp::RdkitMorganFpSize::B512,
        1024 => chematic_fp::RdkitMorganFpSize::B1024,
        2048 => chematic_fp::RdkitMorganFpSize::B2048,
        other => {
            return Err(PyValueError::new_err(format!(
                "unsupported nbits {other} -- rdkit_ecfp_config only supports 128, 256, 512, \
                 1024, or 2048"
            )));
        }
    };
    Ok(chematic_fp::RdkitMorganConfig {
        radius,
        fp_size,
        include_chirality,
    })
}

/// Bit-pack a variable-width `BitVecN` into `bit_width()/8` bytes, LSB-first --
/// generalizes `formats::bitvec2048_to_bytes` to `rdkit_ecfp_config`'s caller-chosen
/// widths (always a multiple of 8: 128/256/512/1024/2048).
fn bitvecn_to_bytes(fp: &chematic_fp::bitvec::BitVecN) -> Vec<u8> {
    fp.to_le_bytes()
}

/// Expand an LSB-first packed fingerprint into the dense ``uint8`` shape used
/// by NumPy-based ML code.  Keeping this conversion in one place prevents the
/// individual fingerprint bindings from drifting on bit order.
fn packed_fingerprint_to_numpy<'py>(
    py: Python<'py>,
    packed: &[u8],
    n_bits: usize,
) -> Bound<'py, PyArray1<u8>> {
    debug_assert!(n_bits <= packed.len() * 8);
    let bits = (0..n_bits)
        .map(|bit| (packed[bit / 8] >> (bit % 8)) & 1)
        .collect::<Vec<_>>();
    Array1::from_vec(bits).into_pyarray(py)
}

pyo3::create_exception!(
    chematic,
    StereoLossWarning,
    pyo3::exceptions::PyUserWarning,
    "A MOL block written without part of the molecule's stereo (a centre \
     left unwedged, a declared E/Z bond written \"either\"). The block is \
     still returned; pass strict=True to raise instead, or call \
     Mol.to_mol_block_with_report() for the lost atoms and bonds."
);

/// The block, with a [`StereoLossWarning`] naming what it loses, or a
/// `ValueError` instead when `strict`.
fn strict_mol_block(
    py: Python<'_>,
    block: String,
    loss: &chematic_mol::MolStereoLoss,
    strict: bool,
) -> PyResult<String> {
    if loss.is_empty() {
        return Ok(block);
    }
    if strict {
        return Err(PyValueError::new_err(loss.to_string()));
    }
    PyErr::warn(
        py,
        py.get_type::<StereoLossWarning>().as_any(),
        &std::ffi::CString::new(loss.to_string()).unwrap_or_default(),
        1,
    )?;
    Ok(block)
}

fn stereo_loss_dict<'py>(
    py: Python<'py>,
    loss: &chematic_mol::MolStereoLoss,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let atoms = |v: &[chematic_core::AtomIdx]| v.iter().map(|a| a.0).collect::<Vec<u32>>();
    d.set_item("centres", atoms(&loss.centres))?;
    d.set_item(
        "double_bonds",
        loss.double_bonds.iter().map(|b| b.0).collect::<Vec<u32>>(),
    )?;
    d.set_item(
        "non_tetrahedral_centres",
        atoms(&loss.non_tetrahedral_centres),
    )?;
    d.set_item("stereo_groups_dropped", loss.stereo_groups_dropped)?;
    Ok(d)
}

fn rdkit_uff_field(
    mol: &chematic_core::Molecule,
    coords: &[[f64; 3]],
    vdw_thresh: f64,
    ignore_interfrag: bool,
) -> PyResult<chematic_ff::rdkit_uff::RdkitUffField> {
    if coords.len() != mol.atom_count() {
        return Err(PyValueError::new_err(format!(
            "expected {} coordinates, got {}",
            mol.atom_count(),
            coords.len()
        )));
    }
    Ok(chematic_ff::rdkit_uff::RdkitUffField::new(
        mol,
        coords,
        vdw_thresh,
        ignore_interfrag,
    ))
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

fn rdkit_mmff_field(
    mol: &chematic_core::Molecule,
    coords: &[[f64; 3]],
    variant: chematic_ff::rdkit_mmff::MmffVariant,
    non_bonded_thresh: f64,
    ignore_interfrag: bool,
    terms: [bool; 7],
) -> PyResult<chematic_ff::rdkit_mmff::RdkitMmffField> {
    if coords.len() != mol.atom_count() {
        return Err(PyValueError::new_err(format!(
            "expected {} coordinates, got {}",
            mol.atom_count(),
            coords.len()
        )));
    }
    chematic_ff::rdkit_mmff::RdkitMmffField::with_terms(
        mol,
        coords,
        variant,
        non_bonded_thresh,
        ignore_interfrag,
        terms,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))
}
