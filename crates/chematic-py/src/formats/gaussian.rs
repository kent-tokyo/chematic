//! Cartesian Gaussian geometry without hiding coordinates or job charge/spin.
use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;
mod writer;
pub(super) use writer::write_gjf_cartesian;

/// Read a Cartesian GJF geometry, preserving input coordinate values and freeze codes.
///
/// Returns a dict with ``elements``, ``mol``, ``coords``, ``freeze_codes``, ``charge`` and
/// ``multiplicity``. Coordinates retain input numeric units: ``Units=Bohr``
/// is not converted. No bonds or atom formal charges are inferred. This is
/// not a complete job reader: route/title/basis sets are not retained.
/// Z-matrices, ONIOM atom rows, fragment charge/spin pairs, malformed rows
/// and multiple Link1 jobs raise ValueError. Freeze codes are signed integers;
/// absent codes are None, and their optimization semantics are not interpreted.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=16_777_216, max_atoms=100_000))]
pub(super) fn parse_gjf_geometry<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: usize,
    max_atoms: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let geometry = chematic_mol::parse_gjf_geometry_with_limits(
        text,
        &chematic_mol::GaussianParseLimits {
            max_input_bytes,
            max_atoms,
            ..Default::default()
        },
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let result = PyDict::new(py);
    result.set_item(
        "elements",
        geometry
            .mol
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
    )?;
    result.set_item(
        "mol",
        Mol {
            inner: Arc::new(geometry.mol),
            props: Default::default(),
            num_rings_cache: Default::default(),
        },
    )?;
    let coords: Vec<[f64; 3]> = geometry
        .coords
        .into_iter()
        .map(|(x, y, z)| [x, y, z])
        .collect();
    result.set_item("coords", coords)?;
    result.set_item("freeze_codes", geometry.freeze_codes)?;
    result.set_item("charge", geometry.charge)?;
    result.set_item("multiplicity", geometry.multiplicity)?;
    Ok(result)
}

/// Retain a bounded single Cartesian GJF job as exact UTF-8 source.
/// Opaque job sections are preserved, not interpreted. No external files are read.
/// This is source retention, not job execution or complete Gaussian validation.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=16_777_216, max_atoms=100_000))]
pub(super) fn roundtrip_gjf_record(
    text: &str,
    max_input_bytes: usize,
    max_atoms: usize,
) -> PyResult<String> {
    let record = chematic_mol::parse_gjf_source_record_with_limits(
        text,
        &chematic_mol::GaussianParseLimits {
            max_input_bytes,
            max_atoms,
            ..Default::default()
        },
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(chematic_mol::write_gjf_source_record(&record))
}

fn gjf_document(
    text: &str,
    max_input_bytes: usize,
    max_atoms: usize,
    max_jobs: usize,
) -> PyResult<chematic_mol::GjfSourceDocument> {
    chematic_mol::parse_gjf_source_document_with_limits(
        text,
        &chematic_mol::GjfDocumentLimits {
            gaussian: chematic_mol::GaussianParseLimits {
                max_input_bytes,
                max_atoms,
                ..Default::default()
            },
            max_jobs,
        },
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Retain every explicit Cartesian Link1 job; never select the first or last.
/// Atom/byte budgets cover the whole document. Opaque sections are not validated.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=16_777_216, max_atoms=100_000, max_jobs=128))]
pub(super) fn roundtrip_gjf_document(
    text: &str,
    max_input_bytes: usize,
    max_atoms: usize,
    max_jobs: usize,
) -> PyResult<String> {
    gjf_document(text, max_input_bytes, max_atoms, max_jobs)
        .map(|document| chematic_mol::write_gjf_source_document(&document))
}

/// Inspect all explicit Cartesian jobs with ordered geometry and UTF-8 source ranges.
/// No bonds, coordinate-unit conversion or external resources are inferred.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=16_777_216, max_atoms=100_000, max_jobs=128))]
pub(super) fn parse_gjf_document_json(
    text: &str,
    max_input_bytes: usize,
    max_atoms: usize,
    max_jobs: usize,
) -> PyResult<String> {
    let json = gjf_document(text, max_input_bytes, max_atoms, max_jobs)?
        .to_json()
        .to_string();
    if json.len() > 16_777_216 {
        return Err(PyValueError::new_err(
            "GJF document JSON exceeds maximum output size",
        ));
    }
    Ok(json)
}
