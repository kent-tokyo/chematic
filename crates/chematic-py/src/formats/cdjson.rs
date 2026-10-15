use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyListMethods};

/// Parse a checked molecular ChemDoodle JSON root or content document.
/// Preserves exact source text and returns every molecule with coordinates/IDs.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=None, max_molecules=None, max_atoms=None, max_bonds=None))]
pub(super) fn parse_cdjson<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_molecules: Option<usize>,
    max_atoms: Option<usize>,
    max_bonds: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let defaults = chematic_mol::CdjsonParseLimits::default();
    let limits = chematic_mol::CdjsonParseLimits {
        max_input_bytes: max_input_bytes.unwrap_or(defaults.max_input_bytes),
        max_molecules: max_molecules.unwrap_or(defaults.max_molecules),
        max_atoms: max_atoms.unwrap_or(defaults.max_atoms),
        max_bonds: max_bonds.unwrap_or(defaults.max_bonds),
    };
    let document = chematic_mol::parse_cdjson_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mut molecules = Vec::new();
    for item in document.molecules() {
        let value = PyDict::new(py);
        value.set_item("mol", Mol::bare(item.molecule.clone()))?;
        value.set_item("coordinates", &item.coordinates)?;
        value.set_item("id", &item.id)?;
        value.set_item("atom_ids", &item.atom_ids)?;
        value.set_item("bond_ids", &item.bond_ids)?;
        molecules.push(value);
    }
    let out = PyDict::new(py);
    out.set_item("source_text", document.source_text())?;
    out.set_item("molecules", molecules)?;
    Ok(out)
}

/// Export one non-aromatic, non-stereochemical molecular graph in drawing units.
#[pyfunction]
pub(super) fn write_cdjson_graph(mol: &Mol, coordinates: &Bound<'_, PyList>) -> PyResult<String> {
    let hard = chematic_mol::CdjsonParseLimits::default();
    if mol.inner.atom_count() > hard.max_atoms || mol.inner.bond_count() > hard.max_bonds {
        return Err(PyValueError::new_err(
            "CDJSON graph exceeds atom/bond limit",
        ));
    }
    if coordinates.len() != mol.inner.atom_count() {
        return Err(PyValueError::new_err(
            "CDJSON coordinate count must equal atom count",
        ));
    }
    let coordinates = coordinates.extract::<Vec<[f64; 3]>>()?;
    chematic_mol::write_cdjson_graph(&mol.inner, &coordinates)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
