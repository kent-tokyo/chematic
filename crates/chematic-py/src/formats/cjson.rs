use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Parse a ChemicalJSON (.cjson) string.
///
/// Returns ``(mol, coords)`` where ``coords`` is a list of ``[x, y, z]``
/// coordinate triples (Å), one per heavy atom.  ``coords`` is empty when
/// the file has no ``atoms.coords.3d`` field.
///
/// ChemicalJSON is the native format of Avogadro 2 and the MolSSI
/// Open Chemistry toolkit.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, coords = chematic.from_cjson(open("mol.cjson").read())
///     print(mol.smiles)
///     # Round-trip:
///     open("out.cjson", "w").write(mol.to_cjson(coords))
#[pyfunction]
pub(super) fn from_cjson(cjson_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, coords) =
        chematic_mol::parse_cjson(cjson_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let py_coords = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
    Ok((
        Mol {
            num_rings_cache: Default::default(),
            inner: Arc::new(mol),
            props: Default::default(),
        },
        py_coords,
    ))
}
