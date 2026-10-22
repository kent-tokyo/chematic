use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Parse the immutable-record CJSON profile into ``(Mol, coords)``.
/// Rejects malformed known arrays, duplicate decoded keys, invalid/duplicate
/// bonds and unsupported order tokens rather than filling or rounding them.
/// Unknown field semantics, stereo and arbitrary metadata are not interpreted.
/// Limits: 16 MiB input, 100000 atoms plus native CJSON tree/bond limits.
#[pyfunction]
pub(super) fn from_cjson_checked(cjson_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let limits = chematic_mol::CjsonParseLimits {
        max_input_bytes: 16 * 1024 * 1024,
        max_atoms: 100_000,
        ..Default::default()
    };
    let (mol, coords) = chematic_mol::parse_cjson_checked_with_limits(cjson_str, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok((
        Mol {
            num_rings_cache: Default::default(),
            inner: Arc::new(mol),
            props: Default::default(),
        },
        coords.into_iter().map(|(x, y, z)| vec![x, y, z]).collect(),
    ))
}

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
