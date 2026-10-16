use super::*;

/// Write rigid PDBQT with exact coordinate/charge counts and checked fixed widths.
/// Each graph atom, including explicit H, needs coordinates in Å and a partial
/// charge in e. Raises ValueError on missing/excess values, non-finite numbers,
/// overflow after rounding, invalid names or unsupported default AutoDock types.
/// Types retain the legacy heuristic; this is not ligand preparation.
#[pyfunction]
#[pyo3(signature = (mol, coords, charges, name="LIG"))]
pub(super) fn write_pdbqt_checked(
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    charges: Vec<f64>,
    name: &str,
) -> PyResult<String> {
    let coords: Vec<_> = coords.into_iter().map(|p| (p[0], p[1], p[2])).collect();
    chematic_mol::write_pdbqt_checked(&mol.inner, &coords, &charges, name)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
