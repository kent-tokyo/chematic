use super::*;

type PdbqtData = (Mol, Vec<[f64; 3]>, Vec<f64>);

/// Parse PDBQT atom elements only. Coordinates and partial charges are discarded.
/// No bonds are inferred; use :func:`from_pdbqt_with_data` to retain numeric data.
#[pyfunction]
pub(super) fn from_pdbqt(pdbqt_str: &str) -> PyResult<Mol> {
    from_pdbqt_with_data(pdbqt_str).map(|(mol, _, _)| mol)
}

/// Parse PDBQT into ``(mol, coords, partial_charges)`` in source atom order.
/// Coordinates are [x,y,z] arrays in Å and partial charges are in e, retaining
/// the parsed finite binary64 values without output rounding or recalculation.
/// Partial charges do not set formal charges. No bonds are inferred. Source
/// atom identities, AutoDock types, remarks and torsion trees are not retained.
/// Accepts str or UTF-8 bytes; malformed atom records raise ValueError.
#[pyfunction]
pub(super) fn from_pdbqt_with_data(
    #[pyo3(from_py_with = serialized_text)] pdbqt_str: &str,
) -> PyResult<PdbqtData> {
    let (mol, coords, charges) =
        chematic_mol::parse_pdbqt(pdbqt_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok((
        Mol {
            num_rings_cache: Default::default(),
            inner: Arc::new(mol),
            props: Default::default(),
        },
        coords.into_iter().map(|(x, y, z)| [x, y, z]).collect(),
        charges,
    ))
}

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
