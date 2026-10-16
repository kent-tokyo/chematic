use super::*;

/// Project ATOM/HETATM and explicit CONECT into a graph, coordinates and serials.
/// Never adds proximity bonds. Extended target slots require explicit opt-in.
#[pyfunction]
#[pyo3(signature = (pdb_str, *, allow_extended_targets=false))]
pub(super) fn from_pdb_conect<'py>(
    py: Python<'py>,
    pdb_str: &str,
    allow_extended_targets: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let frame = chematic_mol::parse_pdb_conect_with_options(
        pdb_str,
        &chematic_mol::PdbConectOptions {
            allow_extended_targets,
        },
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let out = PyDict::new(py);
    out.set_item(
        "mol",
        Mol {
            inner: Arc::new(frame.molecule),
            props: Default::default(),
        },
    )?;
    out.set_item(
        "coords",
        frame
            .coords
            .points
            .iter()
            .map(|p| vec![p.x, p.y, p.z])
            .collect::<Vec<_>>(),
    )?;
    out.set_item("serials", frame.serials)?;
    Ok(out)
}
