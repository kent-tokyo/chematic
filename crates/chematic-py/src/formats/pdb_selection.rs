use super::*;

#[pyfunction]
#[pyo3(signature=(text, atom_indices=None, *, allow_metadata_loss=false, max_input_bytes=None, max_atoms=None, max_line_len=None, max_records=None))]
#[allow(clippy::too_many_arguments)]
pub(super) fn select_pdb_atoms<'py>(
    py: Python<'py>,
    text: &str,
    atom_indices: Option<Vec<usize>>,
    allow_metadata_loss: bool,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
    max_records: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::PdbSelectionLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    if let Some(v) = max_records {
        limits.max_records = v;
    }
    let r = chematic_mol::select_pdb_atoms_with_limits(
        text,
        atom_indices.as_deref(),
        allow_metadata_loss,
        &limits,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let d = PyDict::new(py);
    d.set_item("pdb", r.pdb)?;
    d.set_item("source_text", r.source_text)?;
    d.set_item("source_atom_count", r.source_atom_count)?;
    d.set_item("selected_source_indices", r.selected_source_indices)?;
    d.set_item("old_to_new", r.old_to_new)?;
    d.set_item("model_num", r.model_num)?;
    d.set_item("retargeted_ter_lines", r.retargeted_ter_lines)?;
    d.set_item("filtered_conect_lines", r.filtered_conect_lines)?;
    d.set_item(
        "atoms",
        r.atoms
            .iter()
            .map(|a| atom(py, a))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    let losses = r
        .losses
        .iter()
        .map(|loss| {
            let d = PyDict::new(py);
            d.set_item("line", loss.line)?;
            d.set_item("record_name", &loss.record_name)?;
            Ok(d)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("losses", losses)?;
    Ok(d)
}
fn atom<'py>(py: Python<'py>, a: &chematic_mol::PdbSelectedAtom) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("source_index", a.source_index)?;
    d.set_item("serial", a.serial)?;
    d.set_item("record_type", &a.record_type)?;
    d.set_item("atom_name", &a.atom_name)?;
    d.set_item("alt_loc", a.alt_loc.to_string())?;
    d.set_item("residue_name", &a.residue_name)?;
    d.set_item("chain_id", a.chain_id.to_string())?;
    d.set_item("residue_number", a.residue_number)?;
    d.set_item("insertion_code", a.insertion_code.to_string())?;
    d.set_item("position_angstrom", a.position_angstrom)?;
    Ok(d)
}
