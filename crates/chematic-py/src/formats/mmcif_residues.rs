use super::*;

fn residue_identity<'py>(
    py: Python<'py>,
    id: &chematic_mol::MmcifResidueIdentity,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("model_num", id.model_num)?;
    d.set_item("label_asym_id", &id.label_asym_id)?;
    d.set_item("auth_asym_id", &id.auth_asym_id)?;
    d.set_item("entity_id", &id.entity_id)?;
    d.set_item("label_seq_id", id.label_seq_id)?;
    d.set_item("auth_seq_id", &id.auth_seq_id)?;
    d.set_item("insertion_code", &id.insertion_code)?;
    d.set_item("label_comp_id", &id.label_comp_id)?;
    d.set_item("auth_comp_id", &id.auth_comp_id)?;
    Ok(d)
}

/// General residue view including glycans, ligands and water. No bonds or
/// glycan connectivity are inferred. Identity fields are not truncated.
#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_mmcif_residues<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::MmcifParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let doc = chematic_mol::parse_mmcif_residues_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    mmcif_residues_to_pydict(py, &doc)
}

pub(super) fn mmcif_residues_to_pydict<'py>(
    py: Python<'py>,
    doc: &chematic_mol::MmcifResidueDocument,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let atoms = doc
        .parsed
        .atoms
        .iter()
        .zip(&doc.atom_ids)
        .enumerate()
        .map(|(i, (a, id))| {
            let row = materials::mmcif_atom_to_pydict(py, a)?;
            row.set_item("label_atom_id", &id.label_atom_id)?;
            row.set_item("auth_atom_id", &id.auth_atom_id)?;
            row.set_item("alt_loc", &id.alt_id)?;
            row.set_item(
                "icode",
                &doc.residues[doc.atom_to_residue[i]].identity.insertion_code,
            )?;
            Ok(row)
        })
        .collect::<PyResult<Vec<_>>>()?;
    let residues = doc
        .residues
        .iter()
        .map(|r| {
            let group = residue_identity(py, &r.identity)?;
            group.set_item("atom_indices", &r.atom_indices)?;
            Ok(group)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("atoms", atoms)?;
    d.set_item("residues", residues)?;
    d.set_item("atom_to_residue", &doc.atom_to_residue)?;
    match &doc.parsed.cell {
        Some(cell) => d.set_item("cell", materials::unit_cell_to_pydict(py, cell)?)?,
        None => d.set_item("cell", py.None())?,
    }
    d.set_item("space_group", &doc.parsed.space_group)?;
    d.set_item("unhandled_columns", &doc.parsed.unhandled_columns)?;
    Ok(d)
}
