use super::*;

fn dmol3_error(e: chematic_mol::Dmol3Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_dmol3<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::Dmol3ParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let g = chematic_mol::parse_dmol3_with_limits(text, &limits).map_err(dmol3_error)?;
    let out = PyDict::new(py);
    out.set_item(
        "elements",
        g.atoms
            .iter()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "positions_bohr",
        g.atoms.iter().map(|a| a.position_bohr).collect::<Vec<_>>(),
    )?;
    out.set_item(
        "positions_angstrom",
        g.positions_angstrom().map_err(dmol3_error)?,
    )?;
    out.set_item("cell_bohr", g.cell_bohr)?;
    out.set_item("cell_angstrom", g.cell_angstrom().map_err(dmol3_error)?)?;
    out.set_item("bohr_angstrom", chematic_mol::DMOL3_BOHR_ANGSTROM)?;
    out.set_item("periodicity", g.periodicity())?;
    out.set_item("source_text", g.source_text)?;
    Ok(out)
}

#[pyfunction]
#[pyo3(signature=(elements, positions_bohr, *, cell_bohr=None))]
pub(super) fn write_dmol3(
    elements: Vec<String>,
    positions_bohr: Vec<[f64; 3]>,
    cell_bohr: Option<[[f64; 3]; 3]>,
) -> PyResult<String> {
    if elements.len() != positions_bohr.len() {
        return Err(PyValueError::new_err(
            "element and coordinate counts must match",
        ));
    }
    let atoms = elements
        .into_iter()
        .zip(positions_bohr)
        .map(|(symbol, position_bohr)| {
            let element = chematic_core::Element::from_symbol(&symbol)
                .filter(|e| e.symbol() == symbol)
                .ok_or_else(|| {
                    PyValueError::new_err(format!("invalid canonical element '{symbol}'"))
                })?;
            Ok(chematic_mol::Dmol3Atom {
                element,
                position_bohr,
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    chematic_mol::write_dmol3(&chematic_mol::Dmol3Geometry {
        atoms,
        cell_bohr,
        source_text: String::new(),
    })
    .map_err(dmol3_error)
}
