use super::*;

/// Read one ORCA-AICCM EXYZ document, retaining lattice, offset and virtual flags.
#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_exyz<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::ExyzParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let g = chematic_mol::parse_exyz_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let out = PyDict::new(py);
    out.set_item("header", g.header)?;
    out.set_item(
        "elements",
        g.atoms
            .iter()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "positions_angstrom",
        g.atoms
            .iter()
            .map(|a| a.position_angstrom)
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "virtual_atoms",
        g.atoms.iter().map(|a| a.virtual_atom).collect::<Vec<_>>(),
    )?;
    out.set_item(
        "lattice_angstrom",
        g.cell.as_ref().map(|c| c.lattice_angstrom),
    )?;
    out.set_item(
        "offset_angstrom",
        g.cell.as_ref().map(|c| c.offset_angstrom),
    )?;
    out.set_item("source_text", g.source_text)?;
    Ok(out)
}

/// Write explicit EXYZ geometry. Adds required header directives when absent.
#[pyfunction]
#[pyo3(signature=(elements, positions_angstrom, *, header="", lattice_angstrom=None, offset_angstrom=None, virtual_atoms=None))]
pub(super) fn write_exyz(
    elements: Vec<String>,
    positions_angstrom: Vec<[f64; 3]>,
    header: &str,
    lattice_angstrom: Option<[[f64; 3]; 3]>,
    offset_angstrom: Option<[f64; 3]>,
    virtual_atoms: Option<Vec<bool>>,
) -> PyResult<String> {
    if elements.len() != positions_angstrom.len()
        || virtual_atoms
            .as_ref()
            .is_some_and(|v| v.len() != elements.len())
    {
        return Err(PyValueError::new_err(
            "element, position and virtual-flag counts must match",
        ));
    }
    if offset_angstrom.is_some() && lattice_angstrom.is_none() {
        return Err(PyValueError::new_err(
            "cell offset requires lattice vectors",
        ));
    }
    let mut header = header.to_string();
    for (enabled, flag) in [
        (lattice_angstrom.is_some(), "%PBC"),
        (
            virtual_atoms.as_ref().is_some_and(|v| v.iter().any(|v| *v)),
            "%VIRTUAL",
        ),
    ] {
        if enabled && !header.split_whitespace().any(|s| s == flag) {
            if !header.is_empty() {
                header.push(' ');
            }
            header.push_str(flag);
        }
    }
    let mut atoms = Vec::new();
    for (i, (symbol, position_angstrom)) in elements.into_iter().zip(positions_angstrom).enumerate()
    {
        let element = chematic_core::Element::from_symbol(&symbol)
            .filter(|e| e.symbol() == symbol)
            .ok_or_else(|| {
                PyValueError::new_err(format!("invalid canonical element symbol '{symbol}'"))
            })?;
        atoms.push(chematic_mol::ExyzAtom {
            element,
            position_angstrom,
            virtual_atom: virtual_atoms.as_ref().is_some_and(|v| v[i]),
        });
    }
    let g = chematic_mol::ExyzGeometry {
        header,
        atoms,
        cell: lattice_angstrom.map(|lattice_angstrom| chematic_mol::ExyzCell {
            lattice_angstrom,
            offset_angstrom: offset_angstrom.unwrap_or([0.0; 3]),
        }),
        source_text: String::new(),
    };
    chematic_mol::write_exyz(&g).map_err(|e| PyValueError::new_err(e.to_string()))
}
