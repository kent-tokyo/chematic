use super::*;
use chematic_mol::{MopacAtom, MopacCartesianGeometry, MopacCartesianRow, MopacOptimizationFlag};

fn mopac_error(e: chematic_mol::MopacError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn flags(rows: &[MopacCartesianRow]) -> Vec<Option<Vec<&'static str>>> {
    rows.iter()
        .map(|row| {
            row.optimization_flags
                .map(|f| f.iter().map(|v| v.code()).collect())
        })
        .collect()
}

#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_mopac_cartesian<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::MopacParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let g = chematic_mol::parse_mopac_cartesian_with_limits(text, &limits).map_err(mopac_error)?;
    let out = PyDict::new(py);
    out.set_item("keywords", &g.keywords)?;
    out.set_item("title", &g.title)?;
    out.set_item("comment", &g.comment)?;
    out.set_item(
        "elements",
        g.atoms
            .iter()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "labels",
        g.atoms.iter().map(|a| &a.label).collect::<Vec<_>>(),
    )?;
    out.set_item(
        "positions_angstrom",
        g.atoms
            .iter()
            .map(|a| a.row.position_angstrom)
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "optimization_flags",
        g.atoms
            .iter()
            .map(|a| {
                a.row
                    .optimization_flags
                    .map(|f| f.iter().map(|v| v.code()).collect::<Vec<_>>())
            })
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "translation_vectors_angstrom",
        g.translation_vectors
            .iter()
            .map(|r| r.position_angstrom)
            .collect::<Vec<_>>(),
    )?;
    out.set_item("translation_flags", flags(&g.translation_vectors))?;
    out.set_item("periodicity", g.periodicity())?;
    out.set_item("source_text", &g.source_text)?;
    Ok(out)
}

fn row(position_angstrom: [f64; 3], codes: Option<[String; 3]>) -> PyResult<MopacCartesianRow> {
    let optimization_flags = codes
        .map(|codes| {
            let convert = |s: &str| {
                MopacOptimizationFlag::from_code(s)
                    .ok_or_else(|| PyValueError::new_err("unsupported MOPAC optimization flag"))
            };
            Ok::<_, PyErr>([
                convert(&codes[0])?,
                convert(&codes[1])?,
                convert(&codes[2])?,
            ])
        })
        .transpose()?;
    Ok(MopacCartesianRow {
        position_angstrom,
        optimization_flags,
    })
}

#[pyfunction]
#[pyo3(signature=(elements, positions_angstrom, *, keywords="", title="", comment="", labels=None, optimization_flags=None, translation_vectors_angstrom=None, translation_flags=None))]
#[allow(clippy::too_many_arguments)]
pub(super) fn write_mopac_cartesian(
    elements: Vec<String>,
    positions_angstrom: Vec<[f64; 3]>,
    keywords: &str,
    title: &str,
    comment: &str,
    labels: Option<Vec<Option<String>>>,
    optimization_flags: Option<Vec<Option<[String; 3]>>>,
    translation_vectors_angstrom: Option<Vec<[f64; 3]>>,
    translation_flags: Option<Vec<Option<[String; 3]>>>,
) -> PyResult<String> {
    let n = elements.len();
    if [
        Some(positions_angstrom.len()),
        labels.as_ref().map(Vec::len),
        optimization_flags.as_ref().map(Vec::len),
    ]
    .into_iter()
    .flatten()
    .any(|len| len != n)
    {
        return Err(PyValueError::new_err("all atom-array lengths must match"));
    }
    let vectors = translation_vectors_angstrom.unwrap_or_default();
    if translation_flags
        .as_ref()
        .is_some_and(|v| v.len() != vectors.len())
    {
        return Err(PyValueError::new_err(
            "translation flag count must match vector count",
        ));
    }
    let atoms = elements
        .into_iter()
        .zip(positions_angstrom)
        .enumerate()
        .map(|(i, (symbol, position))| {
            let element = chematic_core::Element::from_symbol(&symbol)
                .filter(|e| e.symbol() == symbol)
                .ok_or_else(|| {
                    PyValueError::new_err(format!("invalid canonical element '{symbol}'"))
                })?;
            Ok(MopacAtom {
                element,
                label: labels.as_ref().and_then(|v| v[i].clone()),
                row: row(
                    position,
                    optimization_flags.as_ref().and_then(|v| v[i].clone()),
                )?,
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let translation_vectors = vectors
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            row(
                position,
                translation_flags.as_ref().and_then(|v| v[i].clone()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    chematic_mol::write_mopac_cartesian(&MopacCartesianGeometry {
        keywords: keywords.to_string(),
        title: title.to_string(),
        comment: comment.to_string(),
        atoms,
        translation_vectors,
        source_text: String::new(),
    })
    .map_err(mopac_error)
}
