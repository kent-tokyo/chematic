use super::*;
use pyo3::types::{PyBool, PyInt};

fn index(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "reference indices must be integers, excluding bool",
        ));
    }
    value
        .extract::<usize>()
        .map_err(|_| PyValueError::new_err("reference index must be nonnegative and fit usize"))
}
fn strict_u8(value: &Bound<'_, PyAny>) -> PyResult<u8> {
    index(value)?
        .try_into()
        .map_err(|_| PyValueError::new_err("integer must fit u8"))
}
fn strict_i32(value: &Bound<'_, PyAny>) -> PyResult<i32> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "charge must be an integer, excluding bool",
        ));
    }
    value
        .extract::<i32>()
        .map_err(|_| PyValueError::new_err("charge must fit i32"))
}
fn strict_u32(value: &Bound<'_, PyAny>) -> PyResult<u32> {
    index(value)?
        .try_into()
        .map_err(|_| PyValueError::new_err("multiplicity must fit u32"))
}
fn reference_rows(references: &Bound<'_, PyAny>, count: usize) -> PyResult<Vec<Vec<usize>>> {
    let mut refs = Vec::new();
    for row in references.try_iter()? {
        if refs.len() >= count {
            return Err(PyValueError::new_err("too many reference rows"));
        }
        let mut indices = Vec::new();
        for v in row?.try_iter()? {
            if indices.len() >= refs.len().min(3) {
                return Err(PyValueError::new_err("too many references for atom row"));
            }
            indices.push(index(&v?)?);
        }
        refs.push(indices);
    }
    Ok(refs)
}
fn convert(
    elements: Vec<String>,
    positions: Vec<[f64; 3]>,
    references: &Bound<'_, PyAny>,
) -> PyResult<Vec<chematic_mol::ZMatrixRow>> {
    if elements.is_empty() || elements.len() > chematic_mol::gzmat::MAX_ZMATRIX_ATOMS {
        return Err(PyValueError::new_err("atom count must be 1..=100000"));
    }
    let refs = reference_rows(references, elements.len())?;
    let elements = elements
        .iter()
        .map(|s| {
            chematic_core::Element::from_symbol(s)
                .filter(|e| e.symbol() == s)
                .ok_or_else(|| PyValueError::new_err(format!("invalid canonical element '{s}'")))
        })
        .collect::<PyResult<Vec<_>>>()?;
    chematic_mol::cartesian_to_zmatrix(&elements, &positions, &refs)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Convert Cartesian angstroms to caller-referenced internal rows.
/// Each dict contains element, zero-based references, and values (angstrom, degrees).
#[pyfunction]
pub(super) fn cartesian_to_zmatrix<'py>(
    py: Python<'py>,
    elements: Vec<String>,
    positions_angstrom: Vec<[f64; 3]>,
    references: &Bound<'py, PyAny>,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    convert(elements, positions_angstrom, references)?
        .iter()
        .map(|row| {
            let result = PyDict::new(py);
            result.set_item("element", row.element.symbol())?;
            result.set_item("references", &row.references)?;
            result.set_item("values", &row.values)?;
            Ok(result)
        })
        .collect()
}

/// Write caller-referenced internal coordinates with 1..17 significant digits.
/// Explicit charge/spin; no automatic topology, H, isotope or dummy-atom inference.
#[pyfunction]
#[pyo3(signature=(elements,positions_angstrom,references,*,charge,multiplicity,significant_digits=17,title="chematic internal coordinates",route="#"))]
#[allow(clippy::too_many_arguments)]
pub(super) fn write_gzmat(
    elements: Vec<String>,
    positions_angstrom: Vec<[f64; 3]>,
    references: &Bound<'_, PyAny>,
    #[pyo3(from_py_with=strict_i32)] charge: i32,
    #[pyo3(from_py_with=strict_u32)] multiplicity: u32,
    #[pyo3(from_py_with=strict_u8)] significant_digits: u8,
    title: &str,
    route: &str,
) -> PyResult<String> {
    let rows = convert(elements, positions_angstrom, references)?;
    let mut options = chematic_mol::GzmatWriteOptions::new(charge, multiplicity);
    options.significant_digits = significant_digits;
    options.title = title.into();
    options.route = route.into();
    chematic_mol::write_gzmat(&rows, &options).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Write GAMESS-US Cartesian $DATA plus explicit caller-defined $ZMAT internals.
/// A geometry template only: charge/spin are explicit; no basis or method is chosen.
#[pyfunction]
#[pyo3(signature=(elements,positions_angstrom,references,*,charge,multiplicity,significant_digits=17,title="chematic GAMESS internal coordinates"))]
pub(super) fn write_gamess_zmat(
    elements: Vec<String>,
    positions_angstrom: Vec<[f64; 3]>,
    references: &Bound<'_, PyAny>,
    #[pyo3(from_py_with=strict_i32)] charge: i32,
    #[pyo3(from_py_with=strict_u32)] multiplicity: u32,
    #[pyo3(from_py_with=strict_u8)] significant_digits: u8,
    title: &str,
) -> PyResult<String> {
    if elements.is_empty()
        || elements.len() > chematic_mol::gamess_zmat::MAX_GAMESS_ZMAT_ATOMS
        || elements.len() != positions_angstrom.len()
    {
        return Err(PyValueError::new_err(
            "expected matching element/position counts in 1..=1024",
        ));
    }
    let rows = convert(elements, positions_angstrom.clone(), references)?;
    let elements: Vec<_> = rows.iter().map(|row| row.element).collect();
    let references: Vec<_> = rows.iter().map(|row| row.references.clone()).collect();
    if title.len() > 80 {
        return Err(PyValueError::new_err("GAMESS title must fit 80 bytes"));
    }
    let mut options = chematic_mol::GamessZmatWriteOptions::new(charge, multiplicity);
    options.significant_digits = significant_digits;
    options.title = title.into();
    chematic_mol::write_gamess_zmat(&elements, &positions_angstrom, &references, &options)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_pyfunctions!(
        m,
        cartesian_to_zmatrix,
        write_gzmat,
        write_gamess_zmat,
        xyz_to_gamess_zmat
    )
}

/// Convert one plain XYZ frame to GAMESS-US Cartesian $DATA and explicit $ZMAT.
#[pyfunction]
#[pyo3(signature=(text,references,*,charge,multiplicity,significant_digits=17,title="chematic GAMESS internal coordinates"))]
pub(super) fn xyz_to_gamess_zmat(
    text: &str,
    references: &Bound<'_, PyAny>,
    #[pyo3(from_py_with=strict_i32)] charge: i32,
    #[pyo3(from_py_with=strict_u32)] multiplicity: u32,
    #[pyo3(from_py_with=strict_u8)] significant_digits: u8,
    title: &str,
) -> PyResult<String> {
    let refs = reference_rows(references, chematic_mol::gamess_zmat::MAX_GAMESS_ZMAT_ATOMS)?;
    if title.len() > 80 {
        return Err(PyValueError::new_err("GAMESS title must fit 80 bytes"));
    }
    let mut options = chematic_mol::GamessZmatWriteOptions::new(charge, multiplicity);
    options.significant_digits = significant_digits;
    options.title = title.into();
    chematic_mol::xyz_to_gamess_zmat(text, &refs, &options)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
