//! Explicit bounded geometry arrays; Python iterators cannot bypass the atom cap.
use chematic_mol::{GjfCoordinateUnits, GjfWriteOptions};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyInt, PyString};

fn strict_charge(value: &Bound<'_, PyAny>) -> PyResult<i32> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "charge/freeze codes must be integers, excluding bool",
        ));
    }
    value
        .extract::<i32>()
        .map_err(|_| PyValueError::new_err("charge/freeze code must fit i32"))
}
fn strict_multiplicity(value: &Bound<'_, PyAny>) -> PyResult<u32> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "multiplicity must be an integer, excluding bool",
        ));
    }
    value
        .extract::<u32>()
        .map_err(|_| PyValueError::new_err("multiplicity must fit u32"))
}
fn elements(value: &Bound<'_, PyAny>) -> PyResult<Vec<chematic_core::Element>> {
    let mut output = Vec::new();
    for item in value.try_iter()? {
        if output.len() == GjfWriteOptions::MAX_ATOMS {
            return Err(PyValueError::new_err("atom count must be 1..=100000"));
        }
        let item = item?;
        let symbol = item.cast::<PyString>()?.to_str()?;
        let element = (symbol.len() <= 2)
            .then(|| chematic_core::Element::from_symbol(symbol))
            .flatten()
            .filter(|e| e.atomic_number() > 0 && e.symbol() == symbol)
            .ok_or_else(|| PyValueError::new_err("expected a canonical elemental symbol"))?;
        output.push(element);
    }
    Ok(output)
}
fn positions(value: &Bound<'_, PyAny>, count: usize) -> PyResult<Vec<[f64; 3]>> {
    let mut output = Vec::new();
    for row in value.try_iter()? {
        if output.len() == count {
            return Err(PyValueError::new_err("too many coordinate rows"));
        }
        let row = row?;
        let mut fields = row.try_iter()?;
        let mut point = [0.0; 3];
        for component in &mut point {
            let item = fields.next().transpose()?.ok_or_else(|| {
                PyValueError::new_err("coordinate rows must contain exactly three values")
            })?;
            if item.is_instance_of::<PyBool>() {
                return Err(PyValueError::new_err(
                    "coordinates must be numbers, excluding bool",
                ));
            }
            *component = item.extract()?;
        }
        if fields.next().transpose()?.is_some() {
            return Err(PyValueError::new_err(
                "coordinate rows must contain exactly three values",
            ));
        }
        output.push(point);
    }
    Ok(output)
}
fn codes(value: Option<&Bound<'_, PyAny>>, count: usize) -> PyResult<Vec<Option<i32>>> {
    let Some(value) = value else {
        return Ok(vec![None; count]);
    };
    let mut output = Vec::new();
    for item in value.try_iter()? {
        if output.len() == count {
            return Err(PyValueError::new_err("too many freeze codes"));
        }
        let item = item?;
        output.push(if item.is_none() {
            None
        } else {
            Some(strict_charge(&item)?)
        });
    }
    Ok(output)
}

/// Write explicit Cartesian geometry without rounding loss or zero filling.
///
/// Required units are "angstrom" or "bohr"; numbers are unchanged. Absent
/// freeze_codes means no flags, while explicit None/0/-1 remain distinct.
/// Uses 17 significant digits with bitwise numeric readback, including signed
/// zero/subnormals. No topology/isotope/job-source preservation is implied.
/// Header checks are lexical; calculation keywords and execution are not verified.
#[pyfunction]
#[pyo3(signature=(elements, coords, *, charge, multiplicity, units, freeze_codes=None, title="chematic Cartesian geometry", route="#"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_gjf_cartesian(
    elements: &Bound<'_, PyAny>,
    coords: &Bound<'_, PyAny>,
    #[pyo3(from_py_with = strict_charge)] charge: i32,
    #[pyo3(from_py_with = strict_multiplicity)] multiplicity: u32,
    units: &str,
    freeze_codes: Option<&Bound<'_, PyAny>>,
    title: &str,
    route: &str,
) -> PyResult<String> {
    if title.len() > 1024 || route.len() > 1024 {
        return Err(PyValueError::new_err("headers must be at most 1024 bytes"));
    }
    let elements = self::elements(elements)?;
    let positions = positions(coords, elements.len())?;
    let codes = codes(freeze_codes, elements.len())?;
    let units = match units {
        "angstrom" => GjfCoordinateUnits::Angstrom,
        "bohr" => GjfCoordinateUnits::Bohr,
        _ => return Err(PyValueError::new_err("units must be 'angstrom' or 'bohr'")),
    };
    let mut options = GjfWriteOptions::new(charge, multiplicity, units);
    options.title = title.into();
    options.route = route.into();
    chematic_mol::write_gjf_cartesian(&elements, &positions, &codes, &options)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
