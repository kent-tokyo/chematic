use super::*;
use chematic_mol::{FhiaimsAtom, FhiaimsCoordinates, FhiaimsGeometry, FhiaimsLatticeVector};

#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_fhiaims<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::FhiaimsParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let g = chematic_mol::parse_fhiaims_with_limits(text, &limits).map_err(fhiaims_error)?;
    let out = PyDict::new(py);
    out.set_item("periodicity", g.periodicity())?;
    out.set_item(
        "positions_angstrom",
        g.positions_angstrom().map_err(fhiaims_error)?,
    )?;
    out.set_item(
        "elements",
        g.atoms
            .iter()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "coordinate_kinds",
        g.atoms
            .iter()
            .map(|a| match a.coordinates {
                FhiaimsCoordinates::Cartesian(_) => "cartesian",
                FhiaimsCoordinates::Fractional(_) => "fractional",
            })
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "coordinates",
        g.atoms
            .iter()
            .map(|a| a.coordinates.values())
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "constraints",
        g.atoms.iter().map(|a| a.constraints).collect::<Vec<_>>(),
    )?;
    out.set_item(
        "initial_charges",
        g.atoms.iter().map(|a| a.initial_charge).collect::<Vec<_>>(),
    )?;
    out.set_item(
        "initial_moments",
        g.atoms.iter().map(|a| a.initial_moment).collect::<Vec<_>>(),
    )?;
    out.set_item(
        "lattice_angstrom",
        g.lattice_vectors
            .iter()
            .map(|v| v.vector_angstrom)
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "lattice_constraints",
        g.lattice_vectors
            .iter()
            .map(|v| v.constraints)
            .collect::<Vec<_>>(),
    )?;
    out.set_item("source_text", g.source_text)?;
    Ok(out)
}

fn fhiaims_error(e: chematic_mol::FhiaimsError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

#[pyfunction]
#[pyo3(signature=(elements, coordinates, *, coordinate_kinds=None, lattice_angstrom=None, constraints=None, initial_charges=None, initial_moments=None, lattice_constraints=None))]
#[allow(clippy::too_many_arguments)]
pub(super) fn write_fhiaims(
    elements: Vec<String>,
    coordinates: Vec<[f64; 3]>,
    coordinate_kinds: Option<Vec<String>>,
    lattice_angstrom: Option<Vec<[f64; 3]>>,
    constraints: Option<Vec<[bool; 3]>>,
    initial_charges: Option<Vec<Option<f64>>>,
    initial_moments: Option<Vec<Option<f64>>>,
    lattice_constraints: Option<Vec<[bool; 3]>>,
) -> PyResult<String> {
    let n = elements.len();
    let lengths = [
        Some(coordinates.len()),
        coordinate_kinds.as_ref().map(Vec::len),
        constraints.as_ref().map(Vec::len),
        initial_charges.as_ref().map(Vec::len),
        initial_moments.as_ref().map(Vec::len),
    ];
    if lengths.into_iter().flatten().any(|len| len != n) {
        return Err(PyValueError::new_err("all atom-array lengths must match"));
    }
    let lattice = lattice_angstrom.unwrap_or_default();
    if lattice_constraints
        .as_ref()
        .is_some_and(|v| v.len() != lattice.len())
    {
        return Err(PyValueError::new_err(
            "lattice constraint count must match vector count",
        ));
    }
    let mut atoms = Vec::new();
    for (i, (symbol, values)) in elements.into_iter().zip(coordinates).enumerate() {
        let element = chematic_core::Element::from_symbol(&symbol)
            .filter(|e| e.symbol() == symbol)
            .ok_or_else(|| {
                PyValueError::new_err(format!("invalid canonical element '{symbol}'"))
            })?;
        let kind = coordinate_kinds
            .as_ref()
            .map_or("cartesian", |v| v[i].as_str());
        let coordinates = match kind {
            "cartesian" => FhiaimsCoordinates::Cartesian(values),
            "fractional" => FhiaimsCoordinates::Fractional(values),
            _ => {
                return Err(PyValueError::new_err(
                    "coordinate kind must be cartesian or fractional",
                ));
            }
        };
        atoms.push(FhiaimsAtom {
            element,
            coordinates,
            constraints: constraints.as_ref().map_or([false; 3], |v| v[i]),
            initial_charge: initial_charges.as_ref().and_then(|v| v[i]),
            initial_moment: initial_moments.as_ref().and_then(|v| v[i]),
        });
    }
    let lattice_vectors = lattice
        .into_iter()
        .enumerate()
        .map(|(i, vector_angstrom)| FhiaimsLatticeVector {
            vector_angstrom,
            constraints: lattice_constraints.as_ref().map_or([false; 3], |v| v[i]),
        })
        .collect();
    chematic_mol::write_fhiaims(&FhiaimsGeometry {
        atoms,
        lattice_vectors,
        source_text: String::new(),
    })
    .map_err(fhiaims_error)
}
