//! Plain CIF atomic-position bindings.
use super::*;

/// Parse a CIF (Crystallographic Information File) string and return a dict.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "cell": dict | None}``
///     where ``cell`` has keys ``a, b, c, alpha, beta, gamma``.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     result = chematic.parse_cif(open("structure.cif").read())
///     mol = result["mol"]
///     coords = result["coords"]
#[pyfunction]
pub(super) fn parse_cif<'py>(
    py: Python<'py>,
    cif_str: &str,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let result =
        chematic_mol::parse_cif(cif_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mol = Mol {
        inner: Arc::new(result.mol),
        props: Default::default(),
    };
    let coords: Vec<Vec<f64>> = result
        .coords
        .iter()
        .map(|&(x, y, z)| vec![x, y, z])
        .collect();
    let d = pyo3::types::PyDict::new(py);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    if let Some(cell) = result.cell {
        let cd = pyo3::types::PyDict::new(py);
        cd.set_item("a", cell.a)?;
        cd.set_item("b", cell.b)?;
        cd.set_item("c", cell.c)?;
        cd.set_item("alpha", cell.alpha)?;
        cd.set_item("beta", cell.beta)?;
        cd.set_item("gamma", cell.gamma)?;
        d.set_item("cell", cd)?;
    } else {
        d.set_item("cell", py.None())?;
    }
    Ok(d)
}

/// Write a minimal P1 CIF from Cartesian coordinates in Å.
///
/// ``cell`` is None or a dict with all six keys: a, b, c, alpha, beta, gamma
/// (lengths in Å, angles in degrees). Without a cell, writes Cartesian columns.
/// Raises ValueError for invalid cells, coordinate counts, non-finite values
/// or conversion overflow. This projects atomic positions only: bonds,
/// occupancies, symmetry operations and other metadata are not serialized.
#[pyfunction]
#[pyo3(signature = (mol, coords, cell=None))]
pub(super) fn write_cif_checked(
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    cell: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let cell = cell
        .map(|d| -> PyResult<chematic_mol::UnitCell> {
            Ok(chematic_mol::UnitCell {
                a: dict_get_required(d, "a")?,
                b: dict_get_required(d, "b")?,
                c: dict_get_required(d, "c")?,
                alpha: dict_get_required(d, "alpha")?,
                beta: dict_get_required(d, "beta")?,
                gamma: dict_get_required(d, "gamma")?,
            })
        })
        .transpose()?;
    let coords: Vec<_> = coords.into_iter().map(|[x, y, z]| (x, y, z)).collect();
    chematic_mol::write_cif_checked(&mol.inner, &coords, cell.as_ref())
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
