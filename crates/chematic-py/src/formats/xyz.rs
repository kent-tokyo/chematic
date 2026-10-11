use super::*;

/// Parse an XYZ string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances.
///
///     mol, coords = chematic.from_xyz(open("molecule.xyz").read())
#[pyfunction]
pub(super) fn from_xyz(xyz_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, c3d) =
        chematic_3d::parse_xyz(xyz_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let coords = c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
        coords,
    ))
}

/// Build the Python dict returned by [`from_extxyz`]/[`from_extxyz_all`] for
/// one parsed [`chematic_mol::XyzFrame`].
pub(crate) fn extxyz_frame_to_pydict<'py>(
    py: Python<'py>,
    frame: &chematic_mol::XyzFrame,
) -> PyResult<Bound<'py, PyDict>> {
    let mol = Mol {
        inner: Arc::new(frame.to_molecule()),
        props: Default::default(),
    };
    let coords: Vec<Vec<f64>> = frame
        .coords()
        .into_iter()
        .map(|(x, y, z)| vec![x, y, z])
        .collect();

    let d = PyDict::new(py);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    match frame.lattice {
        Some(l) => d.set_item("lattice", l.to_vec())?,
        None => d.set_item("lattice", py.None())?,
    }

    let properties = PyDict::new(py);
    for prop in &frame.properties {
        match prop.kind {
            chematic_mol::XyzPropertyKind::Real => {
                let v: Vec<Vec<f64>> = extxyz_property_rows(prop, |x| match x {
                    chematic_mol::XyzValue::Real(r) => *r,
                    _ => unreachable!("XyzProperty.kind invariant: all rows share prop.kind"),
                });
                properties.set_item(&prop.name, v)?;
            }
            chematic_mol::XyzPropertyKind::Integer => {
                let v: Vec<Vec<i64>> = extxyz_property_rows(prop, |x| match x {
                    chematic_mol::XyzValue::Integer(i) => *i,
                    _ => unreachable!("XyzProperty.kind invariant: all rows share prop.kind"),
                });
                properties.set_item(&prop.name, v)?;
            }
            chematic_mol::XyzPropertyKind::String => {
                let v: Vec<Vec<String>> = extxyz_property_rows(prop, |x| match x {
                    chematic_mol::XyzValue::Str(s) => s.clone(),
                    _ => unreachable!("XyzProperty.kind invariant: all rows share prop.kind"),
                });
                properties.set_item(&prop.name, v)?;
            }
            chematic_mol::XyzPropertyKind::Logical => {
                let v: Vec<Vec<bool>> = extxyz_property_rows(prop, |x| match x {
                    chematic_mol::XyzValue::Logical(b) => *b,
                    _ => unreachable!("XyzProperty.kind invariant: all rows share prop.kind"),
                });
                properties.set_item(&prop.name, v)?;
            }
        }
    }
    d.set_item("properties", properties)?;

    let info = PyDict::new(py);
    for (k, v) in &frame.info {
        info.set_item(k, v)?;
    }
    d.set_item("info", info)?;

    Ok(d)
}

pub(super) fn extxyz_property_rows<T>(
    prop: &chematic_mol::XyzProperty,
    extract: impl Fn(&chematic_mol::XyzValue) -> T,
) -> Vec<Vec<T>> {
    prop.values
        .iter()
        .map(|row| row.iter().map(&extract).collect())
        .collect()
}

/// Parse an Extended XYZ (extxyz) frame -- ASE's cell/per-atom-property
/// superset of plain XYZ -- and return a dict describing it. A plain XYZ
/// file (free-form comment, no ``Lattice=``/``Properties=``) parses too,
/// with ``lattice=None`` and empty ``properties``/``info``.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]],
///     "lattice": list[float] | None (9 numbers, row-major cell matrix),
///     "properties": dict[str, list[list[float | int | str | bool]]]
///     (per-atom columns beyond position, e.g. ``"forces"``, ``"charge"``),
///     "info": dict[str, str]}`` (other frame metadata, e.g. ``"energy"``).
///
/// Raises:
///     ValueError: on malformed input (bad ``Lattice=``/``Properties=``,
///     wrong atom-row column count, non-finite values, ...).
///
/// Example::
///
///     result = chematic.from_extxyz(open("frame.xyz").read())
///     forces = result["properties"].get("forces")
#[pyfunction]
pub(super) fn from_extxyz<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let frame =
        chematic_mol::parse_extxyz(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    extxyz_frame_to_pydict(py, &frame)
}

/// Parse every frame of a multi-frame extxyz trajectory. See
/// :func:`from_extxyz` for the shape of each returned dict.
///
/// Raises:
///     ValueError: on the first parse failure (stops there).
#[pyfunction]
pub(super) fn from_extxyz_all<'py>(
    py: Python<'py>,
    text: &str,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let frames =
        chematic_mol::parse_extxyz_all(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    frames
        .iter()
        .map(|f| extxyz_frame_to_pydict(py, f))
        .collect()
}

/// Write a molecule + coordinates as an Extended XYZ (extxyz) frame.
///
/// Args:
///     mol: The molecule to write; atom order is preserved as the extxyz
///         atom order.
///     coords: Cartesian coordinates as ``[[x, y, z], ...]`` (Å), same
///         order and length as `mol`'s atoms.
///     lattice: Optional 9-number row-major cell matrix
///         (``[ax, ay, az, bx, by, bz, cx, cy, cz]``).
///     properties: Optional ``dict[str, list[list[float]]]`` of extra
///         real-valued per-atom columns (e.g.
///         ``{"forces": [[fx, fy, fz], ...]}``), one row per atom, same
///         row length within a column. Only real-valued (``R``) columns
///         are supported from Python; build a
///         ``chematic_mol::XyzProperty`` directly from Rust for
///         integer/string/logical columns.
///     info: Optional ``dict[str, str]`` of extra frame metadata (e.g.
///         ``{"energy": "-76.4"}``), written in dict iteration order.
///
/// Returns:
///     str: extxyz file contents (one frame).
///
/// Raises:
///     ValueError: if `coords`' length doesn't match `mol`'s atom count,
///     or if a `properties` column's row count doesn't match it.
#[pyfunction]
#[pyo3(signature = (mol, coords, lattice=None, properties=None, info=None))]
pub(super) fn to_extxyz<'py>(
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    lattice: Option<[f64; 9]>,
    properties: Option<Bound<'py, PyDict>>,
    info: Option<Bound<'py, PyDict>>,
) -> PyResult<String> {
    let n = mol.inner.atom_count();
    if coords.len() != n {
        return Err(PyValueError::new_err(format!(
            "coords has {} row(s), mol has {n} atom(s)",
            coords.len()
        )));
    }
    let atoms: Vec<chematic_mol::XyzAtom> = (0..n)
        .map(|i| {
            let element = mol.inner.atom(chematic_core::AtomIdx(i as u32)).element;
            let [x, y, z] = coords[i];
            chematic_mol::XyzAtom { element, x, y, z }
        })
        .collect();

    let mut xyz_properties = Vec::new();
    if let Some(properties) = properties {
        for (key, value) in properties.iter() {
            let name: String = key.extract()?;
            let rows: Vec<Vec<f64>> = value.extract().map_err(|_| {
                PyValueError::new_err(format!(
                    "properties['{name}'] must be a list of lists of float"
                ))
            })?;
            if rows.len() != n {
                return Err(PyValueError::new_err(format!(
                    "properties['{name}'] has {} row(s), mol has {n} atom(s)",
                    rows.len()
                )));
            }
            let count = rows.first().map_or(0, |r| r.len());
            if rows.iter().any(|r| r.len() != count) {
                return Err(PyValueError::new_err(format!(
                    "properties['{name}'] rows have inconsistent lengths"
                )));
            }
            let values = rows
                .into_iter()
                .map(|r| r.into_iter().map(chematic_mol::XyzValue::Real).collect())
                .collect();
            xyz_properties.push(chematic_mol::XyzProperty {
                name,
                kind: chematic_mol::XyzPropertyKind::Real,
                count,
                values,
            });
        }
    }

    let mut xyz_info = Vec::new();
    if let Some(info) = info {
        for (key, value) in info.iter() {
            xyz_info.push((key.extract()?, value.extract()?));
        }
    }

    let frame = chematic_mol::XyzFrame {
        atoms,
        comment: String::new(),
        lattice,
        properties: xyz_properties,
        info: xyz_info,
    };
    chematic_mol::write_extxyz(&frame).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Hill-notation molecular formula string into an element count dictionary.
///
/// Returns a ``dict[str, int]`` mapping element symbol → atom count.
/// Mirrors the API of PyPI libraries **chemparse** and **chemformula**.
///
/// Supported syntax:
///   - Simple formulas: ``"H2O"``, ``"C6H12O6"``
///   - Parentheses with multipliers: ``"Ca(OH)2"`` → ``{"Ca":1,"O":2,"H":2}``
///   - SMILES-style brackets: ``"[NH4]+"`` → ``{"N":1,"H":4}``
///   - Trailing charge signs are ignored: ``"NH4+"`` → same as ``"NH4"``
///
/// Raises:
///     ValueError: on empty formula or unbalanced parentheses.
///
///     chematic.parse_formula("C6H12O6")  # {"C": 6, "H": 12, "O": 6}
///     chematic.parse_formula("Ca(OH)2")  # {"Ca": 1, "O": 2, "H": 2}
///     chematic.parse_formula("[NH4]+")   # {"N": 1, "H": 4}
#[pyfunction]
pub(super) fn parse_formula<'py>(formula: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let counts =
        chematic_chem::parse_formula(formula).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let d = PyDict::new(py);
    let mut sorted: Vec<(String, u32)> = counts.into_iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    for (elem, cnt) in sorted {
        d.set_item(elem, cnt)?;
    }
    Ok(d)
}
