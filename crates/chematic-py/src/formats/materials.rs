use super::*;

// ---------------------------------------------------------------------------
// Small dict-extraction helpers shared by mmCIF/PQR/ORCA below.
//
// QCSchema (further down) deliberately does *not* use these -- it routes
// through Python's own `json` module instead (see that section's docs).
// ---------------------------------------------------------------------------

pub(crate) fn dict_get<'py, T>(d: &Bound<'py, PyDict>, key: &str, default: T) -> PyResult<T>
where
    T: for<'a> FromPyObject<'a, 'py>,
{
    match d.get_item(key)? {
        Some(v) if !v.is_none() => v.extract::<T>().map_err(Into::into),
        _ => Ok(default),
    }
}

pub(crate) fn dict_get_opt<'py, T>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<Option<T>>
where
    T: for<'a> FromPyObject<'a, 'py>,
{
    match d.get_item(key)? {
        Some(v) if !v.is_none() => Ok(Some(v.extract::<T>().map_err(Into::into)?)),
        _ => Ok(None),
    }
}

pub(crate) fn dict_get_required<'py, T>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<T>
where
    T: for<'a> FromPyObject<'a, 'py>,
{
    match d.get_item(key)? {
        Some(v) if !v.is_none() => v.extract::<T>().map_err(Into::into),
        _ => Err(PyValueError::new_err(format!(
            "missing required key '{key}'"
        ))),
    }
}

/// A `str`-typed dict value that must be exactly one character (Python has
/// no separate `char` type). Rejects a multi-character string rather than
/// silently truncating it to its first character.
pub(crate) fn dict_get_opt_char(d: &Bound<PyDict>, key: &str) -> PyResult<Option<char>> {
    match dict_get_opt::<String>(d, key)? {
        None => Ok(None),
        Some(s) => {
            let mut chars = s.chars();
            let c = chars.next().ok_or_else(|| {
                PyValueError::new_err(format!(
                    "'{key}' must be a single character, got an empty string"
                ))
            })?;
            if chars.next().is_some() {
                return Err(PyValueError::new_err(format!(
                    "'{key}' must be a single character, got {s:?}"
                )));
            }
            Ok(Some(c))
        }
    }
}

pub(super) fn element_from_symbol(sym: &str) -> PyResult<chematic_core::Element> {
    chematic_core::Element::from_symbol(sym)
        .ok_or_else(|| PyValueError::new_err(format!("unknown element symbol: {sym:?}")))
}

// ---------------------------------------------------------------------------
// mmCIF
// ---------------------------------------------------------------------------

pub(super) fn mmcif_atom_to_pydict<'py>(
    py: Python<'py>,
    a: &chematic_mol::MmcifAtomRecord,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("group_pdb", &a.group_pdb)?;
    d.set_item("serial", a.serial)?;
    d.set_item("element", a.element.symbol())?;
    d.set_item("atom_name", &a.atom_name)?;
    d.set_item("alt_loc", a.alt_loc.map(|c| c.to_string()))?;
    d.set_item("res_name", &a.res_name)?;
    d.set_item("chain_id", &a.chain_id)?;
    d.set_item("res_seq", a.res_seq)?;
    d.set_item("label_seq_id", a.label_seq_id)?;
    d.set_item("icode", a.icode.map(|c| c.to_string()))?;
    d.set_item("x", a.x)?;
    d.set_item("y", a.y)?;
    d.set_item("z", a.z)?;
    d.set_item("occupancy", a.occupancy)?;
    d.set_item("b_iso", a.b_iso)?;
    d.set_item("formal_charge", a.formal_charge)?;
    d.set_item("entity_id", &a.entity_id)?;
    d.set_item("model_num", a.model_num)?;
    Ok(d)
}

pub(super) fn pydict_to_mmcif_atom(d: &Bound<PyDict>) -> PyResult<chematic_mol::MmcifAtomRecord> {
    let element: String = dict_get_required(d, "element")?;
    Ok(chematic_mol::MmcifAtomRecord {
        group_pdb: dict_get(d, "group_pdb", "ATOM".to_string())?,
        serial: dict_get(d, "serial", 1i64)?,
        element: element_from_symbol(&element)?,
        atom_name: dict_get(d, "atom_name", String::new())?,
        alt_loc: dict_get_opt_char(d, "alt_loc")?,
        res_name: dict_get(d, "res_name", "UNL".to_string())?,
        chain_id: dict_get(d, "chain_id", "A".to_string())?,
        res_seq: dict_get(d, "res_seq", 1i64)?,
        label_seq_id: dict_get_opt(d, "label_seq_id")?,
        icode: dict_get_opt_char(d, "icode")?,
        x: dict_get_required(d, "x")?,
        y: dict_get_required(d, "y")?,
        z: dict_get_required(d, "z")?,
        occupancy: dict_get(d, "occupancy", 1.0)?,
        b_iso: dict_get(d, "b_iso", 0.0)?,
        formal_charge: dict_get_opt(d, "formal_charge")?,
        entity_id: dict_get_opt(d, "entity_id")?,
        model_num: dict_get(d, "model_num", 1i32)?,
    })
}

pub(super) fn unit_cell_to_pydict<'py>(
    py: Python<'py>,
    cell: &chematic_mol::UnitCell,
) -> PyResult<Bound<'py, PyDict>> {
    let cd = PyDict::new(py);
    cd.set_item("a", cell.a)?;
    cd.set_item("b", cell.b)?;
    cd.set_item("c", cell.c)?;
    cd.set_item("alpha", cell.alpha)?;
    cd.set_item("beta", cell.beta)?;
    cd.set_item("gamma", cell.gamma)?;
    Ok(cd)
}

pub(super) fn pydict_to_unit_cell(d: &Bound<PyDict>) -> PyResult<chematic_mol::UnitCell> {
    Ok(chematic_mol::UnitCell {
        a: dict_get_required(d, "a")?,
        b: dict_get_required(d, "b")?,
        c: dict_get_required(d, "c")?,
        alpha: dict_get_required(d, "alpha")?,
        beta: dict_get_required(d, "beta")?,
        gamma: dict_get_required(d, "gamma")?,
    })
}

fn protein_projection_json_value(
    projection: &chematic_mol::ProteinProjection,
) -> serde_json::Value {
    serde_json::json!({
        "model_count": projection.model_count(),
        "chain_count": projection.chain_count(),
        "residue_count": projection.residue_count(),
        "atom_count": projection.atom_count(),
        "chains": projection.chains.iter().map(|chain| serde_json::json!({
            "model_number": chain.model_number,
            "id": chain.id,
            "sequence": chain.sequence(),
            "residues": chain.residues.iter().map(|residue| serde_json::json!({
                "name": residue.name,
                "sequence_number": residue.sequence_number,
                "insertion_code": residue.insertion_code.map(|value| value.to_string()),
                "canonical_one_letter_code": residue.canonical_one_letter_code.to_string(),
                "modified": residue.modified,
                "atoms": residue.atoms.iter().map(|atom| serde_json::json!({
                    "serial": atom.serial,
                    "name": atom.name,
                    "element": atom.element.symbol(),
                    "position": atom.position,
                    "occupancy": atom.occupancy,
                    "b_iso": atom.b_iso,
                    "alt_loc": atom.alt_loc.map(|value| value.to_string()),
                })).collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
    })
}

/// Project PDB text to an amino-acid-only model/chain/residue hierarchy.
///
/// Waters, ligands, ions, and nucleic acids are intentionally excluded. The
/// result is JSON so the same stable schema can be consumed by Python and
/// WebAssembly clients.
#[pyfunction]
pub(super) fn protein_from_pdb_json(text: &str) -> PyResult<String> {
    let projection = chematic_mol::project_protein_from_pdb(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    serde_json::to_string(&protein_projection_json_value(&projection))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Project mmCIF text to the same amino-acid-only hierarchy as
/// :func:`protein_from_pdb_json`.
#[pyfunction]
pub(super) fn protein_from_mmcif_json(text: &str) -> PyResult<String> {
    let projection = chematic_mol::project_protein_from_mmcif(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    serde_json::to_string(&protein_projection_json_value(&projection))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Parse an mmCIF (macromolecular CIF, `_atom_site.*` loop) string.
///
/// Unlike :func:`parse_cif` (small-molecule CIF), this reads the
/// `_atom_site` loop's full per-atom detail (residue/chain/occupancy/
/// B-factor/model number) rather than just element+position.
///
/// Args:
///     text: mmCIF file contents.
///     max_input_bytes: byte-size limit (default 128 MiB).
///     max_atoms: `_atom_site` row-count limit, across all models (default
///         2,000,000).
///     max_line_len: per-line byte limit (default 8192).
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "atoms":
///     list[dict], "cell": dict | None, "space_group": str | None,
///     "unhandled_columns": list[str]}``. ``mol``/``coords`` include every
///     model's atoms (no bonds -- mmCIF carries no bond table); filter
///     ``atoms`` by its ``"model_num"`` key first if only one model is
///     wanted. Each atom dict has keys ``group_pdb``, ``serial``,
///     ``element``, ``atom_name``, ``alt_loc``, ``res_name``, ``chain_id``,
///     ``res_seq``, ``label_seq_id``, ``icode``, ``x``, ``y``, ``z``,
///     ``occupancy``, ``b_iso``, ``formal_charge``, ``entity_id``,
///     ``model_num``.
///
/// Raises:
///     ValueError: on parse failure (no ``_atom_site`` loop, missing
///     required column, unresolvable element, non-finite number, or a
///     limit exceeded).
///
/// Example::
///
///     result = chematic.parse_mmcif(open("structure.cif").read())
///     print(result["mol"].formula, len(result["atoms"]))
#[pyfunction]
#[pyo3(signature = (text, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_mmcif<'py>(
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
    let result = chematic_mol::parse_mmcif_with_limits(text, &limits).map_err(|e| {
        use chematic_mol::MmcifError as E;
        let code = match &e {
            E::InputTooLarge { .. } => Some("input_too_large"),
            E::LineTooLong { .. } => Some("line_too_long"),
            E::TooManyAtoms { .. } => Some("too_many_atoms"),
            _ => None,
        };
        match code {
            Some(code) => {
                crate::errors::input_error("resource_limit", code, "mmcif", e.to_string())
            }
            None => crate::errors::malformed("mmcif", e.to_string()),
        }
    })?;
    let (mol, coords) = result.to_molecule();

    let d = PyDict::new(py);
    d.set_item(
        "mol",
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
    )?;
    let py_coords: Vec<Vec<f64>> = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
    d.set_item("coords", py_coords)?;
    let atoms: Vec<Bound<PyDict>> = result
        .atoms
        .iter()
        .map(|a| mmcif_atom_to_pydict(py, a))
        .collect::<PyResult<_>>()?;
    d.set_item("atoms", atoms)?;
    match &result.cell {
        Some(cell) => d.set_item("cell", unit_cell_to_pydict(py, cell)?)?,
        None => d.set_item("cell", py.None())?,
    }
    d.set_item("space_group", &result.space_group)?;
    d.set_item("unhandled_columns", &result.unhandled_columns)?;
    Ok(d)
}

/// Write an mmCIF file from atom records.
///
/// Args:
///     atoms: list of dicts, same shape as :func:`parse_mmcif`'s
///         ``"atoms"`` entries. Every key is optional except ``element``,
///         ``x``, ``y``, ``z`` -- see :func:`parse_mmcif` for the full key
///         list and their defaults (``group_pdb`` defaults to ``"ATOM"``,
///         ``chain_id`` to ``"A"``, ``occupancy`` to ``1.0``, etc.).
///     cell: optional ``{"a", "b", "c", "alpha", "beta", "gamma"}`` dict,
///         same shape as :func:`parse_cif`'s ``"cell"``.
///     space_group: optional space-group name (``_symmetry.space_group_name_H-M``).
///     data_block_name: the CIF ``data_<name>`` block name (default
///         ``"chematic"``).
///
/// Raises:
///     ValueError: on an unknown element symbol or a malformed
///         ``alt_loc``/``icode`` (must be a single character).
#[pyfunction]
#[pyo3(signature = (atoms, cell=None, space_group=None, data_block_name="chematic"))]
pub(super) fn write_mmcif(
    atoms: Vec<Bound<PyDict>>,
    cell: Option<Bound<PyDict>>,
    space_group: Option<&str>,
    data_block_name: &str,
) -> PyResult<String> {
    const MAX_ATOMS: usize = 100_000;
    if atoms.len() > MAX_ATOMS {
        return Err(PyValueError::new_err(format!(
            "mmCIF output exceeds maximum atom count ({MAX_ATOMS})"
        )));
    }
    let records: Vec<chematic_mol::MmcifAtomRecord> = atoms
        .iter()
        .map(pydict_to_mmcif_atom)
        .collect::<PyResult<_>>()?;
    let cell = cell.as_ref().map(pydict_to_unit_cell).transpose()?;
    Ok(chematic_mol::write_mmcif(
        &records,
        cell.as_ref(),
        space_group,
        data_block_name,
    ))
}

// ---------------------------------------------------------------------------
// PQR
// ---------------------------------------------------------------------------

pub(super) fn pqr_atom_to_pydict<'py>(
    py: Python<'py>,
    a: &chematic_mol::PqrAtomRecord,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("group_pdb", &a.group_pdb)?;
    d.set_item("serial", a.serial)?;
    d.set_item("atom_name", &a.atom_name)?;
    d.set_item("res_name", &a.res_name)?;
    d.set_item("chain_id", &a.chain_id)?;
    d.set_item("res_seq", a.res_seq)?;
    d.set_item("icode", a.icode.map(|c| c.to_string()))?;
    d.set_item("x", a.x)?;
    d.set_item("y", a.y)?;
    d.set_item("z", a.z)?;
    d.set_item("charge", a.charge)?;
    d.set_item("radius", a.radius)?;
    d.set_item("element", a.element.symbol())?;
    Ok(d)
}

/// Build a [`chematic_mol::PqrAtomRecord`] from a Python dict. `element` is
/// optional -- when omitted, it is inferred the same way [`parse_pqr`]
/// infers it for a real PQR file, via [`chematic_mol::infer_element`],
/// since [`chematic_mol::write_pqr`] never reads a record's `element` field
/// back out anyway (PQR has no element column at all -- see that module's
/// docs). Callers shouldn't have to supply a field the format itself
/// doesn't store.
pub(super) fn pydict_to_pqr_atom(d: &Bound<PyDict>) -> PyResult<chematic_mol::PqrAtomRecord> {
    let group_pdb: String = dict_get(d, "group_pdb", "ATOM".to_string())?;
    let atom_name: String = dict_get(d, "atom_name", String::new())?;
    let res_name: String = dict_get(d, "res_name", "UNL".to_string())?;
    let element = match dict_get_opt::<String>(d, "element")? {
        Some(sym) => element_from_symbol(&sym)?,
        None => chematic_mol::infer_element(&group_pdb, &res_name, &atom_name).ok_or_else(|| {
            PyValueError::new_err(format!(
                "could not infer an element from atom name {atom_name:?} (pass element= explicitly)"
            ))
        })?,
    };
    Ok(chematic_mol::PqrAtomRecord {
        group_pdb,
        serial: dict_get(d, "serial", 1i64)?,
        atom_name,
        res_name,
        chain_id: dict_get_opt(d, "chain_id")?,
        res_seq: dict_get(d, "res_seq", 1i64)?,
        icode: dict_get_opt_char(d, "icode")?,
        x: dict_get_required(d, "x")?,
        y: dict_get_required(d, "y")?,
        z: dict_get_required(d, "z")?,
        charge: dict_get(d, "charge", 0.0)?,
        radius: dict_get(d, "radius", 0.0)?,
        element,
    })
}

/// Parse a PQR (PDB-like ATOM/HETATM records with per-atom charge/radius)
/// string.
///
/// Args:
///     text: PQR file contents.
///     max_input_bytes: byte-size limit (default 64 MiB).
///     max_atoms: atom-count limit (default 2,000,000).
///     max_line_len: per-line byte limit (default 1024).
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "atoms":
///     list[dict]}`` (no bonds -- PQR carries no connectivity). Each atom
///     dict has keys ``group_pdb``, ``serial``, ``atom_name``, ``res_name``,
///     ``chain_id`` (``None`` if the file omits the chain column), ``res_seq``,
///     ``icode``, ``x``, ``y``, ``z``, ``charge``, ``radius``, ``element``
///     (inferred from the atom name -- see :func:`infer_element`).
///
/// Raises:
///     ValueError: on parse failure (no ATOM/HETATM records, wrong field
///     count, unresolvable element, non-finite value, or a limit exceeded).
///
/// Example::
///
///     result = chematic.parse_pqr(open("structure.pqr").read())
///     charges = [a["charge"] for a in result["atoms"]]
#[pyfunction]
#[pyo3(signature = (text, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_pqr<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::PqrParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let result = chematic_mol::parse_pqr_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (mol, coords) = result.to_molecule();

    let d = PyDict::new(py);
    d.set_item(
        "mol",
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
    )?;
    let py_coords: Vec<Vec<f64>> = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
    d.set_item("coords", py_coords)?;
    let atoms: Vec<Bound<PyDict>> = result
        .atoms
        .iter()
        .map(|a| pqr_atom_to_pydict(py, a))
        .collect::<PyResult<_>>()?;
    d.set_item("atoms", atoms)?;
    Ok(d)
}

/// Write a PQR file from atom records.
///
/// Args:
///     atoms: list of dicts, same shape as :func:`parse_pqr`'s ``"atoms"``
///         entries. Every key is optional except ``x``, ``y``, ``z`` --
///         ``element`` is inferred via :func:`infer_element` when omitted
///         (PQR itself has no element column, so this never affects the
///         written file).
///
/// Example::
///
///     text = chematic.write_pqr([
///         {"atom_name": "N", "res_name": "ALA", "res_seq": 1,
///          "x": -0.966, "y": 1.523, "z": 1.412, "charge": -0.4, "radius": 1.5},
///     ])
#[pyfunction]
pub(super) fn write_pqr(atoms: Vec<Bound<PyDict>>) -> PyResult<String> {
    const MAX_ATOMS: usize = 100_000;
    if atoms.len() > MAX_ATOMS {
        return Err(PyValueError::new_err(format!(
            "PQR output exceeds maximum atom count ({MAX_ATOMS})"
        )));
    }
    let records: Vec<chematic_mol::PqrAtomRecord> = atoms
        .iter()
        .map(pydict_to_pqr_atom)
        .collect::<PyResult<_>>()?;
    Ok(chematic_mol::write_pqr(&records))
}

/// Infer an element symbol from a PQR/PDB atom name (PQR has no element
/// column of its own -- see module docs on ``chematic_mol::pqr``).
///
/// Returns:
///     str | None: the element symbol, or ``None`` if no element could be
///     inferred (e.g. an atom name with no alphabetic characters at all).
///
/// Example::
///
///     chematic.infer_element("ATOM", "ALA", "CA")     # "C" (not "Ca")
///     chematic.infer_element("HETATM", "ZN", "ZN")     # "Zn"
#[pyfunction]
pub(super) fn infer_element(group_pdb: &str, res_name: &str, atom_name: &str) -> Option<String> {
    chematic_mol::infer_element(group_pdb, res_name, atom_name).map(|e| e.symbol().to_string())
}

// ---------------------------------------------------------------------------
// ORCA input
// ---------------------------------------------------------------------------

pub(super) fn orca_block_to_pydict<'py>(
    py: Python<'py>,
    b: &chematic_mol::OrcaBlock,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("name", &b.name)?;
    d.set_item("raw", &b.raw)?;
    d.set_item("has_end", b.has_end)?;
    Ok(d)
}

pub(super) fn pydict_to_orca_block(d: &Bound<PyDict>) -> PyResult<chematic_mol::OrcaBlock> {
    Ok(chematic_mol::OrcaBlock {
        name: dict_get_required(d, "name")?,
        raw: dict_get(d, "raw", String::new())?,
        has_end: dict_get(d, "has_end", true)?,
    })
}

pub(super) fn orca_atom_to_pydict<'py>(
    py: Python<'py>,
    a: &chematic_mol::OrcaAtom,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("element", a.element.symbol())?;
    d.set_item("x", a.x)?;
    d.set_item("y", a.y)?;
    d.set_item("z", a.z)?;
    d.set_item("frozen", a.frozen.to_vec())?;
    d.set_item("extra", &a.extra)?;
    Ok(d)
}

pub(super) fn pydict_to_orca_atom(d: &Bound<PyDict>) -> PyResult<chematic_mol::OrcaAtom> {
    let sym: String = dict_get_required(d, "element")?;
    let frozen: Vec<bool> = dict_get(d, "frozen", vec![false, false, false])?;
    if frozen.len() != 3 {
        return Err(PyValueError::new_err(
            "'frozen' must have exactly 3 entries (x, y, z)",
        ));
    }
    Ok(chematic_mol::OrcaAtom {
        element: element_from_symbol(&sym)?,
        x: dict_get_required(d, "x")?,
        y: dict_get_required(d, "y")?,
        z: dict_get_required(d, "z")?,
        frozen: [frozen[0], frozen[1], frozen[2]],
        extra: dict_get_opt(d, "extra")?,
    })
}

pub(super) fn orca_coords_to_pydict<'py>(
    py: Python<'py>,
    c: &chematic_mol::OrcaCoords,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    match c {
        chematic_mol::OrcaCoords::Xyz {
            charge,
            multiplicity,
            atoms,
        } => {
            d.set_item("kind", "xyz")?;
            d.set_item("charge", charge)?;
            d.set_item("multiplicity", multiplicity)?;
            let py_atoms: Vec<Bound<PyDict>> = atoms
                .iter()
                .map(|a| orca_atom_to_pydict(py, a))
                .collect::<PyResult<_>>()?;
            d.set_item("atoms", py_atoms)?;
        }
        chematic_mol::OrcaCoords::XyzFile {
            charge,
            multiplicity,
            filename,
        } => {
            d.set_item("kind", "xyzfile")?;
            d.set_item("charge", charge)?;
            d.set_item("multiplicity", multiplicity)?;
            d.set_item("filename", filename)?;
        }
        chematic_mol::OrcaCoords::GzmtFile {
            charge,
            multiplicity,
            filename,
        } => {
            d.set_item("kind", "gzmtfile")?;
            d.set_item("charge", charge)?;
            d.set_item("multiplicity", multiplicity)?;
            d.set_item("filename", filename)?;
        }
        chematic_mol::OrcaCoords::Internal {
            charge,
            multiplicity,
            raw,
        } => {
            d.set_item("kind", "int")?;
            d.set_item("charge", charge)?;
            d.set_item("multiplicity", multiplicity)?;
            d.set_item("raw", raw)?;
        }
    }
    Ok(d)
}

pub(super) fn pydict_to_orca_coords(d: &Bound<PyDict>) -> PyResult<chematic_mol::OrcaCoords> {
    let kind: String = dict_get_required(d, "kind")?;
    let charge: i32 = dict_get(d, "charge", 0)?;
    let multiplicity: u32 = dict_get(d, "multiplicity", 1)?;
    match kind.as_str() {
        "xyz" => {
            let atoms_list: Vec<Bound<PyDict>> = dict_get(d, "atoms", Vec::new())?;
            let atoms = atoms_list
                .iter()
                .map(pydict_to_orca_atom)
                .collect::<PyResult<_>>()?;
            Ok(chematic_mol::OrcaCoords::Xyz {
                charge,
                multiplicity,
                atoms,
            })
        }
        "xyzfile" => Ok(chematic_mol::OrcaCoords::XyzFile {
            charge,
            multiplicity,
            filename: dict_get_required(d, "filename")?,
        }),
        "gzmtfile" => Ok(chematic_mol::OrcaCoords::GzmtFile {
            charge,
            multiplicity,
            filename: dict_get_required(d, "filename")?,
        }),
        "int" => Ok(chematic_mol::OrcaCoords::Internal {
            charge,
            multiplicity,
            raw: dict_get(d, "raw", String::new())?,
        }),
        other => Err(PyValueError::new_err(format!(
            "unknown coords 'kind' {other:?}, expected 'xyz'/'xyzfile'/'gzmtfile'/'int'"
        ))),
    }
}

/// Parse an ORCA input file (`.inp`).
///
/// Returns:
///     dict: ``{"comments": list[str], "keywords": list[str], "blocks":
///     list[dict], "coords": dict | None}``.
///
///     Each block dict is ``{"name": str, "raw": str, "has_end": bool}``
///     (``%name ... end`` blocks; ``has_end=False`` for a single-line
///     directive like ``%maxcore 3000``).
///
///     ``coords`` (if the input has a ``* ... *`` block) is tagged by a
///     ``"kind"`` key:
///
///     - ``"xyz"``: ``{"kind": "xyz", "charge": int, "multiplicity": int,
///       "atoms": list[dict]}``, each atom dict
///       ``{"element": str, "x": float, "y": float, "z": float,
///       "frozen": [bool, bool, bool], "extra": str | None}``.
///     - ``"xyzfile"`` / ``"gzmtfile"``: ``{"kind": ..., "charge": int,
///       "multiplicity": int, "filename": str}`` (external geometry file).
///     - ``"int"``: ``{"kind": "int", "charge": int, "multiplicity": int,
///       "raw": str}`` (internal/Z-matrix coordinates, preserved verbatim
///       -- not semantically parsed).
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     result = chematic.parse_orca_input(open("job.inp").read())
///     print(result["keywords"])
#[pyfunction]
pub(super) fn parse_orca_input<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let input =
        chematic_mol::parse_orca_input(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let chematic_mol::OrcaInput {
        comments,
        keywords,
        blocks,
        coords,
    } = input;

    let d = PyDict::new(py);
    d.set_item("comments", comments)?;
    d.set_item("keywords", keywords)?;
    let py_blocks: Vec<Bound<PyDict>> = blocks
        .iter()
        .map(|b| orca_block_to_pydict(py, b))
        .collect::<PyResult<_>>()?;
    d.set_item("blocks", py_blocks)?;
    match &coords {
        Some(c) => d.set_item("coords", orca_coords_to_pydict(py, c)?)?,
        None => d.set_item("coords", py.None())?,
    }
    Ok(d)
}

/// Write an ORCA input file (`.inp`) from a dict, same shape as
/// :func:`parse_orca_input`'s return value (every key optional -- an empty
/// dict writes an empty file).
///
/// Example::
///
///     text = chematic.write_orca_input({
///         "keywords": ["B3LYP", "def2-SVP", "Opt"],
///         "coords": {"kind": "xyz", "charge": 0, "multiplicity": 1,
///                    "atoms": [{"element": "O", "x": 0.0, "y": 0.0, "z": 0.0}]},
///     })
#[pyfunction]
pub(super) fn write_orca_input(input: &Bound<PyDict>) -> PyResult<String> {
    let comments: Vec<String> = dict_get(input, "comments", Vec::new())?;
    let keywords: Vec<String> = dict_get(input, "keywords", Vec::new())?;
    let blocks_list: Vec<Bound<PyDict>> = dict_get(input, "blocks", Vec::new())?;
    let blocks = blocks_list
        .iter()
        .map(pydict_to_orca_block)
        .collect::<PyResult<_>>()?;
    let coords_dict: Option<Bound<PyDict>> = dict_get_opt(input, "coords")?;
    let coords = coords_dict
        .as_ref()
        .map(pydict_to_orca_coords)
        .transpose()?;
    let orca_input = chematic_mol::OrcaInput {
        comments,
        keywords,
        blocks,
        coords,
    };
    Ok(chematic_mol::write_orca_input(&orca_input))
}

// ---------------------------------------------------------------------------
// ORCA output
// ---------------------------------------------------------------------------

/// Parse an ORCA output file (`.out` / `.log`).
///
/// Never raises on truncated/crashed output -- see ``"termination"``.
///
/// Returns:
///     dict: ``{"charge": int | None, "multiplicity": int | None,
///     "final_energy_hartree": float | None, "trajectory": list[dict],
///     "frequencies_cm1": list[float], "termination": dict,
///     "optimization_convergence": str}``.
///
///     ``trajectory`` is one ``{"mol": Mol, "coords": list[list[float]]}``
///     dict per ``CARTESIAN COORDINATES (ANGSTROEM)`` block found, in file
///     order (the geometry trajectory for an optimization job, or a single
///     frame for a single-point/frequency job); ``trajectory[-1]`` is the
///     final geometry.
///
///     ``termination`` is ``{"kind": "normal" | "error" | "incomplete",
///     "message": str | None}`` (``message`` is the verbatim error line,
///     only set for ``"error"``). Per the ORCA manual, ``"normal"`` does
///     **not** by itself mean a requested geometry optimization converged
///     -- check ``"optimization_convergence"`` separately.
///
///     ``optimization_convergence`` is one of ``"not_requested"``
///     (not an optimization job), ``"converged"``, ``"not_converged"``, or
///     ``"unknown"`` (truncated mid-optimization).
///
/// Raises:
///     ValueError: on a non-finite energy/coordinate/frequency value, or
///     an oversized input.
///
/// Example::
///
///     result = chematic.parse_orca_output(open("job.out").read())
///     if result["termination"]["kind"] == "normal":
///         print(result["final_energy_hartree"])
#[pyfunction]
pub(super) fn parse_orca_output<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let output =
        chematic_mol::parse_orca_output(text).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let d = PyDict::new(py);
    d.set_item("charge", output.charge)?;
    d.set_item("multiplicity", output.multiplicity)?;
    d.set_item("final_energy_hartree", output.final_energy_hartree)?;

    let trajectory: Vec<Bound<PyDict>> = output
        .trajectory
        .into_iter()
        .map(|f| {
            let fd = PyDict::new(py);
            fd.set_item(
                "mol",
                Mol {
                    inner: Arc::new(f.mol),
                    props: Default::default(),
                },
            )?;
            let coords: Vec<Vec<f64>> = f.coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
            fd.set_item("coords", coords)?;
            Ok::<_, PyErr>(fd)
        })
        .collect::<PyResult<_>>()?;
    d.set_item("trajectory", trajectory)?;
    d.set_item("frequencies_cm1", output.frequencies_cm1)?;

    let (term_kind, term_msg): (&str, Option<String>) = match output.termination {
        chematic_mol::OrcaTermination::Normal => ("normal", None),
        chematic_mol::OrcaTermination::Error(msg) => ("error", Some(msg)),
        chematic_mol::OrcaTermination::Incomplete => ("incomplete", None),
    };
    let term_dict = PyDict::new(py);
    term_dict.set_item("kind", term_kind)?;
    term_dict.set_item("message", term_msg)?;
    d.set_item("termination", term_dict)?;

    let convergence = match output.optimization_convergence {
        chematic_mol::OrcaOptConvergence::NotRequested => "not_requested",
        chematic_mol::OrcaOptConvergence::Converged => "converged",
        chematic_mol::OrcaOptConvergence::NotConverged => "not_converged",
        chematic_mol::OrcaOptConvergence::Unknown => "unknown",
    };
    d.set_item("optimization_convergence", convergence)?;

    Ok(d)
}

// ---------------------------------------------------------------------------
// QCSchema (MolSSI Quantum Chemistry Schema)
//
// `chematic_mol`'s own `parse_*`/`write_*` functions already speak
// canonical JSON *text* end to end. Rather than hand-mapping every field of
// `QcMolecule`/`AtomicInput`/`AtomicResult` (23+ fields on `QcMolecule`
// alone, several of them open JSON bags -- `extras`/`unknown_fields`/
// `keywords`/`protocols`/`native_files`/`wavefunction` -- that are part of
// the spec's own extensibility design, see `qcschema.rs`'s module docs)
// onto bespoke `PyDict` field-by-field code, this binding routes through
// Python's own `json` module (stdlib) as the dict<->text boundary:
// `json.loads`/`json.dumps` on the canonical text `chematic_mol`'s own
// `write_*` functions already produce/accept. This is strictly *more*
// faithful than a hand-rolled mapping (an open bag field can never be
// silently dropped by an oversight it wasn't written to handle) and needs
// zero new dependencies -- `serde_json` is not a `chematic-py` dependency,
// and this deliberately never needs to name it (`chematic_mol::JsonObject`
// is never touched directly here).
//
// `write_*` strips the `"mol"`/`"coords"` convenience keys `parse_*` adds
// before re-serializing (a bare Python `Mol` object isn't JSON-serializable
// at all), so `chematic.write_qcschema_molecule(chematic.parse_qcschema_molecule(text))`
// round-trips without the caller needing to strip those themselves.
// ---------------------------------------------------------------------------

pub(super) fn json_loads<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let value = py.import("json")?.call_method1("loads", (text,))?;
    value
        .cast::<PyDict>()
        .map_err(|e| PyValueError::new_err(e.to_string()))
        .cloned()
}

/// `json.dumps` a dict, first stripping the `"mol"`/`"coords"` convenience
/// keys the `parse_*` functions in this section add (see section docs).
pub(super) fn json_dumps_stripped(d: &Bound<PyDict>) -> PyResult<String> {
    let py = d.py();
    let clean = d.copy()?;
    let _ = clean.del_item("mol");
    let _ = clean.del_item("coords");
    py.import("json")?
        .call_method1("dumps", (clean,))?
        .extract()
}

/// `(Mol, coords)` from an owned [`chematic_mol::ChematicMoleculeView`],
/// `coords` in Ångström as `[[x, y, z], ...]`.
pub(super) fn qc_view_into_py(view: chematic_mol::ChematicMoleculeView) -> (Mol, Vec<Vec<f64>>) {
    let coords = view
        .coords
        .points
        .iter()
        .map(|p| vec![p.x, p.y, p.z])
        .collect();
    let mol = Mol {
        inner: Arc::new(view.molecule),
        props: Default::default(),
    };
    (mol, coords)
}

/// Parse a QCSchema `Molecule` JSON document (`schema_name:
/// "qcschema_molecule"`).
///
/// Returns a dict with every QCSchema `Molecule` field (``symbols``,
/// ``geometry`` in Bohr, ``molecular_charge``, ``connectivity``, ...; see
/// <https://molssi.github.io/QCElemental/model_molecule.html> for the full
/// field reference) plus two chematic-side convenience keys:
///
/// - ``"mol"``: a :class:`Mol` built from ``symbols``/``connectivity`` (see
///   :func:`qc_molecule_to_chematic` for exactly what is gained/lost).
/// - ``"coords"``: the same geometry, converted to Ångström, as
///   ``[[x, y, z], ...]``.
///
/// :func:`write_qcschema_molecule` ignores/strips both convenience keys, so
/// this dict round-trips through it directly.
///
/// Raises:
///     ValueError: on malformed JSON, a schema violation, or (for the
///     ``"mol"``/``"coords"`` keys specifically) an unresolvable element
///     symbol or an out-of-range ``connectivity`` atom index.
///
/// Example::
///
///     result = chematic.parse_qcschema_molecule(open("water.json").read())
///     print(result["symbols"], result["mol"].formula)
#[pyfunction]
pub(super) fn parse_qcschema_molecule<'py>(
    py: Python<'py>,
    text: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let qc = chematic_mol::parse_qcschema_molecule(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let canonical = chematic_mol::write_qcschema_molecule(&qc);
    let d = json_loads(py, &canonical)?;
    let view = chematic_mol::qc_molecule_to_chematic(&qc)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (mol, coords) = qc_view_into_py(view);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    Ok(d)
}

/// Serialize a QCSchema `Molecule` dict (as returned by
/// :func:`parse_qcschema_molecule`, or hand-built with the same keys) to
/// canonical QCSchema JSON text.
///
/// Raises:
///     ValueError: if `molecule` doesn't validate as a QCSchema `Molecule`
///     document (missing required key, wrong type, non-finite number, ...).
#[pyfunction]
pub(super) fn write_qcschema_molecule(molecule: &Bound<PyDict>) -> PyResult<String> {
    let text = json_dumps_stripped(molecule)?;
    let qc = chematic_mol::parse_qcschema_molecule(&text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_qcschema_molecule(&qc))
}

/// Convert a chematic :class:`Mol` + coordinates into a QCSchema `Molecule`
/// dict (geometry in Bohr, per spec).
///
/// Args:
///     mol: the molecule (topology/bond orders).
///     coords: Cartesian coordinates (Å), ``[[x, y, z], ...]``, same order
///         and length as `mol`'s atoms.
///     molecular_charge: total charge (default ``0.0``).
///     molecular_multiplicity: spin multiplicity (default ``1``).
///
/// Raises:
///     ValueError: if `coords`' length doesn't match `mol`'s atom count.
///
/// Example::
///
///     qc = chematic.chematic_to_qc_molecule(mol, coords)
///     open("mol.json", "w").write(chematic.write_qcschema_molecule(qc))
#[pyfunction]
#[pyo3(signature = (mol, coords, molecular_charge=0.0, molecular_multiplicity=1))]
pub(super) fn chematic_to_qc_molecule<'py>(
    py: Python<'py>,
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    molecular_charge: f64,
    molecular_multiplicity: i64,
) -> PyResult<Bound<'py, PyDict>> {
    let c3d = flat_to_coords3d(&coords);
    let qc = chematic_mol::chematic_to_qc_molecule(
        &mol.inner,
        &c3d,
        molecular_charge,
        molecular_multiplicity,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let text = chematic_mol::write_qcschema_molecule(&qc);
    json_loads(py, &text)
}

/// Convert a QCSchema `Molecule` dict into a chematic ``(Mol, coords,
/// molecular_charge, molecular_multiplicity)`` tuple.
///
/// ``coords`` is Ångström, ``[[x, y, z], ...]``. See
/// ``chematic_mol::qc_molecule_to_chematic``'s Rust docs for exactly what
/// is gained/lost by this conversion (QCSchema has no aromaticity/stereo
/// concept; ``connectivity`` bond orders map onto the nearest chematic
/// bond order).
///
/// Raises:
///     ValueError: unresolvable element symbol, or a ``connectivity`` entry
///     referencing an out-of-range atom index.
#[pyfunction]
pub(super) fn qc_molecule_to_chematic(
    molecule: &Bound<PyDict>,
) -> PyResult<(Mol, Vec<Vec<f64>>, f64, i64)> {
    let text = json_dumps_stripped(molecule)?;
    let qc = chematic_mol::parse_qcschema_molecule(&text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let view = chematic_mol::qc_molecule_to_chematic(&qc)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let charge = view.molecular_charge;
    let multiplicity = view.molecular_multiplicity;
    let (mol, coords) = qc_view_into_py(view);
    Ok((mol, coords, charge, multiplicity))
}

/// Parse a QCSchema `AtomicInput` JSON document (a molecule plus what to
/// compute on it: `driver`, `model`, `keywords`).
///
/// Returns a dict with every `AtomicInput` field (``"molecule"`` -- itself
/// shaped like :func:`parse_qcschema_molecule`'s return value minus the
/// ``"mol"``/``"coords"`` convenience keys --, ``"driver"``, ``"model"``,
/// ``"keywords"``, ``"protocols"``, ``"extras"``, ``"provenance"``, ...)
/// plus two top-level convenience keys built from ``"molecule"``:
///
/// - ``"mol"``: :class:`Mol`.
/// - ``"coords"``: ``[[x, y, z], ...]`` in Ångström.
///
/// :func:`write_atomic_input` strips both before re-serializing.
///
/// Raises:
///     ValueError: on malformed JSON or a schema violation.
///
/// Example::
///
///     result = chematic.parse_atomic_input(open("job.json").read())
///     print(result["driver"], result["model"]["method"])
#[pyfunction]
pub(super) fn parse_atomic_input<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let input =
        chematic_mol::parse_atomic_input(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let canonical = chematic_mol::write_atomic_input(&input);
    let d = json_loads(py, &canonical)?;
    let view = chematic_mol::qc_molecule_to_chematic(&input.molecule)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (mol, coords) = qc_view_into_py(view);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    Ok(d)
}

/// Serialize an `AtomicInput` dict (as returned by
/// :func:`parse_atomic_input`) to canonical QCSchema JSON text.
///
/// Raises:
///     ValueError: if `input` doesn't validate as a QCSchema `AtomicInput`
///     document.
#[pyfunction]
pub(super) fn write_atomic_input(input: &Bound<PyDict>) -> PyResult<String> {
    let text = json_dumps_stripped(input)?;
    let a = chematic_mol::parse_atomic_input(&text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_atomic_input(&a))
}

/// Parse a QCSchema `AtomicResult` JSON document (an `AtomicInput` plus the
/// computed result: `return_result`/`properties`/`success`/`error`/
/// `provenance`).
///
/// Same shape as :func:`parse_atomic_input`'s return value, plus
/// ``"return_result"`` (``float`` for an ``"energy"`` driver, ``list[float]``
/// for ``"gradient"``/``"hessian"``, or a ``dict`` for ``"properties"`` --
/// shaped by ``"driver"``, only present when ``"success"`` is ``True``),
/// ``"success"`` (bool), ``"error"`` (``{"error_type": str,
/// "error_message": str, ...}`` dict, only present when ``"success"`` is
/// ``False``), ``"properties"``, ``"provenance"``, ``"stdout"``/``"stderr"``.
///
/// Raises:
///     ValueError: on malformed JSON or a schema violation (including the
///     spec's own ``success``/``return_result``/``error`` consistency
///     rule).
///
/// Example::
///
///     result = chematic.parse_atomic_result(open("result.json").read())
///     if result["success"]:
///         print(result["return_result"])
#[pyfunction]
pub(super) fn parse_atomic_result<'py>(
    py: Python<'py>,
    text: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let result = chematic_mol::parse_atomic_result(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let canonical = chematic_mol::write_atomic_result(&result);
    let d = json_loads(py, &canonical)?;
    let view = chematic_mol::qc_molecule_to_chematic(&result.molecule)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (mol, coords) = qc_view_into_py(view);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    Ok(d)
}

/// Serialize an `AtomicResult` dict (as returned by
/// :func:`parse_atomic_result`) to canonical QCSchema JSON text.
///
/// Raises:
///     ValueError: if `result` doesn't validate as a QCSchema `AtomicResult`
///     document.
#[pyfunction]
pub(super) fn write_atomic_result(result: &Bound<PyDict>) -> PyResult<String> {
    let text = json_dumps_stripped(result)?;
    let r = chematic_mol::parse_atomic_result(&text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_atomic_result(&r))
}
