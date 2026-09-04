//! Molecule format parsers/writers exposed as free `#[pyfunction]`s (SMILES/SDF/MOL/CML/CJSON/MolJSON/CDXML/MOL2/PDBQT/GJF/CIF/InChI/PDB/XYZ/RXN) plus small serialization helpers.

use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

type SdfWithCoords = Vec<(Mol, String, Vec<Vec<f64>>)>;

/// Stable, typed reason string for a [`chematic_perception::StereoDiagnostic`]
/// -- never a free-form message (see `docs/rfcs/stereo2d_reader_integration_rfc.md`).
pub(crate) fn stereo_reason_str(
    reason: chematic_perception::StereoRejectionReason,
) -> &'static str {
    use chematic_perception::StereoRejectionReason::*;
    match reason {
        ContradictoryWedges => "contradictory_wedges",
        MissingCoordinate => "missing_coordinate",
        DegenerateGeometry => "degenerate_geometry",
        UnsupportedCoordination => "unsupported_coordination",
    }
}

pub(crate) fn stereo_diagnostics_to_py<'py>(
    py: Python<'py>,
    diagnostics: &[chematic_perception::StereoDiagnostic],
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    diagnostics
        .iter()
        .map(|d| {
            let dict = PyDict::new(py);
            dict.set_item("atom_idx", d.atom.0)?;
            dict.set_item("reason", stereo_reason_str(d.reason))?;
            Ok(dict)
        })
        .collect()
}

pub(crate) fn bitvec2048_to_bytes(fp: &chematic_fp::bitvec::BitVec2048) -> Vec<u8> {
    (0..256usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                if fp.get(byte_idx * 8 + bit) {
                    byte |= 1 << bit;
                }
            }
            byte
        })
        .collect()
}

pub(crate) fn flat_to_coords3d(coords: &[[f64; 3]]) -> chematic_3d::Coords3D {
    chematic_3d::Coords3D {
        points: coords
            .iter()
            .map(|c| chematic_3d::Point3::new(c[0], c[1], c[2]))
            .collect(),
    }
}

// ---------------------------------------------------------------------------
// Helper: convert MoleculeReport to Python dict
// ---------------------------------------------------------------------------

/// Parse a SMILES string and return a Mol.
///
/// Raises ``ValueError`` on invalid SMILES.
#[pyfunction]
fn from_smiles(smiles: &str) -> PyResult<Mol> {
    chematic_smiles::parse(smiles)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a CXSMILES string and return the molecule with CX metadata.
///
/// Returns a 2-tuple ``(mol, cx)`` where ``cx`` is a dict with:
///
/// - ``atom_labels``: list of atom label strings (or ``None`` per atom)
/// - ``atom_props``: list of ``{"atom_idx", "key", "value"}`` dicts
/// - ``atom_radicals``: list of radical class integers (or ``None`` per atom)
///
/// Raises ``ValueError`` on parse failure. CXSMILES without a CX extension
/// block behaves like :func:`from_smiles` (all CX fields are empty).
///
///     mol, cx = chematic.from_cxsmiles("CC |$R1;R2$|")
///     print(cx['atom_labels'])   # ['R1', 'R2']
#[pyfunction]
fn from_cxsmiles<'py>(py: Python<'py>, s: &str) -> PyResult<(Mol, Bound<'py, PyDict>)> {
    let cx =
        chematic_smiles::parse_cxsmiles(s).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let d = PyDict::new(py);
    let labels: Vec<Option<&str>> = cx.atom_labels.iter().map(|l| l.as_deref()).collect();
    d.set_item("atom_labels", labels)?;
    let props: Vec<Bound<'py, PyDict>> = cx
        .atom_props
        .iter()
        .map(|p| {
            let pd = PyDict::new(py);
            pd.set_item("atom_idx", p.atom.0 as usize).unwrap();
            pd.set_item("key", &p.key).unwrap();
            pd.set_item("value", &p.value).unwrap();
            pd
        })
        .collect();
    d.set_item("atom_props", props)?;
    let radicals: Vec<Option<u8>> = cx.atom_radicals.clone();
    d.set_item("atom_radicals", radicals)?;

    Ok((
        Mol {
            inner: Arc::new(cx.mol),
            props: Default::default(),
        },
        d,
    ))
}

/// Parse a MOL/SDF block and return a Mol.
///
/// Parse a condensed molecular formula (e.g., ``"CH3OH"``, ``"C6H12O6"``) into a :class:`Mol`.
///
/// Returns ``None`` if the formula is unknown or ambiguous.
/// Unlike SMILES parsing, condensed formulas may not encode connectivity uniquely;
/// this function uses a built-in formula→structure dictionary.
///
/// Equivalent to chempy's condensed formula support.
///
///     mol = chematic.from_condensed("CH3OH")  # methanol
///     if mol:
///         print(mol.smiles)  # CO
///
///     mol = chematic.from_condensed("C6H12O6")  # glucose
#[pyfunction]
fn from_condensed(formula: &str) -> Option<Mol> {
    chematic_chem::parse_condensed(formula).ok().map(|mol| Mol {
        inner: Arc::new(mol),
        props: Default::default(),
    })
}

/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol_block(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V2000 block and return the molecule with its 2D layout coordinates.
///
/// Returns a 3-tuple ``(mol, name, coords_2d)`` where:
///
/// - ``mol``: :class:`Mol` object
/// - ``name``: molecule name from the MOL header (may be empty)
/// - ``coords_2d``: list of ``[x, y]`` pairs (one per heavy atom, Å)
///
/// Raises ``ValueError`` on parse failure.
///
/// Use :func:`from_mol_block` if you only need the molecule graph.
/// Use this function when you want to preserve the 2D layout for display or
/// round-trip back to MOL format via :meth:`Mol.to_mol_block_2d`.
///
///     mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
///     new_block = mol.to_mol_block_2d(coords_2d, name=name)
#[pyfunction]
fn from_mol_block_with_coords(block: &str) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_with_coords(block)
        .map(|(mol, meta, coords)| {
            let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                meta.name,
                py_coords,
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V2000 block, returning stereo-perception diagnostics
/// alongside the molecule.
///
/// Returns a 4-tuple ``(mol, name, coords_2d, stereo_diagnostics)``.
/// ``stereo_diagnostics`` is a list of ``{"atom_idx": int, "reason": str}``
/// dicts, one per wedge/hash center that could not be resolved -- reason is
/// one of ``"contradictory_wedges"``, ``"missing_coordinate"``,
/// ``"degenerate_geometry"``, or ``"unsupported_coordination"``. Empty
/// unless a wedge/hash bond was present at some center and got rejected; an
/// atom with no wedge/hash bond at all never produces an entry.
///
/// Local tetrahedral parity (``Atom.chirality``) is always perceived
/// automatically -- this function differs from
/// :func:`from_mol_block_with_coords` only in also surfacing *why* any
/// center was rejected.
///
/// Raises ``ValueError`` on parse failure.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn from_mol_block_with_diagnostics<'py>(
    py: Python<'py>,
    block: &str,
) -> PyResult<(Mol, String, Vec<Vec<f64>>, Vec<Bound<'py, PyDict>>)> {
    let report = chematic_mol::read_mol_with_diagnostics(block)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let py_coords: Vec<Vec<f64>> = report.coords.iter().map(|(x, y)| vec![*x, *y]).collect();
    let diagnostics = stereo_diagnostics_to_py(py, &report.stereo_diagnostics)?;
    Ok((
        Mol {
            inner: Arc::new(report.mol),
            props: Default::default(),
        },
        report.metadata.name,
        py_coords,
        diagnostics,
    ))
}

/// Parse a multi-record SDF string and return all molecules with their 2D layout coordinates.
///
/// Returns a list of 3-tuples ``(mol, name, coords_2d)`` — one per SDF record.
/// Invalid records are silently skipped (same behaviour as :func:`iter_sdf`).
/// Resource-limit failures are raised as ``ValueError`` rather than skipped.
///
/// This is the batch equivalent of :func:`from_mol_block_with_coords`.
///
/// Note: this function has no diagnostics-returning variant. It inherits 2D
/// wedge/hash stereo *perception* automatically (same core as every other
/// MOL/SDF reader), but rejected-stereo diagnostics aren't surfaced here —
/// use :func:`iter_sdf`/:func:`iter_sdf_batched`/:func:`iter_sdf_str` and
/// ``SdfRecord.stereo_diagnostics()`` if you need them.
///
///     with open("library.sdf") as f:
///         records = chematic.parse_sdf_with_coords(f.read())
///     for mol, name, coords_2d in records:
///         new_block = mol.to_mol_block_2d(coords_2d, name=name)
#[pyfunction]
fn parse_sdf_with_coords(text: &str) -> PyResult<SdfWithCoords> {
    // Delegates to SdfRecordReader (line-anchored $$$$ scanning, and no
    // longer eats a legitimately blank MOL name line -- issue #171) instead
    // of a hand-rolled splitter, so this stays in sync with the one fixed
    // core implementation rather than drifting from it.
    let mut records = Vec::new();
    for result in chematic_mol::SdfRecordReader::new(text) {
        let rec = match result {
            Ok(rec) => rec,
            Err(chematic_mol::MolParseError::ResourceLimit {
                resource,
                actual,
                limit,
            }) => {
                return Err(PyValueError::new_err(format!(
                    "SDF {resource} exceeds limit {limit} (got {actual})"
                )));
            }
            Err(_) => continue,
        };
        records.push({
            let py_coords: Vec<Vec<f64>> = rec.coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(rec.mol),
                    props: Default::default(),
                },
                rec.meta.name,
                py_coords,
            )
        });
    }
    Ok(records)
}

/// Parse an MRV block and return the molecule with its 2D/3D layout coordinates.
///
/// Returns a 3-tuple ``(mol, coords_2d, coords_3d)`` -- either list is
/// empty if the source file didn't carry that dimensionality.
///
/// Use :func:`from_mrv_block` if you only need the molecule graph.
/// Use this function to preserve the layout for round-trip via
/// :meth:`Mol.to_mrv_block_with_coords`.
///
///     mol, coords_2d, coords_3d = chematic.from_mrv_block_with_coords(block)
///     new_block = mol.to_mrv_block_with_coords(coords_2d, coords_3d)
#[pyfunction]
#[allow(clippy::type_complexity)]
fn from_mrv_block_with_coords(mrv_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    chematic_mol::parse_mrv(mrv_str)
        .map(|rec| {
            let coords_2d: Vec<Vec<f64>> = rec
                .coordinates_2d
                .unwrap_or_default()
                .iter()
                .map(|c| vec![c[0], c[1]])
                .collect();
            let coords_3d: Vec<Vec<f64>> = rec
                .coordinates_3d
                .unwrap_or_default()
                .iter()
                .map(|c| vec![c[0], c[1], c[2]])
                .collect();
            (
                Mol {
                    inner: Arc::new(rec.mol),
                    props: Default::default(),
                },
                coords_2d,
                coords_3d,
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Chemical Markup Language (CML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cml") as f:
///         mol = chematic.from_cml(f.read())
#[pyfunction]
fn from_cml(cml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cml(cml_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a ChemAxon Marvin (.mrv) string into a ``Mol`` object.
///
/// S-groups, polymers, reactions, multicenter bonds, query atoms/bonds,
/// R-groups, enhanced stereo groups, and embedded/compressed data are
/// deliberately unsupported and raise ``ValueError`` (see
/// ``chematic_mol::mrv`` for the full scope boundary), same as any other
/// parse failure.
///
///     with open("molecule.mrv") as f:
///         mol = chematic.from_mrv_block(f.read())
#[pyfunction]
fn from_mrv_block(mrv_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_mrv(mrv_str)
        .map(|rec| Mol {
            inner: Arc::new(rec.mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a ChemicalJSON (.cjson) string.
///
/// Returns ``(mol, coords)`` where ``coords`` is a list of ``[x, y, z]``
/// coordinate triples (Å), one per heavy atom.  ``coords`` is empty when
/// the file has no ``atoms.coords.3d`` field.
///
/// ChemicalJSON is the native format of Avogadro 2 and the MolSSI
/// Open Chemistry toolkit.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, coords = chematic.from_cjson(open("mol.cjson").read())
///     print(mol.smiles)
///     # Round-trip:
///     open("out.cjson", "w").write(mol.to_cjson(coords))
#[pyfunction]
fn from_cjson(cjson_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, coords) =
        chematic_mol::parse_cjson(cjson_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let py_coords = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
        py_coords,
    ))
}

/// Parse a MolJSON string into a ``Mol`` object.
///
/// MolJSON is a JSON-based molecular representation designed for LLM
/// (large language model) compatibility.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol = chematic.from_moljson(open("mol.json").read())
///     # Round-trip:
///     json_str = mol.to_moljson()
///     mol2 = chematic.from_moljson(json_str)
#[pyfunction]
fn from_moljson(json_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_moljson(json_str)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a ChemDraw XML (CDXML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cdxml") as f:
///         mol = chematic.from_cdxml(f.read())
#[pyfunction]
fn from_cdxml(cdxml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cdxml(cdxml_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse CDXML page/presentation structure without dropping unknown objects.
#[pyfunction]
fn parse_cdxml_document_json(cdxml_str: &str) -> PyResult<String> {
    let document = chematic_mol::CdxmlDocument::parse(cdxml_str)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    serde_json::to_string(&document.to_json()).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Apply a loss-preserving page/presentation edit to a CDXML document.
///
/// ``edit_json`` is a ``CdxmlEdit`` command object, for example
/// ``{"kind":"set_page_attribute", "page_id":"p1", "key":"title", "value":"Page 1"}``.
/// Unknown presentation XML remains intact.
#[pyfunction]
fn edit_cdxml_document_json(cdxml_str: &str, edit_json: &str) -> PyResult<String> {
    let document = chematic_mol::CdxmlDocument::parse(cdxml_str)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let edited = document
        .apply_json_edit(edit_json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(edited.write())
}

/// Validate and normalize a typed Markush/polymer semantic model JSON.
#[pyfunction]
fn semantic_model_json(model_json: &str) -> PyResult<String> {
    let value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| PyValueError::new_err(format!("invalid semantic JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&value)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    serde_json::to_string(&model.to_json()).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Apply an explicit Markush selection command to a semantic model JSON.
#[pyfunction]
fn semantic_apply_json_command(model_json: &str, command_json: &str) -> PyResult<String> {
    let model_value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| PyValueError::new_err(format!("invalid semantic JSON: {e}")))?;
    let command: serde_json::Value = serde_json::from_str(command_json)
        .map_err(|e| PyValueError::new_err(format!("invalid semantic command JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&model_value)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let next = model
        .apply_json_command(&command)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    serde_json::to_string(&next.to_json()).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Expand a validated semantic model against a base SMILES.
#[pyfunction]
fn semantic_expand_json(base_smiles: &str, model_json: &str) -> PyResult<String> {
    let value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| PyValueError::new_err(format!("invalid semantic JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&value)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let base =
        chematic_smiles::parse(base_smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let expanded = model
        .expand(&base)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    serde_json::to_string(&expanded.to_json()).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse an MDL MOL V3000 (``V3000``) block into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("ligand_v3000.mol") as f:
///         mol = chematic.from_mol_v3000(f.read())
#[pyfunction]
fn from_mol_v3000(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol_v3000(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V3000 block and return the molecule with its 2D layout coordinates.
///
/// Returns a 3-tuple ``(mol, name, coords_2d)`` identical to
/// :func:`from_mol_block_with_coords` but for V3000 input.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, name, coords_2d = chematic.from_mol_v3000_with_coords(block)
///     new_block = mol.to_mol_v3000(coords_2d, name=name)
#[pyfunction]
fn from_mol_v3000_with_coords(block: &str) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_v3000_with_coords(block)
        .map(|(mol, meta, coords)| {
            let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                meta.name,
                py_coords,
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V3000 block, returning stereo-perception diagnostics
/// alongside the molecule.
///
/// Same shape as :func:`from_mol_block_with_diagnostics` but for V3000
/// input -- see that function for the diagnostics dict shape.
///
/// Raises ``ValueError`` on parse failure.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn from_mol_v3000_with_diagnostics<'py>(
    py: Python<'py>,
    block: &str,
) -> PyResult<(Mol, String, Vec<Vec<f64>>, Vec<Bound<'py, PyDict>>)> {
    let report = chematic_mol::read_mol_v3000_with_diagnostics(block)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let py_coords: Vec<Vec<f64>> = report.coords.iter().map(|(x, y)| vec![*x, *y]).collect();
    let diagnostics = stereo_diagnostics_to_py(py, &report.stereo_diagnostics)?;
    Ok((
        Mol {
            inner: Arc::new(report.mol),
            props: Default::default(),
        },
        report.metadata.name,
        py_coords,
        diagnostics,
    ))
}

/// Parse a Tripos MOL2 string into a ``Mol`` object.
///
/// Example::
///
///     with open("ligand.mol2") as f:
///         mol = chematic.from_mol2(f.read())
///     print(mol.mw)
///
/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol2(mol2_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol2(mol2_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse an AutoDock PDBQT string and return a :class:`Mol`.
///
/// Only the molecular graph (elements and bonds) is extracted; 3D coordinates
/// and partial charges are discarded.  To retain them, use the lower-level
/// :func:`chematic_mol.parse_pdbqt` Rust API directly.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     with open("ligand.pdbqt") as f:
///         mol = chematic.from_pdbqt(f.read())
///     print(mol.mw)
#[pyfunction]
fn from_pdbqt(pdbqt_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_pdbqt(pdbqt_str)
        .map(|(mol, _coords, _charges)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Gaussian input file (`.gjf` / `.com`) and return a :class:`Mol`.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     mol = chematic.from_gjf(open("mol.gjf").read())
#[pyfunction]
fn from_gjf(gjf_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_gjf(gjf_str)
        .map(|(mol, _coords, _charge, _mult)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Gaussian output file (`.log` / `.out`) and return a dict with
/// ``mol``, ``coords`` and ``scf_energy`` fields.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "scf_energy": float | None}``
///
/// Raises:
///     ValueError: when no `Standard orientation:` block is found.
#[pyfunction]
fn parse_gaussian_log<'py>(
    py: Python<'py>,
    log_str: &str,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let result = chematic_mol::parse_gaussian_log(log_str)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
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
    d.set_item("scf_energy", result.scf_energy)?;
    Ok(d)
}

/// Generate a Gaussian input file (`.gjf`) string from a molecule.
///
/// Args:
///     mol: The molecule to write.
///     coords: Atomic coordinates as ``[[x, y, z], ...]`` in Ångströms.
///     charge: Formal charge (default 0).
///     multiplicity: Spin multiplicity (default 1).
///     method: Route section keywords (default ``"B3LYP/6-31G* opt"``).
///     title: Job title comment (default ``"chematic"``).
///
/// Returns:
///     str: GJF file contents.
#[pyfunction]
#[pyo3(signature = (mol, coords, charge=0, multiplicity=1, method="B3LYP/6-31G* opt", title="chematic"))]
fn write_gjf(
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    charge: i32,
    multiplicity: u32,
    method: &str,
    title: &str,
) -> String {
    let c: Vec<(f64, f64, f64)> = coords.into_iter().map(|[x, y, z]| (x, y, z)).collect();
    chematic_mol::write_gjf(&mol.inner, &c, charge, multiplicity, method, title)
}

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
fn parse_cif<'py>(py: Python<'py>, cif_str: &str) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
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

/// Return True if the SMILES can be parsed without error.
#[pyfunction]
fn is_valid_smiles(smiles: &str) -> bool {
    chematic_smiles::parse(smiles).is_ok()
}

/// Return ``True`` if ``smarts`` is a valid SMARTS pattern, ``False`` otherwise.
///
/// Mirrors :func:`is_valid_smiles` for SMARTS pattern validation.
/// Useful for validating user-supplied SMARTS before calling
/// :func:`smarts_match` or :func:`smarts_find`.
///
///     chematic.is_valid_smarts("c1ccccc1")  # True
///     chematic.is_valid_smarts("[invalid")  # False
///     chematic.is_valid_smarts("[#6]-[#7]") # True
#[pyfunction]
fn is_valid_smarts(smarts: &str) -> bool {
    chematic_smarts::parse_smarts(smarts).is_ok()
}

/// Parse an InChI string and return a Mol.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol = chematic.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
#[pyfunction]
fn from_inchi(inchi: &str) -> PyResult<Mol> {
    chematic_inchi::parse_inchi(inchi)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a PDB string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances; only ATOM/HETATM records
/// are used.  Returns ``ValueError`` when no atoms are found.
///
///     mol, coords = chematic.from_pdb(open("ligand.pdb").read())
#[pyfunction]
fn from_pdb(pdb_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let atoms = chematic_3d::parse_pdb_atoms(pdb_str);
    if atoms.is_empty() {
        return Err(PyValueError::new_err(
            "no ATOM/HETATM records found in PDB input",
        ));
    }
    let (mol, c3d) = chematic_3d::pdb_to_molecule(&atoms);
    let coords = c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
        coords,
    ))
}

/// Parse an XYZ string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances.
///
///     mol, coords = chematic.from_xyz(open("molecule.xyz").read())
#[pyfunction]
fn from_xyz(xyz_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
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
fn extxyz_frame_to_pydict<'py>(
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

fn extxyz_property_rows<T>(
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
fn from_extxyz<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
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
fn from_extxyz_all<'py>(py: Python<'py>, text: &str) -> PyResult<Vec<Bound<'py, PyDict>>> {
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
fn to_extxyz<'py>(
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
fn parse_formula<'py>(formula: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
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

/// Compute atom economy of a reaction (green chemistry metric).
/// E-factor (Environmental Factor) — waste-to-product mass ratio.
///
/// E-factor = waste_mass / product_mass.  Lower is greener.
/// Fine chemicals typically E=5–50; pharmaceuticals E=25–100.
///
///     ef = chematic.e_factor(waste_kg=90.0, product_kg=10.0)  # → 9.0
/// Fast structural hash for deduplication — one int per molecule.
///
/// Molecules with the same canonical graph return identical hashes.
/// This is a fast screening hash only. Hash collisions and canonical residuals
/// are possible; do not treat it as a definitive identity key.
///
/// Equivalent to RDKit's ``rdMolHash.MolHash()``.
///
///     seen = set()
///     unique = [m for m in mols if (h := chematic.mol_hash(m)) not in seen and not seen.add(h)]
#[pyfunction]
fn mol_hash(mol: &Mol) -> u64 {
    chematic_chem::mol_hash(&mol.inner)
}

/// Check whether two molecules have the same current canonical representation.
///
/// This is not a fail-closed graph-isomorphism proof. Use the stable-key API
/// when an indeterminate result must be distinguishable from ``False``.
///
///     assert chematic.are_identical(
///         chematic.from_smiles("c1ccccc1"),
///         chematic.from_smiles("C1=CC=CC=C1"),  # kekulé form
///     )
#[pyfunction]
fn are_identical(mol1: &Mol, mol2: &Mol) -> bool {
    chematic_chem::are_identical(&mol1.inner, &mol2.inner)
}

/// Compare molecules using chematic's fail-closed canonical identity key.
///
/// Returns ``True`` or ``False`` for converged keys and ``None`` when either
/// molecule belongs to a known canonical residual class.
#[pyfunction]
fn stable_are_identical(mol1: &Mol, mol2: &Mol) -> Option<bool> {
    chematic_chem::stable_are_identical(&mol1.inner, &mol2.inner)
}

/// Normalize and re-serialize a reaction SMILES.
///
/// Parses the reaction SMILES and writes it back in canonical form.
/// Useful for standardizing reaction data before storing or comparing.
///
///     canon = chematic.write_reaction("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn write_reaction(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

///
/// atom_economy = MW(desired products) / MW(all reactants) × 100.
/// Parse a MDL RXN V2000 file and return the canonical reaction SMILES.
///
/// Raises ``ValueError`` on parse failure.
///
///     rxn_smiles = chematic.from_rxn_file(text)
///     ae = chematic.atom_economy(rxn_smiles)
#[pyfunction]
fn from_rxn_file(text: &str) -> PyResult<String> {
    let rxn =
        chematic_mol::parse_rxn_file(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

/// Convert a reaction SMILES string to MDL RXN V2000 format.
///
/// Raises ``ValueError`` on invalid reaction SMILES.
///
///     block = chematic.to_rxn_file("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn to_rxn_file(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_rxn_file(&rxn))
}

/// Parse an MDL RXN V2000 file into the typed, loss-aware reaction-document
/// JSON contract shared with the WASM/Node bindings.
#[pyfunction]
fn from_rxn_document_json(text: &str) -> PyResult<String> {
    let document =
        chematic_mol::parse_rxn_document(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    serde_json::to_string(&document).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Write typed reaction-document JSON as MDL RXN V2000.
///
/// Rich fields with no RXN V2000 representation raise ``ValueError`` instead
/// of being silently discarded.
#[pyfunction]
fn to_rxn_document_json(document_json: &str) -> PyResult<String> {
    let document: chematic_rxn::ReactionDocument = serde_json::from_str(document_json)
        .map_err(|e| PyValueError::new_err(format!("invalid reaction document JSON: {e}")))?;
    chematic_mol::write_rxn_document(&document).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Find top-K nearest neighbors from precomputed fingerprint byte arrays.
///
/// More efficient than :func:`top_k_similar_fp` when the same ``db_fps`` list is
/// reused across multiple queries (fingerprints computed only once).
///
///     db_fps = [mol.ecfp4() for mol in library]   # compute once
///     for query in queries:
///         hits = chematic.nearest_neighbors_from_fp(query.ecfp4(), db_fps, k=10)
///         for idx, score in hits:
///             print(library_smiles[idx], score)
#[pyfunction]
#[pyo3(signature = (query_fp, db_fps, k = 10))]
fn nearest_neighbors_from_fp(query_fp: &[u8], db_fps: Vec<Vec<u8>>, k: usize) -> Vec<(usize, f64)> {
    let qa: u32 = query_fp.iter().map(|b| b.count_ones()).sum();
    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .enumerate()
        .filter_map(|(i, fp)| {
            if fp.len() != query_fp.len() {
                return None;
            }
            let and: u32 = query_fp
                .iter()
                .zip(fp.iter())
                .map(|(a, b)| (a & b).count_ones())
                .sum();
            let db_cnt: u32 = fp.iter().map(|b| b.count_ones()).sum();
            let or = qa + db_cnt - and;
            if or == 0 {
                return None;
            }
            Some((i, and as f64 / or as f64))
        })
        .collect();
    scores.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    scores.truncate(k);
    scores
}

/// Parse a ``.smi`` file (tab/space-separated SMILES + name) into (Mol, name) pairs.
///
/// Each line is ``SMILES[<tab>name]``. Lines with invalid SMILES are silently skipped.
/// Resource-limit failures are raised as ``ValueError`` rather than skipped.
/// Comment lines starting with ``#`` and blank lines are ignored.
/// Equivalent to RDKit's ``Chem.SmilesMolSupplier``.
///
///     records = chematic.parse_smi_file(open("library.smi").read())
///     for mol, name in records:
///         print(name, mol.mw)
#[pyfunction]
fn parse_smi_file(content: &str) -> PyResult<Vec<(Mol, String)>> {
    let mut records = Vec::new();
    for result in chematic_smiles::parse_smi_file(content) {
        let (mol, name) = match result {
            Ok(record) => record,
            Err(chematic_smiles::SmilesError::ResourceLimit {
                resource,
                actual,
                limit,
            }) => {
                return Err(PyValueError::new_err(format!(
                    "SMILES {resource} exceeds limit {limit} (got {actual})"
                )));
            }
            Err(_) => continue,
        };
        records.push({
            (
                Mol {
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                name,
            )
        });
    }
    Ok(records)
}

/// Write (Mol, name) pairs to ``.smi`` format.
///
/// Output format: ``SMILES<TAB>name<NEWLINE>`` per record (name omitted if empty).
/// Equivalent to RDKit's ``Chem.SmilesWriter``.
///
///     text = chematic.write_smi_file([(mol1, "cpd1"), (mol2, "cpd2")])
///     with open("output.smi", "w") as f:
///         f.write(text)
#[pyfunction]
fn write_smi_file(records: Vec<(Mol, String)>) -> PyResult<String> {
    const MAX_RECORDS: usize = 100_000;
    const MAX_OUTPUT_BYTES: usize = 16 * 1024 * 1024;
    if records.len() > MAX_RECORDS {
        return Err(PyValueError::new_err(format!(
            "SMILES output exceeds maximum record count ({})",
            MAX_RECORDS
        )));
    }
    let mut out = String::new();
    for (mol, name) in &records {
        let smiles = chematic_smiles::canonical_smiles(&mol.inner);
        let line_bytes = smiles
            .len()
            .checked_add(name.len())
            .and_then(|n| n.checked_add(if name.is_empty() { 1 } else { 2 }))
            .ok_or_else(|| PyValueError::new_err("SMILES output size overflow"))?;
        if out
            .len()
            .checked_add(line_bytes)
            .is_none_or(|size| size > MAX_OUTPUT_BYTES)
        {
            return Err(PyValueError::new_err(format!(
                "SMILES output exceeds maximum size ({} bytes)",
                MAX_OUTPUT_BYTES
            )));
        }
        if name.is_empty() {
            out.push_str(&smiles);
        } else {
            out.push_str(&smiles);
            out.push('\t');
            out.push_str(name);
        }
        out.push('\n');
    }
    Ok(out)
}

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

fn element_from_symbol(sym: &str) -> PyResult<chematic_core::Element> {
    chematic_core::Element::from_symbol(sym)
        .ok_or_else(|| PyValueError::new_err(format!("unknown element symbol: {sym:?}")))
}

// ---------------------------------------------------------------------------
// mmCIF
// ---------------------------------------------------------------------------

fn mmcif_atom_to_pydict<'py>(
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

fn pydict_to_mmcif_atom(d: &Bound<PyDict>) -> PyResult<chematic_mol::MmcifAtomRecord> {
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

fn unit_cell_to_pydict<'py>(
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

fn pydict_to_unit_cell(d: &Bound<PyDict>) -> PyResult<chematic_mol::UnitCell> {
    Ok(chematic_mol::UnitCell {
        a: dict_get_required(d, "a")?,
        b: dict_get_required(d, "b")?,
        c: dict_get_required(d, "c")?,
        alpha: dict_get_required(d, "alpha")?,
        beta: dict_get_required(d, "beta")?,
        gamma: dict_get_required(d, "gamma")?,
    })
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
fn parse_mmcif<'py>(
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
    let result = chematic_mol::parse_mmcif_with_limits(text, &limits)
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
fn write_mmcif(
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

fn pqr_atom_to_pydict<'py>(
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
fn pydict_to_pqr_atom(d: &Bound<PyDict>) -> PyResult<chematic_mol::PqrAtomRecord> {
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
fn parse_pqr<'py>(
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
fn write_pqr(atoms: Vec<Bound<PyDict>>) -> PyResult<String> {
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
fn infer_element(group_pdb: &str, res_name: &str, atom_name: &str) -> Option<String> {
    chematic_mol::infer_element(group_pdb, res_name, atom_name).map(|e| e.symbol().to_string())
}

// ---------------------------------------------------------------------------
// ORCA input
// ---------------------------------------------------------------------------

fn orca_block_to_pydict<'py>(
    py: Python<'py>,
    b: &chematic_mol::OrcaBlock,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("name", &b.name)?;
    d.set_item("raw", &b.raw)?;
    d.set_item("has_end", b.has_end)?;
    Ok(d)
}

fn pydict_to_orca_block(d: &Bound<PyDict>) -> PyResult<chematic_mol::OrcaBlock> {
    Ok(chematic_mol::OrcaBlock {
        name: dict_get_required(d, "name")?,
        raw: dict_get(d, "raw", String::new())?,
        has_end: dict_get(d, "has_end", true)?,
    })
}

fn orca_atom_to_pydict<'py>(
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

fn pydict_to_orca_atom(d: &Bound<PyDict>) -> PyResult<chematic_mol::OrcaAtom> {
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

fn orca_coords_to_pydict<'py>(
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

fn pydict_to_orca_coords(d: &Bound<PyDict>) -> PyResult<chematic_mol::OrcaCoords> {
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
fn parse_orca_input<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
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
fn write_orca_input(input: &Bound<PyDict>) -> PyResult<String> {
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
fn parse_orca_output<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
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

fn json_loads<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let value = py.import("json")?.call_method1("loads", (text,))?;
    value
        .cast::<PyDict>()
        .map_err(|e| PyValueError::new_err(e.to_string()))
        .cloned()
}

/// `json.dumps` a dict, first stripping the `"mol"`/`"coords"` convenience
/// keys the `parse_*` functions in this section add (see section docs).
fn json_dumps_stripped(d: &Bound<PyDict>) -> PyResult<String> {
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
fn qc_view_into_py(view: chematic_mol::ChematicMoleculeView) -> (Mol, Vec<Vec<f64>>) {
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
fn parse_qcschema_molecule<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
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
fn write_qcschema_molecule(molecule: &Bound<PyDict>) -> PyResult<String> {
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
fn chematic_to_qc_molecule<'py>(
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
fn qc_molecule_to_chematic(molecule: &Bound<PyDict>) -> PyResult<(Mol, Vec<Vec<f64>>, f64, i64)> {
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
fn parse_atomic_input<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
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
fn write_atomic_input(input: &Bound<PyDict>) -> PyResult<String> {
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
fn parse_atomic_result<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
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
fn write_atomic_result(result: &Bound<PyDict>) -> PyResult<String> {
    let text = json_dumps_stripped(result)?;
    let r = chematic_mol::parse_atomic_result(&text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_atomic_result(&r))
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(from_cxsmiles, m)?)?;
    m.add_function(wrap_pyfunction!(from_condensed, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block_with_diagnostics, m)?)?;
    m.add_function(wrap_pyfunction!(parse_sdf_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_cml, m)?)?;
    m.add_function(wrap_pyfunction!(from_mrv_block, m)?)?;
    m.add_function(wrap_pyfunction!(from_mrv_block_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_cjson, m)?)?;
    m.add_function(wrap_pyfunction!(from_moljson, m)?)?;
    m.add_function(wrap_pyfunction!(from_cdxml, m)?)?;
    m.add_function(wrap_pyfunction!(parse_cdxml_document_json, m)?)?;
    m.add_function(wrap_pyfunction!(edit_cdxml_document_json, m)?)?;
    m.add_function(wrap_pyfunction!(semantic_model_json, m)?)?;
    m.add_function(wrap_pyfunction!(semantic_apply_json_command, m)?)?;
    m.add_function(wrap_pyfunction!(semantic_expand_json, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000_with_diagnostics, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol2, m)?)?;
    m.add_function(wrap_pyfunction!(from_pdbqt, m)?)?;
    m.add_function(wrap_pyfunction!(from_gjf, m)?)?;
    m.add_function(wrap_pyfunction!(parse_gaussian_log, m)?)?;
    m.add_function(wrap_pyfunction!(write_gjf, m)?)?;
    m.add_function(wrap_pyfunction!(parse_cif, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(from_inchi, m)?)?;
    m.add_function(wrap_pyfunction!(from_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(from_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(from_extxyz, m)?)?;
    m.add_function(wrap_pyfunction!(from_extxyz_all, m)?)?;
    m.add_function(wrap_pyfunction!(to_extxyz, m)?)?;
    m.add_function(wrap_pyfunction!(parse_formula, m)?)?;
    m.add_function(wrap_pyfunction!(mol_hash, m)?)?;
    m.add_function(wrap_pyfunction!(are_identical, m)?)?;
    m.add_function(wrap_pyfunction!(stable_are_identical, m)?)?;
    m.add_function(wrap_pyfunction!(write_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(from_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(to_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(from_rxn_document_json, m)?)?;
    m.add_function(wrap_pyfunction!(to_rxn_document_json, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_neighbors_from_fp, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smi_file, m)?)?;
    m.add_function(wrap_pyfunction!(write_smi_file, m)?)?;
    m.add_function(wrap_pyfunction!(parse_mmcif, m)?)?;
    m.add_function(wrap_pyfunction!(write_mmcif, m)?)?;
    m.add_function(wrap_pyfunction!(parse_pqr, m)?)?;
    m.add_function(wrap_pyfunction!(write_pqr, m)?)?;
    m.add_function(wrap_pyfunction!(infer_element, m)?)?;
    m.add_function(wrap_pyfunction!(parse_orca_input, m)?)?;
    m.add_function(wrap_pyfunction!(write_orca_input, m)?)?;
    m.add_function(wrap_pyfunction!(parse_orca_output, m)?)?;
    m.add_function(wrap_pyfunction!(parse_qcschema_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(write_qcschema_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(chematic_to_qc_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(qc_molecule_to_chematic, m)?)?;
    m.add_function(wrap_pyfunction!(parse_atomic_input, m)?)?;
    m.add_function(wrap_pyfunction!(write_atomic_input, m)?)?;
    m.add_function(wrap_pyfunction!(parse_atomic_result, m)?)?;
    m.add_function(wrap_pyfunction!(write_atomic_result, m)?)?;
    Ok(())
}
