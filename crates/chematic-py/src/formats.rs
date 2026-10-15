//! Molecule format parsers/writers exposed as free `#[pyfunction]`s (SMILES/SDF/MOL/CML/CJSON/MolJSON/CDXML/MOL2/PDBQT/GJF/CIF/InChI/PDB/XYZ/RXN) plus small serialization helpers.

use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

mod hin;
use hin::*;

mod ccd;
use ccd::*;

mod cjson;
use cjson::*;
mod cif;
use cif::*;
mod cdjson;
use cdjson::*;
mod fhiaims;
use fhiaims::*;
mod dmol3;
use dmol3::*;
mod mopac;
use mopac::*;
mod exyz;
use exyz::*;
mod gro;
mod gzmat;
use gzmat::*;
mod rdfile;
use rdfile::*;
mod turbomole;
use turbomole::*;
mod interchange;
mod materials;
mod mmcif_pdb;
mod pdb_selection;
use pdb_selection::*;
mod mmcif_residues;
use mmcif_pdb::*;
use mmcif_residues::*;
mod pqr_pdbqt;
use pqr_pdbqt::*;
mod xyz;
use gro::parse_gro;
use interchange::*;
use materials::*;
pub(crate) use materials::{dict_get, dict_get_opt, dict_get_required};
pub(crate) use xyz::extxyz_frame_to_pydict;
use xyz::*;

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
    fp.to_le_bytes().to_vec()
}

pub(crate) fn flat_to_coords3d(coords: &[[f64; 3]]) -> chematic_3d::Coords3D {
    chematic_3d::Coords3D {
        points: coords
            .iter()
            .map(|c| chematic_3d::Point3::new(c[0], c[1], c[2]))
            .collect(),
    }
}

/// Convert caller-provided coordinates after validating the public Python
/// boundary.  Most 3D kernels index coordinates by molecule atom index, so a
/// short vector would otherwise surface as a Rust panic instead of a Python
/// input error.
pub(crate) fn checked_coords3d(
    mol: &chematic_core::Molecule,
    coords: &[[f64; 3]],
) -> PyResult<chematic_3d::Coords3D> {
    if coords.len() != mol.atom_count() {
        return Err(PyValueError::new_err(format!(
            "expected {} coordinates, got {}",
            mol.atom_count(),
            coords.len()
        )));
    }
    if coords.iter().flatten().any(|value| !value.is_finite()) {
        return Err(PyValueError::new_err("coordinates must all be finite"));
    }
    Ok(flat_to_coords3d(coords))
}

// ---------------------------------------------------------------------------
// Helper: convert MoleculeReport to Python dict
// ---------------------------------------------------------------------------

/// Parse a SMILES string and return a Mol.
///
/// Raises :class:`ChematicInputError` (a ``ValueError``; ``category``
/// ``"malformed"``, ``code`` ``"smiles_parse"``) on invalid SMILES, and
/// with ``code`` ``"smiles_valence"`` when the graph gives a neutral oxygen
/// or fluorine more bonds than it permits (``CO(C)C``).
#[pyfunction]
fn from_smiles(smiles: &str) -> PyResult<Mol> {
    chematic_smiles::parse(smiles)
        .map(|mol| Mol {
            // Descriptor and ring preparation stays lazy.  Calling one of
            // those APIs memoizes its result on the molecule, so repeated
            // workflows still reuse the work without charging parse-only or
            // single-property callers for unrelated perception.
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| match e {
            chematic_smiles::SmilesError::InvalidValence { .. } => {
                crate::errors::input_error("malformed", "smiles_valence", "smiles", e.to_string())
            }
            e => crate::errors::malformed("smiles", e.to_string()),
        })
}

/// Canonicalize a list of SMILES without aborting on an invalid record.
///
/// Returns a complete JSON manifest containing one terminal record per input,
/// its original index, and explicit accepted/rejected/refused/skipped
/// accounting. The Rust batch API remains the source of truth for ordering
/// and parser limits.
#[pyfunction]
fn canonicalize_smiles_batch_json(smiles: Vec<String>) -> PyResult<String> {
    let canonicalizer = chematic_smiles::SmilesBatchCanonicalizer::default();
    let result = canonicalizer.canonicalize_with_result(smiles.iter());
    let all_succeeded = result.all_succeeded();
    let records = result
        .records
        .into_iter()
        .map(|record| match record.result {
            chematic_smiles::BatchCanonicalization::Accepted { canonical_smiles } => {
                serde_json::json!({
                    "input_index": record.input_index,
                    "input": record.input,
                    "status": "accepted",
                    "canonical_smiles": canonical_smiles,
                    "error": null,
                    "error_stage": null,
                })
            }
            chematic_smiles::BatchCanonicalization::Rejected { error } => serde_json::json!({
                "input_index": record.input_index,
                "input": record.input,
                "status": "rejected",
                "error": error,
                "error_stage": "parse",
            }),
        });
    let records = records.collect::<Vec<_>>();
    serde_json::to_string(&serde_json::json!({
        "schema_version": 1,
        "operation": "canonicalize_smiles",
        "status": "complete",
        "record_count": records.len(),
        "input_count": result.input_count,
        "accepted_count": result.accepted_count,
        "rejected_count": result.rejected_count,
        "refused_count": result.refused_count,
        "skipped_count": result.skipped_count,
        "all_succeeded": all_succeeded,
        "records": records,
    }))
    .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Stateful Python adapter for an unknown-length SMILES source.
///
/// `observe` records rows received from a producer, while `process_next_json`
/// gives each observed row a terminal parse result.  `finish_json` is the
/// normal EOF path.  A cancelled or failed producer calls `stop_json`, which
/// retains the observed-but-unprocessed prefix and reports the unread suffix
/// as unknown rather than silently calling it skipped or successful.
#[pyclass(name = "SmilesBatchStream")]
struct PySmilesBatchStream {
    stream: Option<chematic_smiles::SmilesBatchStream>,
}

fn stream_record_json(record: &chematic_smiles::BatchCanonicalRecord) -> serde_json::Value {
    match &record.result {
        chematic_smiles::BatchCanonicalization::Accepted { canonical_smiles } => {
            serde_json::json!({
                "input_index": record.input_index,
                "input": record.input,
                "status": "accepted",
                "canonical_smiles": canonical_smiles,
                "error": null,
                "error_stage": null,
            })
        }
        chematic_smiles::BatchCanonicalization::Rejected { error } => serde_json::json!({
            "input_index": record.input_index,
            "input": record.input,
            "status": "rejected",
            "error": error,
            "error_stage": "parse",
        }),
    }
}

fn stream_manifest_json(result: chematic_smiles::BatchCanonicalStreamResult) -> PyResult<String> {
    let complete = result.stream_complete;
    serde_json::to_string(&serde_json::json!({
        "schema_version": 1,
        "operation": "canonicalize_smiles_stream",
        "input_kind": "unknown_length_stream",
        "status": if complete { "complete" } else { "incomplete" },
        "observed_input_count": result.observed_input_count,
        "completed_count": result.records.len(),
        "accepted_count": result.accepted_count,
        "rejected_count": result.rejected_count,
        "refused_count": 0,
        "skipped_count": 0,
        "unprocessed_observed_count": result.unprocessed_observed_count,
        "unread_input": if complete { serde_json::json!(0) } else { serde_json::json!("unknown") },
        "all_succeeded": result.all_succeeded(),
        "terminal_reason": result.terminal_reason.map(|reason| reason.as_str()),
        "records": result.records.iter().map(stream_record_json).collect::<Vec<_>>(),
    }))
    .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[pymethods]
impl PySmilesBatchStream {
    #[new]
    fn new() -> Self {
        Self {
            stream: Some(chematic_smiles::SmilesBatchCanonicalizer::default().stream()),
        }
    }

    fn observe(&mut self, smiles: String) -> PyResult<usize> {
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| PyValueError::new_err("stream is already terminal"))?;
        Ok(stream.observe(smiles))
    }

    fn process_next_json(&mut self) -> PyResult<String> {
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| PyValueError::new_err("stream is already terminal"))?;
        serde_json::to_string(&stream.process_next().map(stream_record_json))
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn finish_json(&mut self) -> PyResult<String> {
        let stream = self
            .stream
            .take()
            .ok_or_else(|| PyValueError::new_err("stream is already terminal"))?;
        stream_manifest_json(stream.finish())
    }

    fn stop_json(&mut self, terminal_reason: String) -> PyResult<String> {
        let terminal_reason = terminal_reason
            .parse::<chematic_smiles::StreamTerminalReason>()
            .map_err(PyValueError::new_err)?;
        let stream = self
            .stream
            .take()
            .ok_or_else(|| PyValueError::new_err("stream is already terminal"))?;
        stream_manifest_json(stream.stop(terminal_reason))
    }
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

/// Raises :class:`ChematicInputError` (a ``ValueError``; ``code``
/// ``"mol_block_parse"``) on parse failure.
#[pyfunction]
fn from_mol_block(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| crate::errors::malformed("mol_block", e.to_string()))
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
        .map_err(|e| crate::errors::malformed("mol_block", e.to_string()))
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

/// Parse a structurally valid, non-empty CML string into a ``Mol`` object.
///
/// Unlike :func:`from_cml`, this opt-in entry point rejects missing molecules,
/// empty molecules, and mismatched or unclosed XML elements.
#[pyfunction]
fn from_cml_strict(cml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cml_strict(cml_str)
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
    semantic_expand_json_with_limits(
        base_smiles,
        model_json,
        chematic_mol::SemanticExpansionLimits::default().max_atoms,
        chematic_mol::SemanticExpansionLimits::default().max_repeat_count,
    )
}

/// Expand a validated semantic model with explicit finite atom/repeat budgets.
#[pyfunction]
fn semantic_expand_json_with_limits(
    base_smiles: &str,
    model_json: &str,
    max_atoms: usize,
    max_repeat_count: u32,
) -> PyResult<String> {
    let value: serde_json::Value = serde_json::from_str(model_json)
        .map_err(|e| PyValueError::new_err(format!("invalid semantic JSON: {e}")))?;
    let model = chematic_mol::SemanticModel::from_json(&value)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let base =
        chematic_smiles::parse(base_smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let expanded = model
        .expand_with_limits(
            &base,
            &chematic_mol::SemanticExpansionLimits {
                max_atoms,
                max_repeat_count,
            },
        )
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
/// With the IUPAC InChI library built in (the default wheel), the structure
/// is reconstructed like RDKit's ``Chem.MolFromInchi``; otherwise (or if the
/// library rejects the string) the pure-Rust parser is used.
///
/// Raises :class:`ChematicInputError` (a ``ValueError``; ``code``
/// ``"inchi_parse"``) on parse failure.
///
///     mol = chematic.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
#[pyfunction]
fn from_inchi(inchi: &str) -> PyResult<Mol> {
    #[cfg(feature = "native-inchi")]
    if let Ok(mol) = chematic_inchi::rdkit_mol_from_inchi(inchi) {
        return Ok(Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        });
    }
    chematic_inchi::parse_inchi(inchi)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| crate::errors::malformed("inchi", e.to_string()))
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

/// RDKit-compatible ``Chem.MolFromPDBBlock(text, sanitize, removeHs, flavor,
/// proximityBonding)`` (RDKit 2026.03.1): ``(mol, coords)`` in RDKit's atom
/// order, or ``None`` where RDKit returns no molecule (including parse and
/// sanitization failures). Bonds come from
/// CONECT records and RDKit's proximity bonding, chirality from the 3D
/// coordinates (tetrahedral, square-planar, trigonal-bipyramidal and
/// octahedral); ``mol.rdkit_smiles`` is ``Chem.MolToSmiles`` of RDKit's
/// molecule. Raises ``ValueError`` for features the port does not model
/// (zero-order bonds).
#[pyfunction]
#[pyo3(signature = (text, sanitize = true, remove_hs = true, proximity_bonding = true, flavor = 0))]
fn rdkit_from_pdb_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    proximity_bonding: bool,
    flavor: u32,
) -> PyResult<Option<(Mol, Vec<[f64; 3]>)>> {
    let res = match chematic_smiles::rdkit_mol_from_pdb_block(
        text,
        sanitize,
        remove_hs,
        flavor,
        proximity_bonding,
    ) {
        Ok(res) => res,
        // RDKit logs the failure and returns None.
        Err(chematic_smiles::RdkitSmilesError::Sanitization(_)) => None,
        Err(e) => return Err(PyValueError::new_err(e.to_string())),
    };
    Ok(res.map(|r| {
        (
            Mol {
                inner: Arc::new(r.molecule),
                props: Default::default(),
            },
            r.coords,
        )
    }))
}

/// RDKit-compatible ``Chem.MolFromXYZBlock(text)`` (RDKit 2026.03.1):
/// ``(mol, coords)`` with the atoms in file order and no bonds, or ``None``
/// where RDKit's parser fails.
#[pyfunction]
fn rdkit_from_xyz_block(text: &str) -> PyResult<Option<(Mol, Vec<[f64; 3]>)>> {
    match chematic_smiles::rdkit_mol_from_xyz_block(text) {
        Ok((m, coords)) => Ok(Some((
            Mol {
                inner: Arc::new(m),
                props: Default::default(),
            },
            coords,
        ))),
        Err(chematic_smiles::RdkitSmilesError::Sanitization(_)) => Ok(None),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

/// RDKit-compatible ``Chem.MolFromMol2Block(text, sanitize, removeHs,
/// cleanupSubstructures)`` (RDKit 2026.03.1): ``(mol, coords)`` in RDKit's
/// atom order, or ``None`` where RDKit returns no molecule (parse errors,
/// failed Corina substructure cleanup, sanitization failures). Formal
/// charges come from the Tripos atom types (``N.4``, ``O.co2``, ``C.cat``
/// and RDKit's valence-based guess) or ``UNITY_ATOM_ATTR`` records;
/// tetrahedral and double-bond stereo from the 3D coordinates. Raises
/// ``ValueError`` for features the port does not model (Tripos query atoms
/// ``Du``/``ANY``/``HEV``/``HET``/``HAL``, ``du``/``un`` bonds).
#[pyfunction]
#[pyo3(signature = (text, sanitize = true, remove_hs = true, cleanup_substructures = true))]
fn rdkit_from_mol2_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> PyResult<Option<(Mol, Vec<[f64; 3]>)>> {
    match chematic_smiles::rdkit_mol_from_mol2_block(
        text,
        sanitize,
        remove_hs,
        cleanup_substructures,
    ) {
        Ok(r) => Ok(Some((
            Mol {
                inner: Arc::new(r.molecule),
                props: Default::default(),
            },
            r.coords,
        ))),
        Err(chematic_smiles::RdkitSmilesError::Sanitization(_)) => Ok(None),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

/// ``Chem.MolToSmiles(Chem.MolFromMol2Block(text, ...))`` (RDKit
/// 2026.03.1), or ``None`` where RDKit returns no molecule.
#[pyfunction]
#[pyo3(signature = (text, remove_hs = true, cleanup_substructures = true))]
fn rdkit_mol2_block_to_smiles(
    text: &str,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> PyResult<Option<String>> {
    match chematic_smiles::rdkit_mol_from_mol2_block(text, true, remove_hs, cleanup_substructures) {
        Ok(r) => Ok(Some(r.smiles)),
        Err(chematic_smiles::RdkitSmilesError::Sanitization(_)) => Ok(None),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

/// ``Chem.MolToSmiles(Chem.MolFromPDBBlock(text, ...))`` (RDKit 2026.03.1),
/// or ``None`` where RDKit returns no molecule.
#[pyfunction]
#[pyo3(signature = (text, remove_hs = true, proximity_bonding = true, flavor = 0))]
fn rdkit_pdb_block_to_smiles(
    text: &str,
    remove_hs: bool,
    proximity_bonding: bool,
    flavor: u32,
) -> PyResult<Option<String>> {
    match chematic_smiles::rdkit_mol_from_pdb_block(
        text,
        true,
        remove_hs,
        flavor,
        proximity_bonding,
    ) {
        Ok(r) => Ok(r.map(|r| r.smiles)),
        Err(chematic_smiles::RdkitSmilesError::Sanitization(_)) => Ok(None),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

/// Parse a PDB string with fixed-column validation and return ``(Mol, coords)``.
///
/// Unlike :func:`from_pdb`, this opt-in API rejects ATOM/HETATM records with
/// missing, non-numeric, or non-finite serial, residue sequence, or XYZ fields.
/// Non-atom records remain ignored and the same finite resource limits as the
/// Rust core are applied.
#[pyfunction]
fn from_pdb_strict(pdb_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let atoms =
        chematic_3d::parse_pdb_atoms_strict(pdb_str, &chematic_3d::PdbParseLimits::default())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
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

/// Fast 64-bit structural hash for deduplication: FNV-1a of chematic's own
/// canonical SMILES (:attr:`Mol.smiles`).
///
/// Molecules with the same chematic canonical SMILES get the same integer.
/// It is a screening hash, not an identity key: collisions are possible,
/// and so are canonicalization residuals (two depictions of one molecule
/// with different canonical SMILES). It is not RDKit's
/// ``rdMolHash.MolHash`` (which returns strings such as scaffolds or
/// tautomer hashes); for that, use :meth:`Mol.rdkit_mol_hash`.
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
    let document = chematic_rxn::ReactionDocument::from_json_str(document_json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    chematic_mol::write_rxn_document(&document).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Find top-K nearest neighbors from precomputed fingerprint byte arrays.
///
/// More efficient than :func:`top_k_similar_fp` when the same ``db_fps`` list is
/// reused across multiple queries (fingerprints computed only once).
/// A fingerprint length mismatch raises ``ValueError``. Two all-zero
/// fingerprints compare as identical (1.0), matching the Rust API.
///
///     db_fps = [mol.ecfp4() for mol in library]   # compute once
///     for query in queries:
///         hits = chematic.nearest_neighbors_from_fp(query.ecfp4(), db_fps, k=10)
///         for idx, score in hits:
///             print(library_smiles[idx], score)
#[pyfunction]
#[pyo3(signature = (query_fp, db_fps, k = 10))]
fn nearest_neighbors_from_fp(
    query_fp: &[u8],
    db_fps: Vec<Bound<'_, PyAny>>,
    k: usize,
) -> PyResult<Vec<(usize, f64)>> {
    if k == 0 || db_fps.is_empty() {
        return Ok(Vec::new());
    }
    let db_fps = crate::fingerprint_similarity::fp_list(db_fps, "db_fps")?;
    let db_fps: Vec<&[u8]> = db_fps.iter().map(|f| f.as_slice()).collect();

    let query_popcount = crate::fingerprint_similarity::popcount(query_fp);
    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .enumerate()
        .map(|(index, fp)| {
            crate::fingerprint_similarity::tanimoto_bytes_with_counts(
                query_fp,
                fp,
                query_popcount,
                crate::fingerprint_similarity::popcount(fp),
            )
            .map(|score| (index, score))
            .map_err(|error| {
                PyValueError::new_err(format!(
                    "fingerprint length mismatch at db_fps[{index}]: {error}"
                ))
            })
        })
        .collect::<PyResult<_>>()?;
    scores.retain(|(_, score)| *score > 0.0);
    scores.sort_unstable_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    scores.truncate(k);
    Ok(scores)
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
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySmilesBatchStream>()?;
    register_pyfunctions!(
        m,
        from_smiles,
        parse_cdjson,
        write_cdjson_graph,
        canonicalize_smiles_batch_json,
        from_cxsmiles,
        from_condensed,
        from_mol_block,
        from_mol_block_with_coords,
        from_mol_block_with_diagnostics,
        parse_sdf_with_coords,
        from_cml,
        from_cml_strict,
        from_mrv_block,
        from_mrv_block_with_coords,
        from_cjson,
        from_moljson,
        from_cdxml,
        parse_cdxml_document_json,
        edit_cdxml_document_json,
        semantic_model_json,
        semantic_apply_json_command,
        semantic_expand_json,
        semantic_expand_json_with_limits,
        from_mol_v3000,
        from_mol_v3000_with_coords,
        from_mol_v3000_with_diagnostics,
        from_mol2,
        from_hin,
        parse_hin_record_json,
        roundtrip_hin_record,
        parse_cjson_record_json,
        roundtrip_cjson_record,
        parse_mol2_record_json,
        convert_with_report_json,
        roundtrip_mol2_record,
        from_pdbqt,
        from_gjf,
        parse_gaussian_log,
        write_gjf,
        parse_cif,
        write_cif_checked,
        is_valid_smiles,
        is_valid_smarts,
        from_inchi,
        from_pdb,
        from_pdb_strict,
        rdkit_from_pdb_block,
        rdkit_pdb_block_to_smiles,
        rdkit_from_xyz_block,
        rdkit_from_mol2_block,
        rdkit_mol2_block_to_smiles,
        from_xyz,
        from_extxyz,
        from_extxyz_all,
        to_extxyz,
        parse_formula,
        mol_hash,
        are_identical,
        stable_are_identical,
        write_reaction,
        from_rxn_file,
        to_rxn_file,
        from_rxn_document_json,
        to_rxn_document_json,
        nearest_neighbors_from_fp,
        parse_smi_file,
        write_smi_file,
        parse_ccd_json,
        from_ccd,
        parse_mmcif,
        parse_gro,
        parse_rdfile,
        write_rdfile,
        rdfile_to_smiles,
        cartesian_to_zmatrix,
        write_gzmat,
        parse_fhiaims,
        write_fhiaims,
        parse_dmol3,
        write_dmol3,
        select_pdb_atoms,
        parse_mopac_cartesian,
        write_mopac_cartesian,
        parse_exyz,
        write_exyz,
        parse_turbomole,
        write_turbomole,
        parse_mmcif_residues,
        mmcif_to_pdb_checked,
        protein_from_pdb_json,
        protein_from_mmcif_json,
        write_mmcif,
        parse_pqr,
        write_pqr,
        pqr_to_pdbqt,
        infer_element,
        parse_orca_input,
        write_orca_input,
        parse_orca_output,
        parse_qcschema_molecule,
        write_qcschema_molecule,
        chematic_to_qc_molecule,
        qc_molecule_to_chematic,
        parse_atomic_input,
        write_atomic_input,
        parse_atomic_result,
        write_atomic_result,
    )
}
