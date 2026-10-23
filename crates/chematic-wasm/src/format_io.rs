//! WASM bindings for the 7 file-format modules `chematic-mol` gained in
//! v0.17.0: Gaussian Cube, OpenDX, mmCIF, PQR, QCSchema, ORCA input/output,
//! LAMMPS data/dump. Follows [`crate::mol_io`]'s conventions: `Result<T,
//! JsValue>` for fallible operations (`.map_err(|e| JsValue::from_str(&e.to_string()))`
//! or an equivalent hand-built message), structured multi-field results
//! serialized to JSON strings via `serde_json` rather than bespoke
//! `#[wasm_bindgen]` structs.
//!
//! ## None of these formats carry a bond table
//!
//! mmCIF/PQR/ORCA/Cube/QCSchema (without an explicit `connectivity` list)
//! never perceive bonds -- that is a deliberate scope boundary of each
//! `chematic-mol` module (see their own doc comments), not something this
//! binding layer should silently paper over by fabricating one. So instead
//! of embedding a derived canonical-SMILES or MOL-block string into a
//! "molecule + metadata" JSON result (which would require inventing bonds
//! that were never actually perceived), this file mirrors each Rust
//! module's own `to_molecule() -> (Molecule, coords)` split: a `mol_from_*`
//! function returns a topology-only [`MolHandle`] (heavy atoms, no bonds),
//! and a companion `*_coords_json` function returns the coordinates in the
//! exact same atom order -- the same pattern [`crate::mol_io::mol_from_pdb`]
//! / [`crate::mol_io::pdb_coords_json`] already established for PDB (which
//! has the same "coordinates travel separately from topology" shape).
//!
//! ## Parse-limits types are not exposed to JS
//!
//! `chematic-mol`'s `*ParseLimits` types (`CubeParseLimits`,
//! `MmcifParseLimits`, `PqrParseLimits`, `OpenDxParseLimits`) all have a
//! `Default` impl covering the common case, and this crate already enforces
//! its own DoS-oriented input-size cap ([`crate::WASM_MAX_INPUT_BYTES`]) at
//! the WASM boundary before any parse is attempted (same as every other
//! function in [`crate::mol_io`]). Rather than inventing a JS-facing limits
//! object with no precedent anywhere else in this crate, every parse
//! function below simply uses each type's own `Default` limits.
//!
//! ## Grid JSON is a full round trip; typed-array siblings avoid it
//!
//! [`chematic_mol::VolumetricGrid::values`] can be very large (Cube/OpenDX
//! routinely carry hundreds of thousands to millions of voxels). At the time
//! `cube_grid_json`/`opendx_grid_json` were first added, this crate had no
//! `js_sys::Float64Array`/typed-array precedent anywhere (checked:
//! `mol_3d.rs`/`mol_depict.rs`/`mol_descriptors.rs`/`mol_edit.rs`/
//! `mol_fingerprints.rs`/`mol_io.rs`/`mol_reactions.rs` all used none), so
//! those two functions serialize `values` as a plain JSON number array -- a
//! real, disclosed perf cost (a full parse + JSON-string materialization +
//! JS-side `JSON.parse` of every voxel value). The "typed-array accessors"
//! section near the end of this file (`cube_values_f64`/`cube_shape_u32`/
//! `opendx_values_f64`/`opendx_shape_u32`/`lammps_dump_rows_f64`/
//! `lammps_dump_cartesian_positions_f64`) is that deferred follow-up:
//! `js_sys::Float64Array`/`Uint32Array`-returning siblings that avoid the
//! JSON round trip for the specific large, purely-numeric payloads where it
//! matters most, added ADDITIVELY alongside the original JSON functions
//! (which stay, unchanged, for every other field/format). This is now the
//! crate's first `js_sys` typed-array precedent -- `js-sys` was promoted
//! from a transitive dependency (already pulled in by `web-sys`) to a
//! direct one in `Cargo.toml` to allow it.

use crate::{
    MolHandle, WASM_MAX_ATOMS, WASM_MAX_BATCH_ITEMS, WASM_MAX_INPUT_BYTES,
    WASM_MAX_JSON_STRING_BYTES, WASM_MAX_OUTPUT_BYTES,
};
use wasm_bindgen::prelude::*;

mod lammps;
pub use lammps::*;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn check_input_len(label: &str, input: &str) -> Result<(), JsValue> {
    if input.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            input.len()
        )));
    }
    Ok(())
}

fn check_json_len(label: &str, input: &str) -> Result<(), JsValue> {
    if input.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            input.len()
        )));
    }
    Ok(())
}

/// Parse an MDL RXN V2000 file into the typed reaction-document JSON
/// contract shared with the Rust and Python bindings.
#[wasm_bindgen]
pub fn rxn_document_from_rxn(text: &str) -> Result<String, JsValue> {
    check_input_len("RXN input", text)?;
    let document = chematic_mol::parse_rxn_document(text)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    serde_json::to_string(&document).map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Write typed reaction-document JSON as MDL RXN V2000. Unsupported rich
/// fields return an error instead of being silently discarded.
#[wasm_bindgen]
pub fn rxn_document_to_rxn(document_json: &str) -> Result<String, JsValue> {
    check_json_len("RXN document JSON", document_json)?;
    let document = chematic_rxn::ReactionDocument::from_json_str(document_json)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    chematic_mol::write_rxn_document(&document)
        .map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Convert a rich reaction document to legacy RXN V2000 with structured loss
/// reporting. A lossy projection is never returned as if it were complete.
#[wasm_bindgen]
pub fn reaction_document_to_rxn_v1(document_json: &str) -> Result<String, JsValue> {
    check_json_len("reaction document JSON", document_json)?;
    let document =
        chematic_rxn::ReactionDocument::from_json_str(document_json).map_err(|error| {
            JsValue::from_str(&binding_error_json(
                "malformed_input",
                "/",
                error.to_string(),
            ))
        })?;
    chematic_mol::write_rxn_document(&document).map_err(|error| {
        let code = if matches!(
            error,
            chematic_mol::rxn::RxnDocumentError::Document(
                chematic_rxn::ReactionDocumentError::Losses(_)
            )
        ) {
            "lossy_conversion"
        } else {
            "unsupported_construct"
        };
        JsValue::from_str(&binding_error_json(code, "/", error.to_string()))
    })
}

fn binding_error_json(code: &str, path: &str, message: impl Into<String>) -> String {
    serde_json::json!({
        "ok": false,
        "error": { "code": code, "path": path, "message": message.into() }
    })
    .to_string()
}

/// Validate and deterministically serialize a rich reaction document JSON.
///
/// This is the versioned JSON boundary for WASM consumers. It preserves IDs,
/// metadata, and provenance, and rejects malformed documents before returning
/// JSON. Use [`edit_reaction_document_json_v1`] for bounded ID-addressed edits.
#[wasm_bindgen]
pub fn reaction_document_json_v1(document_json: &str) -> Result<String, JsValue> {
    check_json_len("reaction document JSON", document_json)?;
    let document =
        chematic_rxn::ReactionDocument::from_json_str(document_json).map_err(|error| {
            JsValue::from_str(&binding_error_json(
                "malformed_input",
                "/",
                error.to_string(),
            ))
        })?;
    serde_json::to_string(&document).map_err(|error| {
        JsValue::from_str(&binding_error_json(
            "serialization_error",
            "/",
            error.to_string(),
        ))
    })
}

/// Apply a bounded, stable-ID reaction-document edit and return canonical JSON.
#[wasm_bindgen]
pub fn edit_reaction_document_json_v1(
    document_json: &str,
    edit_json: &str,
) -> Result<String, JsValue> {
    check_json_len("reaction document JSON", document_json)?;
    check_json_len("reaction edit JSON", edit_json)?;
    let document =
        chematic_rxn::ReactionDocument::from_json_str(document_json).map_err(|error| {
            JsValue::from_str(&binding_error_json(
                "malformed_input",
                "/",
                error.to_string(),
            ))
        })?;
    let edited = document.apply_json_edit(edit_json).map_err(|error| {
        JsValue::from_str(&binding_error_json(
            "unsupported_construct",
            "/edit",
            error.to_string(),
        ))
    })?;
    serde_json::to_string(&edited).map_err(|error| {
        JsValue::from_str(&binding_error_json(
            "serialization_error",
            "/",
            error.to_string(),
        ))
    })
}

/// Return a loss-preserving, versioned JSON envelope for a CDXML document.
///
/// `source` is retained for exact re-serialization; `document` is the stable
/// structural summary used by editors. Unknown objects remain in the summary
/// and are listed under `document.diagnostics`.
#[wasm_bindgen]
pub fn cdxml_document_json_v1(cdxml: &str) -> Result<String, JsValue> {
    check_input_len("CDXML input", cdxml)?;
    let document = chematic_mol::CdxmlDocument::parse_with_limits(
        cdxml,
        &chematic_mol::CdxmlParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            ..Default::default()
        },
    )
    .map_err(|error| {
        let code = if matches!(error, chematic_mol::CdxmlError::ResourceLimit { .. }) {
            "resource_limit"
        } else {
            "malformed_input"
        };
        JsValue::from_str(&binding_error_json(code, "/", error.to_string()))
    })?;
    let result = serde_json::json!({
        "schema": "chematic.cdxml-binding.v1",
        "source": document.write(),
        "document": document.to_json(),
    });
    let json = result.to_string();
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(&binding_error_json(
            "resource_limit",
            "/",
            "CDXML binding JSON exceeds output limit",
        )));
    }
    Ok(json)
}

/// Project a CDXML document to its first molecular fragment only when no
/// presentation data would be lost. Unknown/presentation objects are reported
/// as an explicit lossy conversion instead of being silently dropped.
#[wasm_bindgen]
pub fn cdxml_document_projection_json_v1(cdxml: &str) -> Result<String, JsValue> {
    check_input_len("CDXML input", cdxml)?;
    let document = chematic_mol::CdxmlDocument::parse_with_limits(
        cdxml,
        &chematic_mol::CdxmlParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            ..Default::default()
        },
    )
    .map_err(|error| {
        let code = if matches!(error, chematic_mol::CdxmlError::ResourceLimit { .. }) {
            "resource_limit"
        } else {
            "malformed_input"
        };
        JsValue::from_str(&binding_error_json(code, "/", error.to_string()))
    })?;
    let diagnostics = document.diagnostics();
    if !diagnostics.is_empty() {
        return Err(JsValue::from_str(&binding_error_json(
            "lossy_conversion",
            "/document/diagnostics",
            "CDXML presentation objects cannot be represented by the molecule projection",
        )));
    }
    let (molecule, _) = chematic_mol::parse_cdxml_with_limits(
        cdxml,
        &chematic_mol::CdxmlParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            ..Default::default()
        },
    )
    .map_err(|error| {
        JsValue::from_str(&binding_error_json(
            "malformed_input",
            "/",
            error.to_string(),
        ))
    })?;
    Ok(serde_json::json!({
        "ok": true,
        "smiles": chematic_smiles::canonical_smiles(&molecule),
    })
    .to_string())
}

/// Validate a versioned CDXML JSON envelope and serialize its exact source.
#[wasm_bindgen]
pub fn cdxml_document_from_json_v1(document_json: &str) -> Result<String, JsValue> {
    check_json_len("CDXML binding JSON", document_json)?;
    let value: serde_json::Value = serde_json::from_str(document_json).map_err(|error| {
        JsValue::from_str(&binding_error_json(
            "malformed_input",
            "/",
            error.to_string(),
        ))
    })?;
    if value.get("schema").and_then(serde_json::Value::as_str) != Some("chematic.cdxml-binding.v1")
    {
        return Err(JsValue::from_str(&binding_error_json(
            "unsupported_construct",
            "/schema",
            "unsupported CDXML binding schema",
        )));
    }
    let source = value
        .get("source")
        .and_then(serde_json::Value::as_str)
        .ok_or_else(|| {
            JsValue::from_str(&binding_error_json(
                "malformed_input",
                "/source",
                "source is required",
            ))
        })?;
    check_input_len("CDXML source", source)?;
    let document = chematic_mol::CdxmlDocument::parse_with_limits(
        source,
        &chematic_mol::CdxmlParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            ..Default::default()
        },
    )
    .map_err(|error| {
        JsValue::from_str(&binding_error_json(
            "malformed_input",
            "/source",
            error.to_string(),
        ))
    })?;
    Ok(document.write())
}

/// Apply a bounded CDXML edit to a versioned JSON envelope and return the
/// updated envelope. The source is reparsed after editing, so paths and
/// diagnostics cannot drift from the returned document summary.
#[wasm_bindgen]
pub fn edit_cdxml_document_json_v1(
    document_json: &str,
    edit_json: &str,
) -> Result<String, JsValue> {
    let source = cdxml_document_from_json_v1(document_json)?;
    let edited = crate::edit_cdxml_document_json(&source, edit_json)?;
    cdxml_document_json_v1(&edited)
}

fn wasm_orca_input_limits() -> chematic_mol::OrcaInputParseLimits {
    chematic_mol::OrcaInputParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_line_bytes: WASM_MAX_INPUT_BYTES,
        max_lines: WASM_MAX_BATCH_ITEMS,
        max_keywords: WASM_MAX_BATCH_ITEMS,
        max_blocks: WASM_MAX_BATCH_ITEMS,
        max_block_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
    }
}

fn wasm_qcschema_limits() -> chematic_mol::QcSchemaParseLimits {
    chematic_mol::QcSchemaParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_json_depth: 64,
        max_array_items: WASM_MAX_BATCH_ITEMS,
        max_string_bytes: WASM_MAX_JSON_STRING_BYTES,
    }
}

fn wasm_cube_limits() -> chematic_mol::CubeParseLimits {
    chematic_mol::CubeParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        max_grid_points: 1_000_000,
    }
}

fn wasm_opendx_limits() -> chematic_mol::OpenDxParseLimits {
    chematic_mol::OpenDxParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_grid_points: 1_000_000,
    }
}

fn wasm_mmcif_limits() -> chematic_mol::MmcifParseLimits {
    chematic_mol::MmcifParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        max_line_len: WASM_MAX_INPUT_BYTES,
    }
}

fn wasm_pqr_limits() -> chematic_mol::PqrParseLimits {
    chematic_mol::PqrParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        max_line_len: WASM_MAX_INPUT_BYTES,
    }
}

fn wasm_lammps_limits() -> chematic_mol::LammpsDataParseLimits {
    chematic_mol::LammpsDataParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_line_bytes: WASM_MAX_INPUT_BYTES,
        max_header_counts: WASM_MAX_BATCH_ITEMS,
        max_masses: WASM_MAX_BATCH_ITEMS,
        max_atoms: WASM_MAX_ATOMS,
        max_velocities: WASM_MAX_ATOMS,
        max_bonds: WASM_MAX_ATOMS,
        max_opaque_section_bytes: WASM_MAX_INPUT_BYTES,
        max_sections: WASM_MAX_BATCH_ITEMS,
    }
}

fn wasm_lammps_dump_limits() -> chematic_mol::LammpsDumpParseLimits {
    chematic_mol::LammpsDumpParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_line_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms_per_frame: WASM_MAX_ATOMS,
        max_columns: WASM_MAX_BATCH_ITEMS,
        max_frames: WASM_MAX_BATCH_ITEMS,
    }
}

fn mol_handle_from_molecule(mol: chematic_core::Molecule) -> Result<MolHandle, JsValue> {
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {WASM_MAX_ATOMS} atoms)"
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

fn coords_tuples_to_json(coords: &[(f64, f64, f64)]) -> serde_json::Value {
    serde_json::Value::Array(
        coords
            .iter()
            .map(|(x, y, z)| serde_json::json!([x, y, z]))
            .collect(),
    )
}

fn to_json_string(v: &serde_json::Value) -> Result<String, JsValue> {
    let output = serde_json::to_string(v).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if output.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "JSON output exceeds maximum size ({} > {WASM_MAX_OUTPUT_BYTES} bytes)",
            output.len()
        )));
    }
    Ok(output)
}

fn parse_json_value(label: &str, text: &str) -> Result<serde_json::Value, JsValue> {
    check_json_len(label, text)?;
    serde_json::from_str(text)
        .map_err(|e| JsValue::from_str(&format!("{label}: invalid JSON: {e}")))
}

// Small `serde_json::Value` field-extraction helpers shared by every
// hand-built record type below (none of `chematic-mol`'s new format structs
// derive `Serialize`/`Deserialize` -- see qcschema.rs's module docs on why
// this crate has no blanket serde derive -- so every JSON<->struct
// conversion here is written by hand against `serde_json::Value`, matching
// `crate::mol_io::extxyz_frame_from_json_args`'s existing precedent for the
// same situation).
fn jget<'a>(v: &'a serde_json::Value, key: &str) -> Result<&'a serde_json::Value, String> {
    v.get(key).ok_or_else(|| format!("missing field '{key}'"))
}
fn jstr(v: &serde_json::Value, key: &str) -> Result<String, String> {
    jget(v, key)?
        .as_str()
        .map(str::to_string)
        .ok_or_else(|| format!("'{key}' must be a string"))
}
fn ji64(v: &serde_json::Value, key: &str) -> Result<i64, String> {
    jget(v, key)?
        .as_i64()
        .ok_or_else(|| format!("'{key}' must be an integer"))
}
fn jf64(v: &serde_json::Value, key: &str) -> Result<f64, String> {
    jget(v, key)?
        .as_f64()
        .ok_or_else(|| format!("'{key}' must be a number"))
}
fn jf64_or(v: &serde_json::Value, key: &str, default: f64) -> f64 {
    v.get(key).and_then(|x| x.as_f64()).unwrap_or(default)
}
fn ji64_or(v: &serde_json::Value, key: &str, default: i64) -> i64 {
    v.get(key).and_then(|x| x.as_i64()).unwrap_or(default)
}
fn jopt_i64(v: &serde_json::Value, key: &str) -> Option<i64> {
    v.get(key).and_then(|x| x.as_i64())
}
fn jopt_str(v: &serde_json::Value, key: &str) -> Option<String> {
    v.get(key).and_then(|x| x.as_str()).map(str::to_string)
}
fn jarr<'a>(v: &'a serde_json::Value, key: &str) -> Result<&'a Vec<serde_json::Value>, String> {
    bounded_array(jget(v, key)?, key)
}

fn bounded_array<'a>(
    value: &'a serde_json::Value,
    label: &str,
) -> Result<&'a Vec<serde_json::Value>, String> {
    let array = value
        .as_array()
        .ok_or_else(|| format!("'{label}' must be an array"))?;
    if array.len() > WASM_MAX_BATCH_ITEMS {
        return Err(format!(
            "{label} exceeds maximum item count ({} > {})",
            array.len(),
            WASM_MAX_BATCH_ITEMS
        ));
    }
    Ok(array)
}

/// A single-character JSON field (or `null`) -- used by mmCIF/PQR's
/// `alt_loc`/`icode` fields.
fn json_to_opt_char(v: &serde_json::Value, field: &str) -> Result<Option<char>, String> {
    match v {
        serde_json::Value::Null => Ok(None),
        serde_json::Value::String(s) => {
            let mut chars = s.chars();
            let c = chars
                .next()
                .ok_or_else(|| format!("'{field}' must be a single character, got empty string"))?;
            if chars.next().is_some() {
                return Err(format!("'{field}' must be a single character"));
            }
            Ok(Some(c))
        }
        _ => Err(format!("'{field}' must be a string or null")),
    }
}
fn jopt_char(v: &serde_json::Value, key: &str) -> Result<Option<char>, String> {
    match v.get(key) {
        Some(x) => json_to_opt_char(x, key),
        None => Ok(None),
    }
}
fn opt_char_json(c: Option<char>) -> serde_json::Value {
    match c {
        Some(c) => serde_json::json!(c.to_string()),
        None => serde_json::Value::Null,
    }
}

fn element_from_symbol(sym: &str) -> Result<chematic_core::Element, String> {
    chematic_core::Element::from_symbol(sym).ok_or_else(|| format!("unknown element '{sym}'"))
}

// ===========================================================================
// Gaussian Cube + OpenDX (share chematic_mol::VolumetricGrid)
// ===========================================================================

fn grid_units_str(u: chematic_mol::GridUnits) -> &'static str {
    match u {
        chematic_mol::GridUnits::Bohr => "bohr",
        chematic_mol::GridUnits::Angstrom => "angstrom",
    }
}

fn grid_units_from_str(s: &str) -> Result<chematic_mol::GridUnits, String> {
    match s {
        "bohr" => Ok(chematic_mol::GridUnits::Bohr),
        "angstrom" => Ok(chematic_mol::GridUnits::Angstrom),
        other => Err(format!(
            "unknown units '{other}' (expected 'bohr' or 'angstrom')"
        )),
    }
}

fn grid_atom_to_json(a: &chematic_mol::GridAtom) -> serde_json::Value {
    serde_json::json!({
        "element": a.element.symbol(),
        "charge": a.charge,
        "position": a.position,
    })
}

fn grid_atom_from_json(v: &serde_json::Value) -> Result<chematic_mol::GridAtom, String> {
    let element = element_from_symbol(&jstr(v, "element")?)?;
    let charge = jf64(v, "charge")?;
    let position: [f64; 3] = serde_json::from_value(jget(v, "position")?.clone())
        .map_err(|e| format!("invalid 'position': {e}"))?;
    Ok(chematic_mol::GridAtom {
        element,
        charge,
        position,
    })
}

fn volumetric_grid_to_json(g: &chematic_mol::VolumetricGrid) -> serde_json::Value {
    serde_json::json!({
        "origin": g.origin,
        "axes": g.axes,
        "shape": g.shape,
        "values": g.values,
        "atoms": g.atoms.iter().map(grid_atom_to_json).collect::<Vec<_>>(),
        "units": grid_units_str(g.units),
    })
}

fn volumetric_grid_from_json_str(text: &str) -> Result<chematic_mol::VolumetricGrid, String> {
    let v: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("invalid grid JSON: {e}"))?;
    let origin: [f64; 3] = serde_json::from_value(jget(&v, "origin")?.clone())
        .map_err(|e| format!("invalid 'origin': {e}"))?;
    let axes: [[f64; 3]; 3] = serde_json::from_value(jget(&v, "axes")?.clone())
        .map_err(|e| format!("invalid 'axes': {e}"))?;
    let shape: [usize; 3] = serde_json::from_value(jget(&v, "shape")?.clone())
        .map_err(|e| format!("invalid 'shape': {e}"))?;
    let values: Vec<f64> = serde_json::from_value(jget(&v, "values")?.clone())
        .map_err(|e| format!("invalid 'values': {e}"))?;
    let atoms = match v.get("atoms") {
        Some(serde_json::Value::Array(arr)) => arr
            .iter()
            .map(grid_atom_from_json)
            .collect::<Result<Vec<_>, _>>()?,
        _ => Vec::new(),
    };
    let units = grid_units_from_str(v.get("units").and_then(|x| x.as_str()).unwrap_or("bohr"))?;
    Ok(chematic_mol::VolumetricGrid {
        origin,
        axes,
        shape,
        values,
        atoms,
        units,
    })
}

/// Parse a Gaussian Cube file and return a `MolHandle` (topology only --
/// element list, no bonds; Cube carries no bond table). Use
/// [`cube_grid_json`] to recover coordinates, the scalar field, and the
/// grid geometry.
#[wasm_bindgen]
pub fn mol_from_cube(text: &str) -> Result<MolHandle, JsValue> {
    check_input_len("cube input", text)?;
    let grid = chematic_mol::parse_cube_with_limits(text, &wasm_cube_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (mol, _coords) = grid.to_molecule();
    mol_handle_from_molecule(mol)
}

/// Parse a Gaussian Cube file and return its full [`chematic_mol::VolumetricGrid`]
/// as JSON: `{"origin":[x,y,z],"axes":[[..],[..],[..]],"shape":[nx,ny,nz],
/// "values":[...flat, row-major third-axis-fastest...],
/// "atoms":[{"element":"C","charge":6.0,"position":[x,y,z]}],
/// "units":"bohr"|"angstrom"}`. See module docs for the perf tradeoff of a
/// full `values` JSON round trip on a large grid.
#[wasm_bindgen]
pub fn cube_grid_json(text: &str) -> Result<String, JsValue> {
    check_input_len("cube input", text)?;
    let grid = chematic_mol::parse_cube_with_limits(text, &wasm_cube_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    to_json_string(&volumetric_grid_to_json(&grid))
}

/// Write a grid (in the JSON shape [`cube_grid_json`] returns) as a
/// Gaussian Cube file.
#[wasm_bindgen]
pub fn write_cube_json(grid_json: &str) -> Result<String, JsValue> {
    check_json_len("cube grid JSON", grid_json)?;
    let grid = volumetric_grid_from_json_str(grid_json).map_err(|e| JsValue::from_str(&e))?;
    chematic_mol::write_cube(&grid).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Parse an OpenDX (APBS scalar-field subset) file and return its full
/// [`chematic_mol::VolumetricGrid`] as JSON (same shape as
/// [`cube_grid_json`]; `atoms` is always empty -- OpenDX has no atom
/// section). No `mol_from_opendx` is provided: an OpenDX grid never carries
/// atoms, so a `MolHandle` from one would always be empty and is not a
/// useful binding.
#[wasm_bindgen]
pub fn opendx_grid_json(text: &str) -> Result<String, JsValue> {
    check_input_len("OpenDX input", text)?;
    let grid = chematic_mol::parse_opendx_with_limits(text, &wasm_opendx_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    to_json_string(&volumetric_grid_to_json(&grid))
}

/// Write a grid as an OpenDX file. Fails closed for a
/// [`chematic_mol::GridUnits::Bohr`]-tagged grid (OpenDX has no unit tag of
/// its own and is universally read back as Ångström -- see
/// `chematic_mol::opendx`'s module docs) and for a grid carrying any atoms
/// (OpenDX has no atom section). Use [`write_opendx_lossy_json`] to opt
/// into an explicit Bohr->Ångström conversion instead of failing.
#[wasm_bindgen]
pub fn write_opendx_json(grid_json: &str) -> Result<String, JsValue> {
    check_json_len("OpenDX grid JSON", grid_json)?;
    let grid = volumetric_grid_from_json_str(grid_json).map_err(|e| JsValue::from_str(&e))?;
    chematic_mol::write_opendx(&grid).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Like [`write_opendx_json`], but a [`chematic_mol::GridUnits::Bohr`] grid
/// has its `origin`/`axes` explicitly converted to Ångström rather than
/// rejected (`values` -- the scalar-field samples themselves -- are never
/// rescaled; see `write_opendx_lossy`'s doc comment). Still fails for a
/// grid carrying any atoms.
#[wasm_bindgen]
pub fn write_opendx_lossy_json(grid_json: &str) -> Result<String, JsValue> {
    check_json_len("OpenDX grid JSON", grid_json)?;
    let grid = volumetric_grid_from_json_str(grid_json).map_err(|e| JsValue::from_str(&e))?;
    chematic_mol::write_opendx_lossy(&grid).map_err(|e| JsValue::from_str(&e.to_string()))
}

// ===========================================================================
// mmCIF
// ===========================================================================

fn unit_cell_to_json(c: &chematic_mol::UnitCell) -> serde_json::Value {
    serde_json::json!({
        "a": c.a, "b": c.b, "c": c.c,
        "alpha": c.alpha, "beta": c.beta, "gamma": c.gamma,
    })
}

fn unit_cell_from_json(v: &serde_json::Value) -> Result<chematic_mol::UnitCell, String> {
    Ok(chematic_mol::UnitCell {
        a: jf64(v, "a")?,
        b: jf64(v, "b")?,
        c: jf64(v, "c")?,
        alpha: jf64(v, "alpha")?,
        beta: jf64(v, "beta")?,
        gamma: jf64(v, "gamma")?,
    })
}

fn mmcif_atom_to_json(a: &chematic_mol::MmcifAtomRecord) -> serde_json::Value {
    serde_json::json!({
        "group_pdb": a.group_pdb,
        "serial": a.serial,
        "element": a.element.symbol(),
        "atom_name": a.atom_name,
        "alt_loc": opt_char_json(a.alt_loc),
        "res_name": a.res_name,
        "chain_id": a.chain_id,
        "res_seq": a.res_seq,
        "label_seq_id": a.label_seq_id,
        "icode": opt_char_json(a.icode),
        "x": a.x, "y": a.y, "z": a.z,
        "occupancy": a.occupancy,
        "b_iso": a.b_iso,
        "formal_charge": a.formal_charge,
        "entity_id": a.entity_id,
        "model_num": a.model_num,
    })
}

fn mmcif_atom_from_json(v: &serde_json::Value) -> Result<chematic_mol::MmcifAtomRecord, String> {
    Ok(chematic_mol::MmcifAtomRecord {
        group_pdb: jstr(v, "group_pdb")?,
        serial: ji64(v, "serial")?,
        element: element_from_symbol(&jstr(v, "element")?)?,
        atom_name: jstr(v, "atom_name")?,
        alt_loc: jopt_char(v, "alt_loc")?,
        res_name: jstr(v, "res_name")?,
        chain_id: jstr(v, "chain_id")?,
        res_seq: ji64(v, "res_seq")?,
        label_seq_id: jopt_i64(v, "label_seq_id"),
        icode: jopt_char(v, "icode")?,
        x: jf64(v, "x")?,
        y: jf64(v, "y")?,
        z: jf64(v, "z")?,
        occupancy: jf64_or(v, "occupancy", 1.0),
        b_iso: jf64_or(v, "b_iso", 0.0),
        formal_charge: jopt_i64(v, "formal_charge").map(|x| x as i32),
        entity_id: jopt_str(v, "entity_id"),
        model_num: ji64_or(v, "model_num", 1) as i32,
    })
}

/// Parse an mmCIF file and return a `MolHandle` (topology only -- element
/// list, no bonds; mmCIF's `_atom_site` category carries no connectivity).
/// Includes every model's atoms if the file has more than one -- use
/// [`mmcif_to_json`] to get each atom's `model_num` for filtering. Use
/// [`mmcif_coords_json`] to recover coordinates in the same atom order.
#[wasm_bindgen]
pub fn mol_from_mmcif(text: &str) -> Result<MolHandle, JsValue> {
    check_input_len("mmCIF input", text)?;
    let result = chematic_mol::parse_mmcif_with_limits(text, &wasm_mmcif_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (mol, _coords) = result.to_molecule();
    mol_handle_from_molecule(mol)
}

/// Cartesian coordinates from an mmCIF file, in the SAME atom order
/// [`mol_from_mmcif`] returns topology for. Returns JSON `[[x,y,z],...]`
/// (Å).
#[wasm_bindgen]
pub fn mmcif_coords_json(text: &str) -> Result<String, JsValue> {
    check_input_len("mmCIF input", text)?;
    let result = chematic_mol::parse_mmcif_with_limits(text, &wasm_mmcif_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (_mol, coords) = result.to_molecule();
    to_json_string(&coords_tuples_to_json(&coords))
}

/// Parse an mmCIF file and return every `_atom_site` field (occupancy,
/// B-factor, chain/residue bookkeeping, formal charge, model number, ...),
/// the unit cell, space group, and any loop column this reader saw but
/// does not model, as JSON: `{"atoms":[{...}],"cell":{...}|null,
/// "space_group":"..."|null,"unhandled_columns":[...]}`. See
/// [`chematic_mol::MmcifAtomRecord`]'s doc comment for each atom field's
/// exact source column and defaulting rule.
#[wasm_bindgen]
pub fn mmcif_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("mmCIF input", text)?;
    let result = chematic_mol::parse_mmcif_with_limits(text, &wasm_mmcif_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let out = serde_json::json!({
        "atoms": result.atoms.iter().map(mmcif_atom_to_json).collect::<Vec<_>>(),
        "cell": result.cell.as_ref().map(unit_cell_to_json),
        "space_group": result.space_group,
        "unhandled_columns": result.unhandled_columns,
    });
    to_json_string(&out)
}

fn protein_projection_to_json(projection: &chematic_mol::ProteinProjection) -> serde_json::Value {
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

/// Amino-acid-only PDB projection as stable JSON.
#[wasm_bindgen]
pub fn protein_from_pdb_json(text: &str) -> Result<String, JsValue> {
    check_input_len("PDB input", text)?;
    let limits = chematic_3d::PdbParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_line_bytes: WASM_MAX_JSON_STRING_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        max_models: 1_000,
    };
    let projection = chematic_mol::project_protein_from_pdb_with_limits(text, &limits)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    to_json_string(&protein_projection_to_json(&projection))
}

/// Amino-acid-only mmCIF projection in the same schema as
/// [`protein_from_pdb_json`].
#[wasm_bindgen]
pub fn protein_from_mmcif_json(text: &str) -> Result<String, JsValue> {
    check_input_len("mmCIF input", text)?;
    let projection =
        chematic_mol::project_protein_from_mmcif_with_limits(text, &wasm_mmcif_limits())
            .map_err(|error| JsValue::from_str(&error.to_string()))?;
    to_json_string(&protein_projection_to_json(&projection))
}

/// Write an mmCIF file from atom records in the JSON shape
/// [`mmcif_to_json`]'s `"atoms"` array uses (a full record per atom, not
/// just element+coordinates -- mmCIF has no equivalent of "build from a
/// bare `MolHandle`", since occupancy/B-factor/chain/residue fields have no
/// source in a plain [`MolHandle`]).
///
/// `cell_json`: `"null"` or `{"a":...,"b":...,"c":...,"alpha":...,"beta":...,"gamma":...}`.
/// `space_group`: pass `""` for none.
#[wasm_bindgen]
pub fn write_mmcif_json(
    records_json: &str,
    cell_json: &str,
    space_group: &str,
    data_block_name: &str,
) -> Result<String, JsValue> {
    let records = parse_json_value("mmCIF records JSON", records_json)?;
    let atoms: Vec<chematic_mol::MmcifAtomRecord> = bounded_array(&records, "records_json")
        .map_err(|e| JsValue::from_str(&e))?
        .iter()
        .map(mmcif_atom_from_json)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| JsValue::from_str(&e))?;

    let cell_value = parse_json_value("cell JSON", cell_json)?;
    let cell = if cell_value.is_null() {
        None
    } else {
        Some(unit_cell_from_json(&cell_value).map_err(|e| JsValue::from_str(&e))?)
    };
    let space_group_opt = if space_group.is_empty() {
        None
    } else {
        Some(space_group)
    };

    Ok(chematic_mol::write_mmcif(
        &atoms,
        cell.as_ref(),
        space_group_opt,
        data_block_name,
    ))
}

// ===========================================================================
// PQR
// ===========================================================================

mod pqr_checked;
pub use pqr_checked::*;

fn pqr_atom_to_json(a: &chematic_mol::PqrAtomRecord) -> serde_json::Value {
    serde_json::json!({
        "group_pdb": a.group_pdb,
        "serial": a.serial,
        "atom_name": a.atom_name,
        "res_name": a.res_name,
        "chain_id": a.chain_id,
        "res_seq": a.res_seq,
        "icode": opt_char_json(a.icode),
        "x": a.x, "y": a.y, "z": a.z,
        "charge": a.charge,
        "radius": a.radius,
        "element": a.element.symbol(),
    })
}

fn pqr_atom_from_json(v: &serde_json::Value) -> Result<chematic_mol::PqrAtomRecord, String> {
    Ok(chematic_mol::PqrAtomRecord {
        group_pdb: jstr(v, "group_pdb")?,
        serial: ji64(v, "serial")?,
        atom_name: jstr(v, "atom_name")?,
        res_name: jstr(v, "res_name")?,
        chain_id: jopt_str(v, "chain_id"),
        res_seq: ji64(v, "res_seq")?,
        icode: jopt_char(v, "icode")?,
        x: jf64(v, "x")?,
        y: jf64(v, "y")?,
        z: jf64(v, "z")?,
        charge: jf64(v, "charge")?,
        radius: jf64(v, "radius")?,
        element: element_from_symbol(&jstr(v, "element")?)?,
    })
}

/// Parse a PQR file and return a `MolHandle` (topology only -- element
/// list inferred per-atom, no bonds; PQR carries no connectivity). Use
/// [`pqr_coords_json`] to recover coordinates in the same atom order.
#[wasm_bindgen]
pub fn mol_from_pqr(text: &str) -> Result<MolHandle, JsValue> {
    check_input_len("PQR input", text)?;
    let result = chematic_mol::parse_pqr_with_limits(text, &wasm_pqr_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (mol, _coords) = result.to_molecule();
    mol_handle_from_molecule(mol)
}

/// Cartesian coordinates from a PQR file, in the SAME atom order
/// [`mol_from_pqr`] returns topology for. Returns JSON `[[x,y,z],...]` (Å).
#[wasm_bindgen]
pub fn pqr_coords_json(text: &str) -> Result<String, JsValue> {
    check_input_len("PQR input", text)?;
    let result = chematic_mol::parse_pqr_with_limits(text, &wasm_pqr_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (_mol, coords) = result.to_molecule();
    to_json_string(&coords_tuples_to_json(&coords))
}

/// Parse a PQR file and return every field (charge, radius, chain,
/// residue, inferred element, ...) as JSON: `{"atoms":[{...}]}`. See
/// [`chematic_mol::PqrAtomRecord`]'s doc comment for each field's meaning.
#[wasm_bindgen]
pub fn pqr_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("PQR input", text)?;
    let result = chematic_mol::parse_pqr_with_limits(text, &wasm_pqr_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let out = serde_json::json!({
        "atoms": result.atoms.iter().map(pqr_atom_to_json).collect::<Vec<_>>(),
    });
    to_json_string(&out)
}

/// Write a PQR file from atom records in the JSON shape [`pqr_to_json`]'s
/// `"atoms"` array uses. Each atom's `chain_id` independently controls
/// whether that line is written with or without the (optional) chain
/// column.
#[wasm_bindgen]
pub fn write_pqr_json(records_json: &str) -> Result<String, JsValue> {
    let records = parse_json_value("PQR records JSON", records_json)?;
    let atoms: Vec<chematic_mol::PqrAtomRecord> = bounded_array(&records, "records_json")
        .map_err(|e| JsValue::from_str(&e))?
        .iter()
        .map(pqr_atom_from_json)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| JsValue::from_str(&e))?;
    Ok(chematic_mol::write_pqr(&atoms))
}

/// Infer an element from a PQR atom name (see
/// [`chematic_mol::infer_element`]'s doc comment for the heuristic).
/// Returns `undefined` (JS) / `None` if no element could be inferred.
#[wasm_bindgen]
pub fn pqr_infer_element(group_pdb: &str, res_name: &str, atom_name: &str) -> Option<String> {
    chematic_mol::infer_element(group_pdb, res_name, atom_name).map(|e| e.symbol().to_string())
}

// ===========================================================================
// QCSchema
// ===========================================================================

/// Parse a QCSchema `qcschema_molecule` JSON document and return a
/// `MolHandle` (topology + `atomic_numbers`-derived isotopes -- no bonds
/// unless the document's optional `connectivity` list is present, in which
/// case those bond orders are mapped onto the nearest
/// [`chematic_core::BondOrder`]). Use [`qcschema_molecule_coords_json`] to
/// recover coordinates (converted Bohr -> Å) plus molecular
/// charge/multiplicity, in the same atom order.
#[wasm_bindgen]
pub fn mol_from_qcschema_molecule(json: &str) -> Result<MolHandle, JsValue> {
    check_input_len("QCSchema molecule JSON", json)?;
    let qc = chematic_mol::parse_qcschema_molecule_with_limits(json, &wasm_qcschema_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let view = chematic_mol::qc_molecule_to_chematic(&qc)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    mol_handle_from_molecule(view.molecule)
}

/// Coordinates (Å) plus molecular charge/multiplicity from a QCSchema
/// `qcschema_molecule` document, in the SAME atom order
/// [`mol_from_qcschema_molecule`] returns topology for. Returns JSON
/// `{"coords":[[x,y,z],...],"molecular_charge":0.0,"molecular_multiplicity":1}`.
#[wasm_bindgen]
pub fn qcschema_molecule_coords_json(json: &str) -> Result<String, JsValue> {
    check_input_len("QCSchema molecule JSON", json)?;
    let qc = chematic_mol::parse_qcschema_molecule_with_limits(json, &wasm_qcschema_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let view = chematic_mol::qc_molecule_to_chematic(&qc)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let coords: Vec<[f64; 3]> = (0..view.molecule.atom_count())
        .map(|i| {
            let p = view.coords.get(chematic_core::AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();
    let out = serde_json::json!({
        "coords": coords,
        "molecular_charge": view.molecular_charge,
        "molecular_multiplicity": view.molecular_multiplicity,
    });
    to_json_string(&out)
}

/// Serialize a `MolHandle` + coordinates (Å) + molecular charge/multiplicity
/// as a QCSchema `qcschema_molecule` JSON document (coordinates converted
/// to Bohr).
///
/// `coords_json`: `[[x,y,z],...]` (Å), same order and length as `mol`'s
/// atoms.
#[wasm_bindgen]
pub fn to_qcschema_molecule_json(
    mol: &MolHandle,
    coords_json: &str,
    charge: f64,
    multiplicity: i64,
) -> Result<String, JsValue> {
    check_json_len("coords_json", coords_json)?;
    let coords_flat: Vec<[f64; 3]> = serde_json::from_str(coords_json)
        .map_err(|e| JsValue::from_str(&format!("invalid coords_json: {e}")))?;
    let n = mol.inner.atom_count();
    if coords_flat.len() != n {
        return Err(JsValue::from_str(&format!(
            "coords_json has {} row(s), mol has {n} atom(s)",
            coords_flat.len()
        )));
    }
    let mut coords3d = chematic_core::Coords3D::new_zeroed(n);
    for (i, [x, y, z]) in coords_flat.into_iter().enumerate() {
        coords3d.set(
            chematic_core::AtomIdx(i as u32),
            chematic_core::Point3::new(x, y, z),
        );
    }
    let qc = chematic_mol::chematic_to_qc_molecule(&mol.inner, &coords3d, charge, multiplicity)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_mol::write_qcschema_molecule(&qc))
}

/// Parse a QCSchema `qcschema_input`/`qc_schema_input` JSON document
/// (molecule + driver + model + keywords) and re-emit it, validating and
/// canonicalizing field defaults in the process (e.g. a missing
/// `schema_name`/`schema_version` is filled in). Job-level fields
/// (`driver`, `model`, `keywords`, `protocols`, `extras`) are round-tripped
/// opaquely -- this binding validates/reformats the document; it does not
/// expose a separate JS-facing accessor for each field (out of scope for
/// this first pass, see module docs' "None of these formats carry a bond
/// table" section for the analogous molecule-centric scope choice made
/// elsewhere in this file).
#[wasm_bindgen]
pub fn qcschema_validate_atomic_input(json: &str) -> Result<String, JsValue> {
    check_input_len("QCSchema AtomicInput JSON", json)?;
    let input = chematic_mol::parse_atomic_input_with_limits(json, &wasm_qcschema_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_mol::write_atomic_input(&input))
}

/// Like [`qcschema_validate_atomic_input`], for a QCSchema
/// `qcschema_output`/`qc_schema_output` (`AtomicResult`) document.
#[wasm_bindgen]
pub fn qcschema_validate_atomic_result(json: &str) -> Result<String, JsValue> {
    check_input_len("QCSchema AtomicResult JSON", json)?;
    let result = chematic_mol::parse_atomic_result_with_limits(json, &wasm_qcschema_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_mol::write_atomic_result(&result))
}

// ===========================================================================
// ORCA input (.inp)
// ===========================================================================

fn orca_atom_to_json(a: &chematic_mol::OrcaAtom) -> serde_json::Value {
    serde_json::json!({
        "element": a.element.symbol(),
        "x": a.x, "y": a.y, "z": a.z,
        "frozen": a.frozen,
        "extra": a.extra,
    })
}

fn orca_atom_from_json(v: &serde_json::Value) -> Result<chematic_mol::OrcaAtom, String> {
    let frozen: [bool; 3] = match v.get("frozen") {
        Some(x) => {
            serde_json::from_value(x.clone()).map_err(|e| format!("invalid 'frozen': {e}"))?
        }
        None => [false, false, false],
    };
    Ok(chematic_mol::OrcaAtom {
        element: element_from_symbol(&jstr(v, "element")?)?,
        x: jf64(v, "x")?,
        y: jf64(v, "y")?,
        z: jf64(v, "z")?,
        frozen,
        extra: jopt_str(v, "extra"),
    })
}

fn orca_coords_to_json(c: &chematic_mol::OrcaCoords) -> serde_json::Value {
    match c {
        chematic_mol::OrcaCoords::Xyz {
            charge,
            multiplicity,
            atoms,
        } => serde_json::json!({
            "type": "xyz",
            "charge": charge,
            "multiplicity": multiplicity,
            "atoms": atoms.iter().map(orca_atom_to_json).collect::<Vec<_>>(),
        }),
        chematic_mol::OrcaCoords::XyzFile {
            charge,
            multiplicity,
            filename,
        } => serde_json::json!({
            "type": "xyzfile",
            "charge": charge,
            "multiplicity": multiplicity,
            "filename": filename,
        }),
        chematic_mol::OrcaCoords::GzmtFile {
            charge,
            multiplicity,
            filename,
        } => serde_json::json!({
            "type": "gzmtfile",
            "charge": charge,
            "multiplicity": multiplicity,
            "filename": filename,
        }),
        chematic_mol::OrcaCoords::Internal {
            charge,
            multiplicity,
            raw,
        } => serde_json::json!({
            "type": "internal",
            "charge": charge,
            "multiplicity": multiplicity,
            "raw": raw,
        }),
    }
}

fn orca_coords_from_json(v: &serde_json::Value) -> Result<chematic_mol::OrcaCoords, String> {
    let charge = ji64(v, "charge")? as i32;
    let multiplicity = ji64(v, "multiplicity")? as u32;
    match jstr(v, "type")?.as_str() {
        "xyz" => {
            let atoms = jarr(v, "atoms")?
                .iter()
                .map(orca_atom_from_json)
                .collect::<Result<Vec<_>, _>>()?;
            Ok(chematic_mol::OrcaCoords::Xyz {
                charge,
                multiplicity,
                atoms,
            })
        }
        "xyzfile" => Ok(chematic_mol::OrcaCoords::XyzFile {
            charge,
            multiplicity,
            filename: jstr(v, "filename")?,
        }),
        "gzmtfile" => Ok(chematic_mol::OrcaCoords::GzmtFile {
            charge,
            multiplicity,
            filename: jstr(v, "filename")?,
        }),
        "internal" => Ok(chematic_mol::OrcaCoords::Internal {
            charge,
            multiplicity,
            raw: jopt_str(v, "raw").unwrap_or_default(),
        }),
        other => Err(format!(
            "unknown coords 'type' '{other}' (expected 'xyz'/'xyzfile'/'gzmtfile'/'internal')"
        )),
    }
}

fn orca_block_to_json(b: &chematic_mol::OrcaBlock) -> serde_json::Value {
    serde_json::json!({ "name": b.name, "raw": b.raw, "has_end": b.has_end })
}

fn orca_block_from_json(v: &serde_json::Value) -> Result<chematic_mol::OrcaBlock, String> {
    Ok(chematic_mol::OrcaBlock {
        name: jstr(v, "name")?,
        raw: jstr(v, "raw")?,
        has_end: v.get("has_end").and_then(|x| x.as_bool()).unwrap_or(true),
    })
}

/// Parse an ORCA input file (`.inp`) and return a `MolHandle` (topology
/// only -- element list, no bonds; ORCA input carries no bond table).
/// Returns a JS error unless the file's coordinate block is an embedded
/// `* xyz ... *` block -- `xyzfile`/`gzmtfile`/`int` (Z-matrix) blocks
/// carry no atom list to convert, or none is present at all. Use
/// [`orca_input_coords_json`] to recover coordinates + charge +
/// multiplicity in the same atom order, or [`orca_input_to_json`] for the
/// full input (comments/keywords/blocks/any coordinate-block kind).
#[wasm_bindgen]
pub fn mol_from_orca_input(text: &str) -> Result<MolHandle, JsValue> {
    check_input_len("ORCA input", text)?;
    let input = chematic_mol::parse_orca_input_with_limits(text, &wasm_orca_input_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let coords = input
        .coords
        .as_ref()
        .and_then(|c| c.to_molecule())
        .ok_or_else(|| {
            JsValue::from_str(
                "ORCA input has no embedded '* xyz ... *' coordinate block (xyzfile/gzmtfile/int/absent)",
            )
        })?;
    let (mol, _coords, _charge, _mult) = coords;
    mol_handle_from_molecule(mol)
}

/// Coordinates + charge + multiplicity from an ORCA input file's embedded
/// `* xyz ... *` block, in the SAME atom order [`mol_from_orca_input`]
/// returns topology for. Returns JSON
/// `{"coords":[[x,y,z],...],"charge":0,"multiplicity":1}`, or a JS error
/// under the same conditions as [`mol_from_orca_input`].
/// This legacy view returns only Angstrom input. Bohr, conflicting units and
/// %coords overrides are refused; use [`orca_cartesian_geometry_json`] instead.
#[wasm_bindgen]
pub fn orca_input_coords_json(text: &str) -> Result<String, JsValue> {
    check_input_len("ORCA input", text)?;
    let input = chematic_mol::parse_orca_cartesian_with_limits(text, &wasm_orca_input_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if input.coordinate_units != chematic_mol::OrcaCoordinateUnits::Angstrom {
        return Err(JsValue::from_str(
            "Bohr coordinates require orca_cartesian_geometry_json with explicit units",
        ));
    }
    let out = serde_json::json!({
        "coords": input.atoms.iter().map(|a| [a.x,a.y,a.z]).collect::<Vec<_>>(),
        "charge": input.charge,
        "multiplicity": input.multiplicity,
    });
    to_json_string(&out)
}

/// Read ordered Cartesian values and explicit Angstrom/Bohr units; no rescaling.
/// Refuses ambiguous/unsupported unit or coordinate declarations and %coords.
/// Atom extras remain opaque; no external geometry resources are opened.
#[wasm_bindgen]
pub fn orca_cartesian_geometry_json(text: &str) -> Result<String, JsValue> {
    let input = chematic_mol::parse_orca_cartesian_with_limits(text, &wasm_orca_input_limits())
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    to_json_string(&input.to_json())
}

/// Parse an ORCA input file and return every field as JSON:
/// `{"comments":[...],"keywords":[...],
/// "blocks":[{"name":"scf","raw":"...","has_end":true},...],
/// "coords":{"type":"xyz"|"xyzfile"|"gzmtfile"|"internal",...}|null}`.
/// See [`chematic_mol::OrcaInput`]'s doc comment for each field's meaning.
#[wasm_bindgen]
pub fn orca_input_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("ORCA input", text)?;
    let input = chematic_mol::parse_orca_input_with_limits(text, &wasm_orca_input_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let out = serde_json::json!({
        "comments": input.comments,
        "keywords": input.keywords,
        "blocks": input.blocks.iter().map(orca_block_to_json).collect::<Vec<_>>(),
        "coords": input.coords.as_ref().map(orca_coords_to_json),
    });
    to_json_string(&out)
}

/// Write an ORCA input file from the JSON shape [`orca_input_to_json`]
/// returns.
#[wasm_bindgen]
pub fn write_orca_input_json(json: &str) -> Result<String, JsValue> {
    let v = parse_json_value("ORCA input JSON", json)?;
    let comments: Vec<String> = match v.get("comments") {
        Some(x) => {
            bounded_array(x, "comments").map_err(|e| JsValue::from_str(&e))?;
            serde_json::from_value(x.clone())
                .map_err(|e| JsValue::from_str(&format!("invalid 'comments': {e}")))?
        }
        None => Vec::new(),
    };
    let keywords: Vec<String> = match v.get("keywords") {
        Some(x) => {
            bounded_array(x, "keywords").map_err(|e| JsValue::from_str(&e))?;
            serde_json::from_value(x.clone())
                .map_err(|e| JsValue::from_str(&format!("invalid 'keywords': {e}")))?
        }
        None => Vec::new(),
    };
    let blocks = match v.get("blocks") {
        Some(x @ serde_json::Value::Array(_)) => bounded_array(x, "blocks")
            .map_err(|e| JsValue::from_str(&e))?
            .iter()
            .map(orca_block_from_json)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| JsValue::from_str(&e))?,
        _ => Vec::new(),
    };
    let coords = match v.get("coords") {
        None | Some(serde_json::Value::Null) => None,
        Some(c) => Some(orca_coords_from_json(c).map_err(|e| JsValue::from_str(&e))?),
    };
    let input = chematic_mol::OrcaInput {
        comments,
        keywords,
        blocks,
        coords,
    };
    Ok(chematic_mol::write_orca_input(&input))
}

// ===========================================================================
// ORCA output (.out / .log) -- read-only, no writer (it's a job log)
// ===========================================================================

fn orca_geometry_frame_to_json(f: &chematic_mol::GeometryFrame) -> serde_json::Value {
    let elements: Vec<&str> = f
        .mol
        .atoms()
        .map(|(_, atom)| atom.element.symbol())
        .collect();
    serde_json::json!({
        "elements": elements,
        "coords": coords_tuples_to_json(&f.coords),
    })
}

fn orca_termination_to_json(t: &chematic_mol::OrcaTermination) -> serde_json::Value {
    match t {
        chematic_mol::OrcaTermination::Normal => serde_json::json!({"kind": "normal"}),
        chematic_mol::OrcaTermination::Error(detail) => {
            serde_json::json!({"kind": "error", "detail": detail})
        }
        chematic_mol::OrcaTermination::Incomplete => serde_json::json!({"kind": "incomplete"}),
    }
}

fn orca_opt_convergence_str(c: chematic_mol::OrcaOptConvergence) -> &'static str {
    match c {
        chematic_mol::OrcaOptConvergence::NotRequested => "not_requested",
        chematic_mol::OrcaOptConvergence::Converged => "converged",
        chematic_mol::OrcaOptConvergence::NotConverged => "not_converged",
        chematic_mol::OrcaOptConvergence::Unknown => "unknown",
    }
}

/// Parse an ORCA output file (`.out`/`.log`) and return every extracted
/// field as JSON: `{"charge":N|null,"multiplicity":N|null,
/// "final_energy_hartree":N|null,
/// "trajectory":[{"elements":[...],"coords":[[x,y,z],...]},...],
/// "frequencies_cm1":[...],
/// "termination":{"kind":"normal"|"error"|"incomplete","detail":"..."?},
/// "optimization_convergence":"not_requested"|"converged"|"not_converged"|"unknown"}`.
/// No writer is provided -- an ORCA output file is a job log, not a
/// document this crate constructs.
#[wasm_bindgen]
pub fn orca_output_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("ORCA output", text)?;
    let orca_limits = chematic_mol::OrcaOutputParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_line_bytes: WASM_MAX_INPUT_BYTES,
        max_geometry_frames: WASM_MAX_BATCH_ITEMS,
        max_geometry_atoms: WASM_MAX_ATOMS,
        ..Default::default()
    };
    let output = chematic_mol::parse_orca_output_with_limits(text, &orca_limits)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let out = serde_json::json!({
        "charge": output.charge,
        "multiplicity": output.multiplicity,
        "final_energy_hartree": output.final_energy_hartree,
        "trajectory": output.trajectory.iter().map(orca_geometry_frame_to_json).collect::<Vec<_>>(),
        "frequencies_cm1": output.frequencies_cm1,
        "termination": orca_termination_to_json(&output.termination),
        "optimization_convergence": orca_opt_convergence_str(output.optimization_convergence),
    });
    to_json_string(&out)
}

#[cfg(test)]
use lammps::{
    cube_shape_vec, cube_values_vec, dump_cartesian_flat, dump_rows_flat, opendx_shape_vec,
    opendx_values_vec,
};

#[cfg(test)]
include!("format_io_tests.rs");
