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
//! ## Grid JSON is a full round trip, not a typed-array view
//!
//! [`chematic_mol::VolumetricGrid::values`] can be very large (Cube/OpenDX
//! routinely carry hundreds of thousands to millions of voxels). This crate
//! has no existing `js_sys::Float64Array`/typed-array precedent to follow
//! (checked: none of `mol_3d.rs`/`mol_depict.rs`/`mol_descriptors.rs`/
//! `mol_edit.rs`/`mol_fingerprints.rs`/`mol_io.rs`/`mol_reactions.rs` use
//! one), so `cube_grid_json`/`opendx_grid_json` serialize `values` as a
//! plain JSON number array -- a real, disclosed perf cost (a full parse +
//! JSON-string materialization + JS-side `JSON.parse` of every voxel value)
//! versus a typed-array accessor view that would avoid it. Acceptable for a
//! first pass; revisit with a `js_sys::Float64Array`-returning accessor
//! struct if a large-grid workload proves this too slow in practice.

use crate::{MolHandle, WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES};
use wasm_bindgen::prelude::*;

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
    serde_json::to_string(v).map_err(|e| JsValue::from_str(&e.to_string()))
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
    jget(v, key)?
        .as_array()
        .ok_or_else(|| format!("'{key}' must be an array"))
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
    let grid = chematic_mol::parse_cube(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let grid = chematic_mol::parse_cube(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let grid = chematic_mol::parse_opendx(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let result = chematic_mol::parse_mmcif(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (mol, _coords) = result.to_molecule();
    mol_handle_from_molecule(mol)
}

/// Cartesian coordinates from an mmCIF file, in the SAME atom order
/// [`mol_from_mmcif`] returns topology for. Returns JSON `[[x,y,z],...]`
/// (Å).
#[wasm_bindgen]
pub fn mmcif_coords_json(text: &str) -> Result<String, JsValue> {
    check_input_len("mmCIF input", text)?;
    let result = chematic_mol::parse_mmcif(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let result = chematic_mol::parse_mmcif(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let out = serde_json::json!({
        "atoms": result.atoms.iter().map(mmcif_atom_to_json).collect::<Vec<_>>(),
        "cell": result.cell.as_ref().map(unit_cell_to_json),
        "space_group": result.space_group,
        "unhandled_columns": result.unhandled_columns,
    });
    to_json_string(&out)
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
    let atoms: Vec<chematic_mol::MmcifAtomRecord> = records
        .as_array()
        .ok_or_else(|| JsValue::from_str("records_json must be a JSON array"))?
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
    let result = chematic_mol::parse_pqr(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (mol, _coords) = result.to_molecule();
    mol_handle_from_molecule(mol)
}

/// Cartesian coordinates from a PQR file, in the SAME atom order
/// [`mol_from_pqr`] returns topology for. Returns JSON `[[x,y,z],...]` (Å).
#[wasm_bindgen]
pub fn pqr_coords_json(text: &str) -> Result<String, JsValue> {
    check_input_len("PQR input", text)?;
    let result = chematic_mol::parse_pqr(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (_mol, coords) = result.to_molecule();
    to_json_string(&coords_tuples_to_json(&coords))
}

/// Parse a PQR file and return every field (charge, radius, chain,
/// residue, inferred element, ...) as JSON: `{"atoms":[{...}]}`. See
/// [`chematic_mol::PqrAtomRecord`]'s doc comment for each field's meaning.
#[wasm_bindgen]
pub fn pqr_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("PQR input", text)?;
    let result = chematic_mol::parse_pqr(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let atoms: Vec<chematic_mol::PqrAtomRecord> = records
        .as_array()
        .ok_or_else(|| JsValue::from_str("records_json must be a JSON array"))?
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
    let qc = chematic_mol::parse_qcschema_molecule(json)
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
    let qc = chematic_mol::parse_qcschema_molecule(json)
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
    let input =
        chematic_mol::parse_atomic_input(json).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(chematic_mol::write_atomic_input(&input))
}

/// Like [`qcschema_validate_atomic_input`], for a QCSchema
/// `qcschema_output`/`qc_schema_output` (`AtomicResult`) document.
#[wasm_bindgen]
pub fn qcschema_validate_atomic_result(json: &str) -> Result<String, JsValue> {
    check_input_len("QCSchema AtomicResult JSON", json)?;
    let result =
        chematic_mol::parse_atomic_result(json).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let input =
        chematic_mol::parse_orca_input(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
#[wasm_bindgen]
pub fn orca_input_coords_json(text: &str) -> Result<String, JsValue> {
    check_input_len("ORCA input", text)?;
    let input =
        chematic_mol::parse_orca_input(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let (_mol, coords, charge, multiplicity) = input
        .coords
        .as_ref()
        .and_then(|c| c.to_molecule())
        .ok_or_else(|| {
            JsValue::from_str(
                "ORCA input has no embedded '* xyz ... *' coordinate block (xyzfile/gzmtfile/int/absent)",
            )
        })?;
    let out = serde_json::json!({
        "coords": coords_tuples_to_json(&coords),
        "charge": charge,
        "multiplicity": multiplicity,
    });
    to_json_string(&out)
}

/// Parse an ORCA input file and return every field as JSON:
/// `{"comments":[...],"keywords":[...],
/// "blocks":[{"name":"scf","raw":"...","has_end":true},...],
/// "coords":{"type":"xyz"|"xyzfile"|"gzmtfile"|"internal",...}|null}`.
/// See [`chematic_mol::OrcaInput`]'s doc comment for each field's meaning.
#[wasm_bindgen]
pub fn orca_input_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("ORCA input", text)?;
    let input =
        chematic_mol::parse_orca_input(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
        Some(x) => serde_json::from_value(x.clone())
            .map_err(|e| JsValue::from_str(&format!("invalid 'comments': {e}")))?,
        None => Vec::new(),
    };
    let keywords: Vec<String> = match v.get("keywords") {
        Some(x) => serde_json::from_value(x.clone())
            .map_err(|e| JsValue::from_str(&format!("invalid 'keywords': {e}")))?,
        None => Vec::new(),
    };
    let blocks = match v.get("blocks") {
        Some(serde_json::Value::Array(arr)) => arr
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
    let output =
        chematic_mol::parse_orca_output(text).map_err(|e| JsValue::from_str(&e.to_string()))?;
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

// ===========================================================================
// LAMMPS data (`read_data` format) + dump/trajectory
// ===========================================================================

fn lammps_box_to_json(b: &chematic_mol::LammpsBox) -> serde_json::Value {
    serde_json::json!({ "lo": b.lo, "hi": b.hi, "tilt": b.tilt })
}

fn lammps_box_from_json(v: &serde_json::Value) -> Result<chematic_mol::LammpsBox, String> {
    let lo: [f64; 3] =
        serde_json::from_value(jget(v, "lo")?.clone()).map_err(|e| format!("invalid 'lo': {e}"))?;
    let hi: [f64; 3] =
        serde_json::from_value(jget(v, "hi")?.clone()).map_err(|e| format!("invalid 'hi': {e}"))?;
    let tilt: Option<[f64; 3]> = match v.get("tilt") {
        None | Some(serde_json::Value::Null) => None,
        Some(x) => {
            Some(serde_json::from_value(x.clone()).map_err(|e| format!("invalid 'tilt': {e}"))?)
        }
    };
    Ok(chematic_mol::LammpsBox { lo, hi, tilt })
}

/// Maps a JS-facing `atom_style` string onto [`chematic_mol::LammpsAtomStyle`].
/// Any value other than `"atomic"`/`"charge"`/`"molecular"`/`"full"` becomes
/// [`chematic_mol::LammpsAtomStyle::Other`], which `parse_lammps_data`
/// itself rejects with a typed `UnsupportedAtomStyle` error -- this
/// function does not pre-validate, matching the Rust API's own fail-closed
/// behavior rather than duplicating it.
fn lammps_atom_style_from_str(s: &str) -> chematic_mol::LammpsAtomStyle {
    match s {
        "atomic" => chematic_mol::LammpsAtomStyle::Atomic,
        "charge" => chematic_mol::LammpsAtomStyle::Charge,
        "molecular" => chematic_mol::LammpsAtomStyle::Molecular,
        "full" => chematic_mol::LammpsAtomStyle::Full,
        other => chematic_mol::LammpsAtomStyle::Other(other.to_string()),
    }
}

fn lammps_atom_style_str(s: &chematic_mol::LammpsAtomStyle) -> String {
    match s {
        chematic_mol::LammpsAtomStyle::Atomic => "atomic".to_string(),
        chematic_mol::LammpsAtomStyle::Charge => "charge".to_string(),
        chematic_mol::LammpsAtomStyle::Molecular => "molecular".to_string(),
        chematic_mol::LammpsAtomStyle::Full => "full".to_string(),
        chematic_mol::LammpsAtomStyle::Other(s) => s.clone(),
    }
}

fn lammps_mass_to_json(m: &chematic_mol::LammpsMass) -> serde_json::Value {
    serde_json::json!({ "atom_type": m.atom_type, "mass": m.mass })
}
fn lammps_mass_from_json(v: &serde_json::Value) -> Result<chematic_mol::LammpsMass, String> {
    Ok(chematic_mol::LammpsMass {
        atom_type: ji64(v, "atom_type")?,
        mass: jf64(v, "mass")?,
    })
}

fn lammps_atom_to_json(a: &chematic_mol::LammpsAtom) -> serde_json::Value {
    serde_json::json!({
        "id": a.id,
        "molecule_id": a.molecule_id,
        "atom_type": a.atom_type,
        "charge": a.charge,
        "x": a.x, "y": a.y, "z": a.z,
        "image": a.image,
    })
}
fn lammps_atom_from_json(v: &serde_json::Value) -> Result<chematic_mol::LammpsAtom, String> {
    let image: Option<[i32; 3]> = match v.get("image") {
        None | Some(serde_json::Value::Null) => None,
        Some(x) => {
            Some(serde_json::from_value(x.clone()).map_err(|e| format!("invalid 'image': {e}"))?)
        }
    };
    Ok(chematic_mol::LammpsAtom {
        id: ji64(v, "id")?,
        molecule_id: jopt_i64(v, "molecule_id"),
        atom_type: ji64(v, "atom_type")?,
        charge: v.get("charge").and_then(|x| x.as_f64()),
        x: jf64(v, "x")?,
        y: jf64(v, "y")?,
        z: jf64(v, "z")?,
        image,
    })
}

fn lammps_velocity_to_json(v: &chematic_mol::LammpsVelocity) -> serde_json::Value {
    serde_json::json!({ "atom_id": v.atom_id, "vx": v.vx, "vy": v.vy, "vz": v.vz })
}
fn lammps_velocity_from_json(
    v: &serde_json::Value,
) -> Result<chematic_mol::LammpsVelocity, String> {
    Ok(chematic_mol::LammpsVelocity {
        atom_id: ji64(v, "atom_id")?,
        vx: jf64(v, "vx")?,
        vy: jf64(v, "vy")?,
        vz: jf64(v, "vz")?,
    })
}

fn lammps_bond_to_json(b: &chematic_mol::LammpsBond) -> serde_json::Value {
    serde_json::json!({ "id": b.id, "bond_type": b.bond_type, "atom1": b.atom1, "atom2": b.atom2 })
}
fn lammps_bond_from_json(v: &serde_json::Value) -> Result<chematic_mol::LammpsBond, String> {
    Ok(chematic_mol::LammpsBond {
        id: ji64(v, "id")?,
        bond_type: ji64(v, "bond_type")?,
        atom1: ji64(v, "atom1")?,
        atom2: ji64(v, "atom2")?,
    })
}

/// Parse a LAMMPS data file (`read_data` format) and return every section
/// as JSON: `{"counts":[["atoms",120],["atom types",4],...],
/// "atom_style":"atomic"|"charge"|"molecular"|"full"|"<other>",
/// "simulation_box":{"lo":[x,y,z],"hi":[x,y,z],"tilt":[xy,xz,yz]|null},
/// "masses":[{"atom_type":N,"mass":N}],
/// "atoms":[{"id":N,"molecule_id":N|null,"atom_type":N,"charge":N|null,"x":N,"y":N,"z":N,"image":[ix,iy,iz]|null}],
/// "velocities":[{"atom_id":N,"vx":N,"vy":N,"vz":N}],
/// "bonds":[{"id":N,"bond_type":N,"atom1":N,"atom2":N}],
/// "unparsed_sections":[["Angles","<raw row text>"],...]}`. `atom_type`
/// must be exactly `"atomic"`/`"charge"`/`"molecular"`/`"full"` -- LAMMPS's
/// atom style is not recoverable from the file itself (see
/// [`chematic_mol::LammpsData`]'s module doc comment); any other value is
/// rejected with a JS error, matching
/// [`chematic_mol::LammpsDataError::UnsupportedAtomStyle`].
///
/// This module has no bond-perception step of its own: `Angles`/
/// `Dihedrals`/`Impropers`/`*Coeffs`/any other section not listed above
/// are preserved verbatim (byte-for-byte, `#` comments included) in
/// `unparsed_sections`, not modeled field-by-field.
#[wasm_bindgen]
pub fn lammps_data_to_json(text: &str, atom_style: &str) -> Result<String, JsValue> {
    check_input_len("LAMMPS data input", text)?;
    let data = chematic_mol::parse_lammps_data(text, lammps_atom_style_from_str(atom_style))
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let counts: Vec<serde_json::Value> = data
        .counts
        .iter()
        .map(|(k, v)| serde_json::json!([k, v]))
        .collect();
    let unparsed: Vec<serde_json::Value> = data
        .unparsed_sections
        .iter()
        .map(|(k, v)| serde_json::json!([k, v]))
        .collect();
    let out = serde_json::json!({
        "counts": counts,
        "atom_style": lammps_atom_style_str(&data.atom_style),
        "simulation_box": lammps_box_to_json(&data.simulation_box),
        "masses": data.masses.iter().map(lammps_mass_to_json).collect::<Vec<_>>(),
        "atoms": data.atoms.iter().map(lammps_atom_to_json).collect::<Vec<_>>(),
        "velocities": data.velocities.iter().map(lammps_velocity_to_json).collect::<Vec<_>>(),
        "bonds": data.bonds.iter().map(lammps_bond_to_json).collect::<Vec<_>>(),
        "unparsed_sections": unparsed,
    });
    to_json_string(&out)
}

/// Write a LAMMPS data file from the JSON shape [`lammps_data_to_json`]
/// returns.
#[wasm_bindgen]
pub fn write_lammps_data_json(json: &str) -> Result<String, JsValue> {
    let v = parse_json_value("LAMMPS data JSON", json)?;
    let counts: Vec<(String, i64)> = jarr(&v, "counts")
        .map_err(|e| JsValue::from_str(&e))?
        .iter()
        .map(|pair| {
            let arr = pair
                .as_array()
                .ok_or_else(|| "each 'counts' entry must be a [label, count] pair".to_string())?;
            let label = arr
                .first()
                .and_then(|x| x.as_str())
                .ok_or_else(|| "counts[i][0] must be a string label".to_string())?
                .to_string();
            let count = arr
                .get(1)
                .and_then(|x| x.as_i64())
                .ok_or_else(|| "counts[i][1] must be an integer".to_string())?;
            Ok((label, count))
        })
        .collect::<Result<Vec<_>, String>>()
        .map_err(|e: String| JsValue::from_str(&e))?;
    let unparsed_sections: Vec<(String, String)> = match v.get("unparsed_sections") {
        Some(serde_json::Value::Array(arr)) => arr
            .iter()
            .map(|pair| {
                let a = pair.as_array().ok_or_else(|| {
                    "each 'unparsed_sections' entry must be a [name, raw] pair".to_string()
                })?;
                let name = a
                    .first()
                    .and_then(|x| x.as_str())
                    .ok_or_else(|| "unparsed_sections[i][0] must be a string".to_string())?
                    .to_string();
                let raw = a
                    .get(1)
                    .and_then(|x| x.as_str())
                    .ok_or_else(|| "unparsed_sections[i][1] must be a string".to_string())?
                    .to_string();
                Ok((name, raw))
            })
            .collect::<Result<Vec<_>, String>>()
            .map_err(|e: String| JsValue::from_str(&e))?,
        _ => Vec::new(),
    };
    let masses = match v.get("masses") {
        Some(serde_json::Value::Array(arr)) => arr
            .iter()
            .map(lammps_mass_from_json)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| JsValue::from_str(&e))?,
        _ => Vec::new(),
    };
    let atoms = match v.get("atoms") {
        Some(serde_json::Value::Array(arr)) => arr
            .iter()
            .map(lammps_atom_from_json)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| JsValue::from_str(&e))?,
        _ => Vec::new(),
    };
    let velocities = match v.get("velocities") {
        Some(serde_json::Value::Array(arr)) => arr
            .iter()
            .map(lammps_velocity_from_json)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| JsValue::from_str(&e))?,
        _ => Vec::new(),
    };
    let bonds = match v.get("bonds") {
        Some(serde_json::Value::Array(arr)) => arr
            .iter()
            .map(lammps_bond_from_json)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| JsValue::from_str(&e))?,
        _ => Vec::new(),
    };
    let simulation_box =
        lammps_box_from_json(jget(&v, "simulation_box").map_err(|e| JsValue::from_str(&e))?)
            .map_err(|e| JsValue::from_str(&e))?;
    let atom_style =
        lammps_atom_style_from_str(&jstr(&v, "atom_style").map_err(|e| JsValue::from_str(&e))?);

    let data = chematic_mol::LammpsData {
        counts,
        atom_style,
        simulation_box,
        masses,
        atoms,
        velocities,
        bonds,
        unparsed_sections,
    };
    Ok(chematic_mol::write_lammps_data(&data))
}

fn lammps_dump_frame_to_json(f: &chematic_mol::LammpsDumpFrame) -> serde_json::Value {
    serde_json::json!({
        "timestep": f.timestep,
        "num_atoms": f.num_atoms,
        "box_bounds": lammps_box_to_json(&f.box_bounds),
        "boundary_flags": f.boundary_flags,
        "column_names": f.column_names,
        "rows": f.rows,
    })
}

fn lammps_dump_frame_from_json(
    v: &serde_json::Value,
) -> Result<chematic_mol::LammpsDumpFrame, String> {
    let boundary_flags: [String; 3] = serde_json::from_value(jget(v, "boundary_flags")?.clone())
        .map_err(|e| format!("invalid 'boundary_flags': {e}"))?;
    let column_names: Vec<String> = serde_json::from_value(jget(v, "column_names")?.clone())
        .map_err(|e| format!("invalid 'column_names': {e}"))?;
    let rows: Vec<Vec<f64>> = serde_json::from_value(jget(v, "rows")?.clone())
        .map_err(|e| format!("invalid 'rows': {e}"))?;
    Ok(chematic_mol::LammpsDumpFrame {
        timestep: ji64(v, "timestep")?,
        num_atoms: ji64_or(v, "num_atoms", rows.len() as i64) as usize,
        box_bounds: lammps_box_from_json(jget(v, "box_bounds")?)?,
        boundary_flags,
        column_names,
        rows,
    })
}

/// Parse a single LAMMPS dump/trajectory frame and return it as JSON:
/// `{"timestep":N,"num_atoms":N,
/// "box_bounds":{"lo":[x,y,z],"hi":[x,y,z],"tilt":[xy,xz,yz]|null},
/// "boundary_flags":["pp","pp","pp"],"column_names":[...],
/// "rows":[[...values, one per column_names entry...],...]}`.
/// `box_bounds` is the dump file's raw bounding-box shape -- see
/// [`chematic_mol::box_bounds_to_true`] if the true (un-shifted) box is
/// needed for a triclinic frame.
#[wasm_bindgen]
pub fn lammps_dump_frame_to_json_str(text: &str) -> Result<String, JsValue> {
    check_input_len("LAMMPS dump input", text)?;
    let frame = chematic_mol::parse_lammps_dump_frame(text)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    to_json_string(&lammps_dump_frame_to_json(&frame))
}

/// Parse every frame of a LAMMPS dump/trajectory file and return them as a
/// JSON array (same per-frame shape as [`lammps_dump_frame_to_json_str`]).
///
/// This reads the whole input, parses it fully, and returns every frame at
/// once -- [`chematic_mol::LammpsDumpReader`]'s per-frame streaming
/// iteration (reading one frame at a time from a `BufRead` without holding
/// the whole trajectory in memory) has no natural equivalent across the
/// JS/WASM boundary in this first pass and is deliberately not exposed
/// here, not silently dropped: a JS caller with a truly large trajectory
/// that needs bounded memory should process it server-side instead.
#[wasm_bindgen]
pub fn lammps_trajectory_to_json(text: &str) -> Result<String, JsValue> {
    check_input_len("LAMMPS dump input", text)?;
    let cursor = std::io::Cursor::new(text.as_bytes());
    let reader = std::io::BufReader::new(cursor);
    let mut frames = Vec::new();
    for frame in chematic_mol::LammpsDumpReader::new(reader) {
        let frame = frame.map_err(|e| JsValue::from_str(&e.to_string()))?;
        frames.push(lammps_dump_frame_to_json(&frame));
    }
    to_json_string(&serde_json::Value::Array(frames))
}

/// Write a single LAMMPS dump frame from the JSON shape
/// [`lammps_dump_frame_to_json_str`] returns.
#[wasm_bindgen]
pub fn write_lammps_dump_frame_json(json: &str) -> Result<String, JsValue> {
    let v = parse_json_value("LAMMPS dump frame JSON", json)?;
    let frame = lammps_dump_frame_from_json(&v).map_err(|e| JsValue::from_str(&e))?;
    Ok(chematic_mol::write_lammps_dump_frame(&frame))
}

/// Write a LAMMPS trajectory (N frames concatenated back to back, matching
/// [`chematic_mol::write_lammps_trajectory`]) from a JSON array of frames
/// in the shape [`lammps_dump_frame_to_json_str`] returns.
#[wasm_bindgen]
pub fn write_lammps_trajectory_json(json: &str) -> Result<String, JsValue> {
    let v = parse_json_value("LAMMPS trajectory JSON", json)?;
    let frames: Vec<chematic_mol::LammpsDumpFrame> = v
        .as_array()
        .ok_or_else(|| JsValue::from_str("trajectory JSON must be an array of frames"))?
        .iter()
        .map(lammps_dump_frame_from_json)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| JsValue::from_str(&e))?;
    Ok(chematic_mol::write_lammps_trajectory(&frames))
}

// ===========================================================================
// Tests
//
// Only happy paths are exercised here: constructing a `JsValue` (as every
// `Err` branch above does) outside a real wasm/JS runtime aborts the native
// test process rather than returning an `Err` -- the same constraint
// documented on `crate::mol_io::extxyz_frame_from_json_args` and the
// `run_reactants` error-path note in `crate::tests`. Error-path coverage of
// each hand-written JSON<->struct mapping is implicit: every `Ok` path
// below round-trips through the exact same field-extraction code an error
// would come from.
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const CUBE_2X2X2: &str = "Water density\n\
Generated for chematic tests\n\
1    0.000000    0.000000    0.000000\n\
2    1.000000    0.000000    0.000000\n\
2    0.000000    1.000000    0.000000\n\
2    0.000000    0.000000    1.000000\n\
8    8.000000    0.500000    0.500000    0.500000\n\
0.0 1.0 2.0 3.0\n\
4.0 5.0 6.0 7.0\n";

    const OPENDX_2X2X2: &str = "object 1 class gridpositions counts 2 2 2\n\
origin -1.0 -1.0 -1.0\n\
delta 0.5 0.0 0.0\n\
delta 0.0 0.5 0.0\n\
delta 0.0 0.0 0.5\n\
object 2 class gridconnections counts 2 2 2\n\
object 3 class array type double rank 0 items 8 data follows\n\
0.0 1.0 2.0\n\
3.0 4.0 5.0\n\
6.0 7.0\n\
attribute \"dep\" string \"positions\"\n\
object \"regular positions regular connections\" class field\n\
component \"positions\" value 1\n\
component \"connections\" value 2\n\
component \"data\" value 3\n";

    #[test]
    fn cube_round_trip_via_json() {
        assert_eq!(mol_from_cube(CUBE_2X2X2).unwrap().atom_count(), 1);
        let json1 = cube_grid_json(CUBE_2X2X2).unwrap();
        let text2 = write_cube_json(&json1).unwrap();
        let json2 = cube_grid_json(&text2).unwrap();
        assert_eq!(json1, json2);
    }

    #[test]
    fn opendx_round_trip_via_json() {
        let json1 = opendx_grid_json(OPENDX_2X2X2).unwrap();
        // parse_opendx always yields Angstrom (no unit tag in the format),
        // so the fail-closed default writer must accept it directly.
        let text2 = write_opendx_json(&json1).unwrap();
        let json2 = opendx_grid_json(&text2).unwrap();
        assert_eq!(json1, json2);
        // The lossy writer must also accept an already-Angstrom grid.
        let text3 = write_opendx_lossy_json(&json1).unwrap();
        assert_eq!(opendx_grid_json(&text3).unwrap(), json1);
    }

    const MMCIF_ATOMS_JSON: &str = r#"[
        {"group_pdb":"ATOM","serial":1,"element":"C","atom_name":"CA","alt_loc":null,
         "res_name":"ALA","chain_id":"A","res_seq":1,"label_seq_id":1,"icode":null,
         "x":1.0,"y":2.0,"z":3.0,"occupancy":1.0,"b_iso":20.0,"formal_charge":null,
         "entity_id":"1","model_num":1},
        {"group_pdb":"HETATM","serial":2,"element":"Zn","atom_name":"ZN","alt_loc":null,
         "res_name":"ZN","chain_id":"A","res_seq":2,"label_seq_id":null,"icode":"A",
         "x":4.0,"y":5.0,"z":6.0,"occupancy":0.5,"b_iso":15.0,"formal_charge":2,
         "entity_id":null,"model_num":1}
    ]"#;

    #[test]
    fn mmcif_round_trip_via_json() {
        let text = write_mmcif_json(
            MMCIF_ATOMS_JSON,
            r#"{"a":10.0,"b":10.0,"c":10.0,"alpha":90.0,"beta":90.0,"gamma":90.0}"#,
            "P1",
            "TEST",
        )
        .unwrap();
        assert_eq!(mol_from_mmcif(&text).unwrap().atom_count(), 2);
        let coords: Vec<[f64; 3]> =
            serde_json::from_str(&mmcif_coords_json(&text).unwrap()).unwrap();
        assert_eq!(coords, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let out = mmcif_to_json(&text).unwrap();
        assert!(out.contains("\"formal_charge\":2"));
        assert!(out.contains("\"space_group\":\"P1\""));
    }

    const PQR_ATOMS_JSON: &str = r#"[
        {"group_pdb":"ATOM","serial":1,"atom_name":"N","res_name":"ALA","chain_id":null,
         "res_seq":1,"icode":null,"x":1.0,"y":2.0,"z":3.0,"charge":-0.4,"radius":1.5,"element":"N"},
        {"group_pdb":"HETATM","serial":2,"atom_name":"ZN","res_name":"ZN","chain_id":"A",
         "res_seq":2,"icode":null,"x":4.0,"y":5.0,"z":6.0,"charge":2.0,"radius":1.09,"element":"Zn"}
    ]"#;

    #[test]
    fn pqr_round_trip_via_json() {
        let text = write_pqr_json(PQR_ATOMS_JSON).unwrap();
        assert_eq!(mol_from_pqr(&text).unwrap().atom_count(), 2);
        let out = pqr_to_json(&text).unwrap();
        assert!(out.contains("\"radius\":1.09"));
        assert!(out.contains("\"chain_id\":\"A\""));
    }

    #[test]
    fn pqr_infer_element_matches_rust_api() {
        assert_eq!(pqr_infer_element("ATOM", "ALA", "CA").as_deref(), Some("C"));
        assert_eq!(
            pqr_infer_element("HETATM", "ZN", "ZN").as_deref(),
            Some("Zn")
        );
    }

    const QCSCHEMA_WATER_JSON: &str = r#"{
        "schema_name": "qcschema_molecule",
        "schema_version": 1,
        "symbols": ["O", "H", "H"],
        "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.7, -0.5],
        "molecular_charge": 0.0,
        "molecular_multiplicity": 1
    }"#;

    #[test]
    fn qcschema_molecule_round_trip() {
        let handle = mol_from_qcschema_molecule(QCSCHEMA_WATER_JSON).unwrap();
        assert_eq!(handle.atom_count(), 3);
        let coords_json = qcschema_molecule_coords_json(QCSCHEMA_WATER_JSON).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&coords_json).unwrap();
        let coords_only = serde_json::to_string(&parsed["coords"]).unwrap();
        let written = to_qcschema_molecule_json(&handle, &coords_only, 0.0, 1).unwrap();
        let handle2 = mol_from_qcschema_molecule(&written).unwrap();
        assert_eq!(handle2.atom_count(), 3);
    }

    #[test]
    fn qcschema_validate_atomic_input_round_trip() {
        let json = format!(
            r#"{{"schema_name":"qcschema_input","schema_version":1,"molecule":{QCSCHEMA_WATER_JSON},
                "driver":"energy","model":{{"method":"b3lyp","basis":"def2-svp"}},"keywords":{{}}}}"#
        );
        let out = qcschema_validate_atomic_input(&json).unwrap();
        assert!(out.contains("b3lyp"));
    }

    const ORCA_INPUT_JSON: &str = r##"{
        "comments": ["# a comment line"],
        "keywords": ["B3LYP", "def2-SVP"],
        "blocks": [{"name":"pal","raw":"nprocs 4","has_end":true}],
        "coords": {"type":"xyz","charge":0,"multiplicity":1,"atoms":[
            {"element":"O","x":0.0,"y":0.0,"z":0.0,"frozen":[false,false,false],"extra":null},
            {"element":"H","x":0.0,"y":0.0,"z":0.96,"frozen":[false,false,false],"extra":null}
        ]}
    }"##;

    #[test]
    fn orca_input_round_trip_via_json() {
        let text = write_orca_input_json(ORCA_INPUT_JSON).unwrap();
        assert_eq!(mol_from_orca_input(&text).unwrap().atom_count(), 2);
        let out = orca_input_to_json(&text).unwrap();
        assert!(out.contains("\"B3LYP\""));
        assert!(out.contains("\"type\":\"xyz\""));
    }

    #[test]
    fn orca_output_to_json_extracts_energy_and_trajectory() {
        let text = "\
Total Charge           Charge          ....    0
Multiplicity           Mult            ....    1
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     0.000000    0.000000    0.000000
  H     0.000000    0.000000    0.960000
-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -76.320145981234
-------------------------   --------------------
****ORCA TERMINATED NORMALLY****
";
        let out = orca_output_to_json(text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["final_energy_hartree"].as_f64(), Some(-76.320145981234));
        assert_eq!(v["trajectory"].as_array().unwrap().len(), 1);
        assert_eq!(
            v["trajectory"][0]["elements"],
            serde_json::json!(["O", "H"])
        );
        assert_eq!(v["termination"]["kind"], "normal");
    }

    const LAMMPS_DATA_JSON: &str = r#"{
        "counts": [["atoms", 2], ["atom types", 1]],
        "atom_style": "atomic",
        "simulation_box": {"lo":[0.0,0.0,0.0],"hi":[10.0,10.0,10.0],"tilt":null},
        "masses": [{"atom_type":1,"mass":12.011}],
        "atoms": [
            {"id":1,"molecule_id":null,"atom_type":1,"charge":null,"x":1.0,"y":1.0,"z":1.0,"image":null},
            {"id":2,"molecule_id":null,"atom_type":1,"charge":null,"x":2.0,"y":2.0,"z":2.0,"image":[1,0,0]}
        ],
        "velocities": [],
        "bonds": [],
        "unparsed_sections": []
    }"#;

    #[test]
    fn lammps_data_round_trip_via_json() {
        let text = write_lammps_data_json(LAMMPS_DATA_JSON).unwrap();
        let out = lammps_data_to_json(&text, "atomic").unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["atoms"].as_array().unwrap().len(), 2);
        assert_eq!(v["atoms"][1]["image"], serde_json::json!([1, 0, 0]));
        assert_eq!(v["atom_style"], "atomic");
    }

    // Note: an "unsupported atom_style" error-path test is omitted here --
    // like `run_reactants`'s error-path tests in `crate::tests`,
    // constructing the resulting `JsValue` outside a real wasm/JS runtime
    // aborts the native test process. `chematic_mol::LammpsDataError::
    // UnsupportedAtomStyle` coverage lives in chematic-mol's own unit tests;
    // `lammps_atom_style_from_str` (this file's only new logic in that
    // path) is a direct 1:1 string match with no branch of its own to miss.

    const LAMMPS_DUMP_FRAME_JSON: &str = r#"{
        "timestep": 100,
        "num_atoms": 2,
        "box_bounds": {"lo":[0.0,0.0,0.0],"hi":[10.0,10.0,10.0],"tilt":null},
        "boundary_flags": ["pp","pp","pp"],
        "column_names": ["id","x","y","z"],
        "rows": [[1.0,1.0,1.0,1.0],[2.0,2.0,2.0,2.0]]
    }"#;

    #[test]
    fn lammps_dump_frame_round_trip_via_json() {
        let text = write_lammps_dump_frame_json(LAMMPS_DUMP_FRAME_JSON).unwrap();
        let out = lammps_dump_frame_to_json_str(&text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["timestep"], 100);
        assert_eq!(v["rows"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn lammps_trajectory_round_trip_multi_frame() {
        let frames_json = format!("[{LAMMPS_DUMP_FRAME_JSON},{LAMMPS_DUMP_FRAME_JSON}]");
        let text = write_lammps_trajectory_json(&frames_json).unwrap();
        let out = lammps_trajectory_to_json(&text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v.as_array().unwrap().len(), 2);
    }
}
