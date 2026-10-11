use super::*;

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
    let data = chematic_mol::parse_lammps_data_with_limits(
        text,
        lammps_atom_style_from_str(atom_style),
        &wasm_lammps_limits(),
    )
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
    let counts = lammps_counts_from_json(&v).map_err(|e| JsValue::from_str(&e))?;
    let unparsed_sections =
        lammps_unparsed_sections_from_json(&v).map_err(|e| JsValue::from_str(&e))?;
    let masses = lammps_array_from_json(&v, "masses", lammps_mass_from_json)?;
    let atoms = lammps_array_from_json(&v, "atoms", lammps_atom_from_json)?;
    let velocities = lammps_array_from_json(&v, "velocities", lammps_velocity_from_json)?;
    let bonds = lammps_array_from_json(&v, "bonds", lammps_bond_from_json)?;
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

fn lammps_counts_from_json(v: &serde_json::Value) -> Result<Vec<(String, i64)>, String> {
    jarr(v, "counts")?
        .iter()
        .map(|pair| {
            let fields = pair
                .as_array()
                .ok_or_else(|| "each 'counts' entry must be a [label, count] pair".to_string())?;
            let label = fields
                .first()
                .and_then(|value| value.as_str())
                .ok_or_else(|| "counts[i][0] must be a string label".to_string())?;
            let count = fields
                .get(1)
                .and_then(|value| value.as_i64())
                .ok_or_else(|| "counts[i][1] must be an integer".to_string())?;
            Ok((label.to_string(), count))
        })
        .collect()
}

fn lammps_unparsed_sections_from_json(
    v: &serde_json::Value,
) -> Result<Vec<(String, String)>, String> {
    let Some(serde_json::Value::Array(entries)) = v.get("unparsed_sections") else {
        return Ok(Vec::new());
    };
    entries
        .iter()
        .map(|entry| {
            let fields = entry.as_array().ok_or_else(|| {
                "each 'unparsed_sections' entry must be a [name, raw] pair".to_string()
            })?;
            let name = fields
                .first()
                .and_then(|value| value.as_str())
                .ok_or_else(|| "unparsed_sections[i][0] must be a string".to_string())?;
            let raw = fields
                .get(1)
                .and_then(|value| value.as_str())
                .ok_or_else(|| "unparsed_sections[i][1] must be a string".to_string())?;
            Ok((name.to_string(), raw.to_string()))
        })
        .collect()
}

fn lammps_array_from_json<T>(
    v: &serde_json::Value,
    key: &str,
    parse: fn(&serde_json::Value) -> Result<T, String>,
) -> Result<Vec<T>, JsValue> {
    let Some(serde_json::Value::Array(entries)) = v.get(key) else {
        return Ok(Vec::new());
    };
    entries
        .iter()
        .map(parse)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| JsValue::from_str(&error))
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
/// `box_bounds` is already the resolved TRUE simulation box (the parser
/// applies [`chematic_mol::box_bounds_to_true`] internally before
/// `LammpsDumpFrame` is ever built) -- not the file's raw
/// `xlo_bound`/`xhi_bound`/... values. `rows` is the raw per-atom column
/// data as declared by `column_names`, which may be `x y z`
/// (already-Cartesian), `xs ys zs` (box-scaled), `xu yu zu` (unwrapped), or
/// any other dump-command column -- use
/// [`lammps_dump_cartesian_positions_json`] to resolve real Cartesian
/// positions from whichever convention is present, rather than
/// reimplementing that resolution/transform in JS.
#[wasm_bindgen]
pub fn lammps_dump_frame_to_json_str(text: &str) -> Result<String, JsValue> {
    check_input_len("LAMMPS dump input", text)?;
    let frame = chematic_mol::parse_lammps_dump_frame_with_limits(text, &wasm_lammps_dump_limits())
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    to_json_string(&lammps_dump_frame_to_json(&frame))
}

/// Real Cartesian positions for a LAMMPS dump frame (in the JSON shape
/// [`lammps_dump_frame_to_json_str`] returns), resolved by delegating
/// directly to [`chematic_mol::LammpsDumpFrame::cartesian_positions`] --
/// this function does not reimplement any part of the box-bounds or
/// scaled-coordinate math itself; that method is the single place this
/// crate gets the (orthogonal or triclinic) transform right, and every
/// WASM caller must go through it rather than re-deriving the transform in
/// JS (the same reasoning behind this crate's OpenDX fail-closed unit
/// handling and QCSchema's single Bohr<->Ångström conversion point).
///
/// - `x y z` columns: passed straight through.
/// - `xs ys zs` columns: transformed through `frame.box_bounds` (including
///   the triclinic shear terms when a tilt is present).
/// - Neither present (including an `xu yu zu`-only frame -- "unwrapped" is
///   a materially different physical quantity from a scaled coordinate,
///   never resolved by this method): returns JSON `null`, not an error and
///   not an empty array, matching
///   [`chematic_mol::LammpsDumpFrame::cartesian_positions`]'s own
///   `Option` semantics exactly.
///
/// Returns JSON `[[x,y,z],...]` on success, in the same atom order as
/// `frame.rows`. See [`lammps_dump_cartesian_positions_f64`] for a flat
/// `Float64Array` sibling -- note its `null` case becomes an `Err` there
/// instead, a disclosed, real API-shape difference (a typed array has no
/// `null`).
#[wasm_bindgen]
pub fn lammps_dump_cartesian_positions_json(frame_json: &str) -> Result<String, JsValue> {
    let v = parse_json_value("LAMMPS dump frame JSON", frame_json)?;
    let frame = lammps_dump_frame_from_json(&v).map_err(|e| JsValue::from_str(&e))?;
    match frame.cartesian_positions() {
        Some(positions) => {
            serde_json::to_string(&positions).map_err(|e| JsValue::from_str(&e.to_string()))
        }
        None => Ok("null".to_string()),
    }
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
    let limits = chematic_mol::LammpsDumpParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_line_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms_per_frame: WASM_MAX_ATOMS,
        max_columns: 256,
        max_frames: WASM_MAX_BATCH_ITEMS,
    };
    let mut frames = Vec::new();
    for frame in chematic_mol::LammpsDumpReader::with_limits(reader, limits) {
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
    let frames: Vec<chematic_mol::LammpsDumpFrame> = bounded_array(&v, "trajectory JSON")
        .map_err(|e| JsValue::from_str(&e))?
        .iter()
        .map(lammps_dump_frame_from_json)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| JsValue::from_str(&e))?;
    Ok(chematic_mol::write_lammps_trajectory(&frames))
}

// ===========================================================================
// Typed-array accessors (Float64Array/Uint32Array)
//
// Additive alongside the JSON-returning functions above (which stay,
// unchanged) -- see this file's module docs, "Grid JSON is a full round
// trip; typed-array siblings avoid it". Every `#[wasm_bindgen]` function
// below is a thin wrapper around a private helper that returns a plain
// `Vec<f64>`/`[u32; 3]`/`String` error rather than `js_sys`/`JsValue`: a
// `js_sys::Float64Array`/`Uint32Array` is an `extern "C"` JS binding that
// aborts the process when called outside a real JS runtime (same
// constraint `mod tests` below already documents for `JsValue`), so the
// actual parsing/flattening logic is kept out of that unreachable-in-tests
// zone and is exercised natively (against the JSON siblings' own output,
// the strongest cross-check available) in `mod tests`; the Node
// `.test.mjs` suite exercises the real typed arrays plus the
// `lammps_dump_cartesian_positions_f64` error path.
//
// Every helper below delegates to the exact same `chematic_mol` call as
// its JSON-returning sibling -- no parsing/math is reimplemented.
// ===========================================================================

fn shape_to_u32_array(shape: [usize; 3], label: &str) -> Result<[u32; 3], String> {
    let mut out = [0u32; 3];
    for (i, &d) in shape.iter().enumerate() {
        out[i] =
            u32::try_from(d).map_err(|_| format!("{label} shape[{i}]={d} does not fit u32"))?;
    }
    Ok(out)
}

pub(super) fn cube_values_vec(text: &str) -> Result<Vec<f64>, String> {
    chematic_mol::parse_cube_with_limits(text, &wasm_cube_limits())
        .map(|g| g.values)
        .map_err(|e| e.to_string())
}

pub(super) fn cube_shape_vec(text: &str) -> Result<[u32; 3], String> {
    let grid = chematic_mol::parse_cube_with_limits(text, &wasm_cube_limits())
        .map_err(|e| e.to_string())?;
    shape_to_u32_array(grid.shape, "cube")
}

pub(super) fn opendx_values_vec(text: &str) -> Result<Vec<f64>, String> {
    chematic_mol::parse_opendx_with_limits(text, &wasm_opendx_limits())
        .map(|g| g.values)
        .map_err(|e| e.to_string())
}

pub(super) fn opendx_shape_vec(text: &str) -> Result<[u32; 3], String> {
    let grid = chematic_mol::parse_opendx_with_limits(text, &wasm_opendx_limits())
        .map_err(|e| e.to_string())?;
    shape_to_u32_array(grid.shape, "OpenDX")
}

fn dump_frame_from_json_str(frame_json: &str) -> Result<chematic_mol::LammpsDumpFrame, String> {
    let v: serde_json::Value =
        serde_json::from_str(frame_json).map_err(|e| format!("invalid frame JSON: {e}"))?;
    lammps_dump_frame_from_json(&v)
}

/// Flattens `rows` (one row per atom, `column_names.len()` values per row)
/// into a single flat `Vec<f64>`, row-major: atom 0's columns, then atom
/// 1's, etc.
pub(super) fn dump_rows_flat(frame_json: &str) -> Result<Vec<f64>, String> {
    let frame = dump_frame_from_json_str(frame_json)?;
    Ok(frame.rows.into_iter().flatten().collect())
}

/// Same resolution as [`lammps_dump_cartesian_positions_json`]
/// (`LammpsDumpFrame::cartesian_positions`), flattened to `[x0,y0,z0,
/// x1,y1,z1,...]`. `Err` (not `None`/`null`) when no recognized coordinate
/// columns are present/resolvable -- see [`lammps_dump_cartesian_positions_f64`]'s
/// doc comment for why.
pub(super) fn dump_cartesian_flat(frame_json: &str) -> Result<Vec<f64>, String> {
    let frame = dump_frame_from_json_str(frame_json)?;
    let positions = frame.cartesian_positions().ok_or_else(|| {
        "no recognized coordinate columns (x/y/z, xs/ys/zs, or xu/yu/zu present but unresolvable)"
            .to_string()
    })?;
    Ok(positions.into_iter().flatten().collect())
}

/// Flat `values` from a Gaussian Cube file's grid, as a `Float64Array` --
/// same data [`cube_grid_json`]'s `"values"` field carries (row-major,
/// third-axis-fastest order -- see `chematic_mol::volumetric`'s module
/// docs for the exact index formula), as a real typed array instead of a
/// JSON number array.
#[wasm_bindgen]
pub fn cube_values_f64(text: &str) -> Result<js_sys::Float64Array, JsValue> {
    check_input_len("cube input", text)?;
    let values = cube_values_vec(text).map_err(|e| JsValue::from_str(&e))?;
    Ok(js_sys::Float64Array::from(values.as_slice()))
}

/// `[nx, ny, nz]` for a Gaussian Cube file's grid, as a `Uint32Array`.
#[wasm_bindgen]
pub fn cube_shape_u32(text: &str) -> Result<js_sys::Uint32Array, JsValue> {
    check_input_len("cube input", text)?;
    let shape = cube_shape_vec(text).map_err(|e| JsValue::from_str(&e))?;
    Ok(js_sys::Uint32Array::from(shape.as_slice()))
}

/// Flat `values` from an OpenDX file's grid, as a `Float64Array` -- same
/// data [`opendx_grid_json`]'s `"values"` field carries.
#[wasm_bindgen]
pub fn opendx_values_f64(text: &str) -> Result<js_sys::Float64Array, JsValue> {
    check_input_len("OpenDX input", text)?;
    let values = opendx_values_vec(text).map_err(|e| JsValue::from_str(&e))?;
    Ok(js_sys::Float64Array::from(values.as_slice()))
}

/// `[nx, ny, nz]` for an OpenDX file's grid, as a `Uint32Array`.
#[wasm_bindgen]
pub fn opendx_shape_u32(text: &str) -> Result<js_sys::Uint32Array, JsValue> {
    check_input_len("OpenDX input", text)?;
    let shape = opendx_shape_vec(text).map_err(|e| JsValue::from_str(&e))?;
    Ok(js_sys::Uint32Array::from(shape.as_slice()))
}

/// Flattens a LAMMPS dump frame's `rows` (JSON shape
/// [`lammps_dump_frame_to_json_str`] returns) into a single flat
/// `Float64Array`, row-major (atom 0's `column_names.len()` values, then
/// atom 1's, ...). The caller already has `column_names` from
/// [`lammps_dump_frame_to_json_str`] and can compute the row length
/// itself (`column_names.length`); no separate row-length accessor is
/// provided here.
#[wasm_bindgen]
pub fn lammps_dump_rows_f64(frame_json: &str) -> Result<js_sys::Float64Array, JsValue> {
    check_json_len("LAMMPS dump frame JSON", frame_json)?;
    let flat = dump_rows_flat(frame_json).map_err(|e| JsValue::from_str(&e))?;
    Ok(js_sys::Float64Array::from(flat.as_slice()))
}

/// Like [`lammps_dump_cartesian_positions_json`], but returns a flat
/// `Float64Array` (`[x0,y0,z0,x1,y1,z1,...]`, 3 values per atom) instead
/// of a JSON `[[x,y,z],...]` array.
///
/// **Behavioral difference from the JSON sibling**: when the frame has no
/// recognized coordinate columns, [`lammps_dump_cartesian_positions_json`]
/// returns JSON `null`; a `Float64Array` has no `null`, so this function
/// returns `Err` instead, with a message naming the columns it looked for.
#[wasm_bindgen]
pub fn lammps_dump_cartesian_positions_f64(
    frame_json: &str,
) -> Result<js_sys::Float64Array, JsValue> {
    check_json_len("LAMMPS dump frame JSON", frame_json)?;
    let flat = dump_cartesian_flat(frame_json).map_err(|e| JsValue::from_str(&e))?;
    Ok(js_sys::Float64Array::from(flat.as_slice()))
}
