//! Molecule I/O format bindings: SDF/MOL, CML, CDXML, MOL2, XYZ, PDB, PDBQT, InChI.

use crate::{
    MolHandle, WASM_MAX_ATOMS, WASM_MAX_BATCH_ITEMS, WASM_MAX_INPUT_BYTES,
    WASM_MAX_JSON_STRING_BYTES, enforce_wasm_molecule_size, escape_json_string,
    parse_smiles_json_array, parse_wasm_string_json_array,
};
use wasm_bindgen::prelude::*;

/// Generate 3D coordinates for the molecule and return a PDB string.
///
/// Coordinates are generated using distance-geometry placement with ring templates.
/// Returns heavy-atom PDB (HETATM records, no explicit H).
#[wasm_bindgen]
pub fn generate_3d_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    chematic_3d::write_pdb(&mol.inner, &coords)
}

/// Generate energy-minimized 3D coordinates and return a PDB string.
///
/// Runs distance-geometry placement followed by gradient-descent force-field
/// minimization.  Geometry quality is better than `generate_3d_pdb` for
/// flexible molecules; the force field is approximate (not MMFF94/UFF).
#[wasm_bindgen]
pub fn generate_3d_minimized_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let minimized = chematic_3d::minimize(&mol.inner, coords);
    chematic_3d::write_pdb(&mol.inner, &minimized)
}

/// Generate 3D coordinates using ETKDG (torsion angle preferences) and return PDB block.
/// ETKDG produces higher-quality conformations than rule-based DG by applying
/// experimental torsion angle preferences to common structural patterns.
#[wasm_bindgen]
pub fn generate_3d_etkdg_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords_etkdg(&mol.inner);
    chematic_3d::write_pdb(&mol.inner, &coords)
}

/// Generate 3D coordinates using ETKDG and minimize with DREIDING force field.
#[wasm_bindgen]
pub fn generate_3d_etkdg_minimized_pdb(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords_etkdg(&mol.inner);
    let minimized = chematic_3d::minimize(&mol.inner, coords);
    chematic_3d::write_pdb(&mol.inner, &minimized)
}

/// Stable, typed reason string for a [`chematic_perception::StereoDiagnostic`]
/// -- matches the vocabulary the Python bindings use (never a free-form message).
fn stereo_reason_str(reason: chematic_perception::StereoRejectionReason) -> &'static str {
    use chematic_perception::StereoRejectionReason::*;
    match reason {
        ContradictoryWedges => "contradictory_wedges",
        MissingCoordinate => "missing_coordinate",
        DegenerateGeometry => "degenerate_geometry",
        UnsupportedCoordination => "unsupported_coordination",
    }
}

/// Format rejected wedge/hash stereocenters as a JSON array of
/// `{"atom_idx":N,"reason":"..."}` objects.
fn stereo_diagnostics_json(diagnostics: &[chematic_perception::StereoDiagnostic]) -> String {
    let parts: Vec<String> = diagnostics
        .iter()
        .map(|d| {
            format!(
                "{{\"atom_idx\":{},\"reason\":\"{}\"}}",
                d.atom.0,
                stereo_reason_str(d.reason)
            )
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Rejected wedge/hash stereocenters for a MOL V2000 block, as JSON.
///
/// Companion to [`mol_from_sdf_block`]/[`mol_block_coords_json`] -- returns
/// `[{"atom_idx":N,"reason":"..."}]`, empty unless a wedge/hash bond was
/// present at some center and got rejected. See
/// `crates/chematic-py/src/formats.rs`'s `from_mol_block_with_diagnostics`
/// for the reason vocabulary (kept identical across bindings).
#[wasm_bindgen]
pub fn mol_block_stereo_diagnostics_json(mol_block: &str) -> Result<String, JsValue> {
    let report = chematic_mol::read_mol_with_diagnostics(mol_block)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(stereo_diagnostics_json(&report.stereo_diagnostics))
}

/// Rejected wedge/hash stereocenters for a MOL V3000 block, as JSON. Same
/// shape as [`mol_block_stereo_diagnostics_json`] but for V3000 input.
#[wasm_bindgen]
pub fn mol_v3000_stereo_diagnostics_json(block: &str) -> Result<String, JsValue> {
    let report = chematic_mol::read_mol_v3000_with_diagnostics(block)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(stereo_diagnostics_json(&report.stereo_diagnostics))
}

/// Parse a MOL V3000 block and return a `MolHandle`.
///
/// Returns a JS error string on parse failure.
#[wasm_bindgen]
pub fn mol_from_v3000_block(block: &str) -> Result<MolHandle, JsValue> {
    if block.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("V3000 block too large"));
    }
    let (mol, _meta) =
        chematic_mol::parse_mol_v3000(block).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {} atoms)",
            WASM_MAX_ATOMS
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Parse a MOL V2000 block and return a `MolHandle`.
///
/// Returns a JS error string on parse failure.
#[wasm_bindgen]
pub fn mol_from_sdf_block(block: &str) -> Result<MolHandle, JsValue> {
    if block.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("MOL block too large"));
    }
    let (mol, _meta) =
        chematic_mol::parse_mol(block).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {} atoms)",
            WASM_MAX_ATOMS
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Serialize a molecule to a MOL V2000 block with 2D coordinates.
///
/// Atom positions are computed via the same layout engine used for SVG depiction
/// and converted to Ångström units (`1.5 Å` per bond).
#[wasm_bindgen]
pub fn to_mol_block(mol: &MolHandle) -> String {
    let meta = chematic_mol::MolMetadata::default();
    let layout = chematic_depict::compute_layout(&mol.inner);
    // SVG bond length is 40 px; standard MOL bond length is 1.5 Å.
    // Negate Y because SVG y-axis points down, MOL y-axis points up.
    const SCALE: f64 = 1.5 / 40.0;
    let coords: Vec<(f64, f64)> = (0..mol.inner.atom_count())
        .map(|i| {
            let p = layout.get(chematic_core::AtomIdx(i as u32));
            (p.x * SCALE, -p.y * SCALE)
        })
        .collect();
    chematic_mol::write_mol_with_coords(&mol.inner, &meta, &coords)
}

/// Parse an SDF string and return a JSON array of canonical SMILES strings.
///
/// Invalid records are represented as `null` in the array.
#[wasm_bindgen]
pub fn sdf_to_smiles_json(sdf: &str) -> String {
    if sdf.len() > WASM_MAX_INPUT_BYTES {
        return format!(
            r#"[{{"error":"SDF input too large ({} bytes)"}}]"#,
            sdf.len()
        );
    }
    let entries: Vec<String> = chematic_mol::SdfReader::new(sdf)
        .take(WASM_MAX_BATCH_ITEMS)
        .map(|r| match r {
            Ok((mol, _)) => {
                let smi = chematic_smiles::canonical_smiles(&mol);
                format!("\"{}\"", smi.replace('"', "\\\""))
            }
            Err(_) => "null".to_string(),
        })
        .collect();
    format!("[{}]", entries.join(","))
}

// ---------------------------------------------------------------------------
// EState free functions (Sprint P)
// ---------------------------------------------------------------------------

/// Serialize a SMILES string directly to a MOL V2000 block with 2D coordinates.
///
/// Returns a JS error on SMILES parse failure.
#[wasm_bindgen]
pub fn mol_block_from_smiles(smiles: &str) -> Result<String, JsValue> {
    let mol = chematic_smiles::parse(smiles)
        .map(|m| MolHandle {
            inner: std::rc::Rc::new(m),
        })
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(to_mol_block(&mol))
}

/// Parse an SDF string and return a JSON array of record objects.
///
/// Each record has the shape:
/// ```json
/// {"smiles":"CC(=O)O","name":"aspirin","properties":{"MW":"180.2","Activity":"high"},"stereo_diagnostics":[]}
/// ```
///
/// Invalid records are represented as `null`.  SD data fields are included in
/// `properties`; multi-line values are joined with `\n`. `stereo_diagnostics`
/// is a list of `{"atom_idx":N,"reason":"..."}` objects, one per rejected
/// wedge/hash center (see [`mol_block_stereo_diagnostics_json`] for the
/// reason vocabulary) -- empty unless a wedge/hash bond was present at some
/// center and got rejected.
#[wasm_bindgen]
pub fn sdf_to_records_json(sdf: &str) -> String {
    if sdf.len() > WASM_MAX_INPUT_BYTES {
        return format!(
            r#"[{{"error":"SDF input too large ({} bytes)"}}]"#,
            sdf.len()
        );
    }
    let entries: Vec<String> = chematic_mol::SdfRecordReader::new(sdf)
        .take(WASM_MAX_BATCH_ITEMS)
        .map(|r| match r {
            Ok(rec) => {
                let smi = chematic_smiles::canonical_smiles(&rec.mol);
                let name = escape_json_string(&rec.meta.name);
                let props: Vec<String> = rec
                    .properties
                    .iter()
                    .map(|(k, v)| {
                        format!(
                            "\"{}\":\"{}\"",
                            escape_json_string(k),
                            escape_json_string(v)
                        )
                    })
                    .collect();
                format!(
                    "{{\"smiles\":\"{}\",\"name\":\"{}\",\"properties\":{{{}}},\"stereo_diagnostics\":{}}}",
                    escape_json_string(&smi),
                    name,
                    props.join(","),
                    stereo_diagnostics_json(&rec.stereo_diagnostics)
                )
            }
            Err(_) => "null".to_string(),
        })
        .collect();
    format!("[{}]", entries.join(","))
}

/// Serialise a JSON array of SMILES to an SDF string.
///
/// Generates 2D coordinates for each molecule.  Property data can be
/// included by using `sdf_from_records_json` instead.
#[wasm_bindgen]
#[allow(clippy::type_complexity)]
pub fn smiles_array_to_sdf(smiles_json: &str) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mut records: Vec<(
        chematic_core::Molecule,
        chematic_mol::MolMetadata,
        Vec<(f64, f64)>,
    )> = Vec::new();

    for smi in &smiles_list {
        let mol = chematic_smiles::parse(smi).map_err(|e| JsValue::from_str(&e.to_string()))?;
        enforce_wasm_molecule_size(&mol)?;
        let layout = chematic_depict::compute_layout(&mol);
        const SCALE: f64 = 1.5 / 40.0;
        let coords: Vec<(f64, f64)> = (0..mol.atom_count())
            .map(|i| {
                let p = layout.get(chematic_core::AtomIdx(i as u32));
                (p.x * SCALE, -p.y * SCALE)
            })
            .collect();
        let meta = chematic_mol::MolMetadata::default();
        records.push((mol, meta, coords));
    }

    let refs: Vec<(
        &chematic_core::Molecule,
        &chematic_mol::MolMetadata,
        &[(f64, f64)],
    )> = records
        .iter()
        .map(|(m, meta, c)| (m, meta, c.as_slice()))
        .collect();

    Ok(chematic_mol::write_sdf(&refs))
}

/// Serialise a `MolHandle` to MOL V3000 format with 2D coordinates.
#[wasm_bindgen]
pub fn to_mol_v3000_block(mol: &MolHandle) -> String {
    let layout = chematic_depict::compute_layout(&mol.inner);
    const SCALE: f64 = 1.5 / 40.0;
    let coords: Vec<(f64, f64)> = (0..mol.inner.atom_count())
        .map(|i| {
            let p = layout.get(chematic_core::AtomIdx(i as u32));
            (p.x * SCALE, -p.y * SCALE)
        })
        .collect();
    let meta = chematic_mol::MolMetadata::default();
    chematic_mol::write_mol_v3000(&mol.inner, &meta, &coords)
}

// ---------------------------------------------------------------------------
// DepictData
// ---------------------------------------------------------------------------

/// Parse a CML string into a `MolHandle`.
///
/// Returns a JS error if the CML is invalid (unknown element, bad bond, etc.).
#[wasm_bindgen]
pub fn mol_from_cml(cml: &str) -> Result<MolHandle, JsValue> {
    if cml.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("CML input too large"));
    }
    let (mol, _coords) =
        chematic_mol::parse_cml(cml).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {} atoms)",
            WASM_MAX_ATOMS
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Serialise a `MolHandle` to a CML string with 2D coordinates.
///
/// Coordinates are generated using the same 2D layout engine as `to_mol_block`.
#[wasm_bindgen]
pub fn to_cml(mol: &MolHandle) -> String {
    let layout = chematic_depict::compute_layout(&mol.inner);
    const SCALE: f64 = 1.5 / 40.0;
    let coords: Vec<(f64, f64)> = (0..mol.inner.atom_count())
        .map(|i| {
            let p = layout.get(chematic_core::AtomIdx(i as u32));
            (p.x * SCALE, -p.y * SCALE)
        })
        .collect();
    chematic_mol::write_cml(&mol.inner, Some(&coords))
}

/// Parse a MolJSON string into a `MolHandle`.
///
/// MolJSON is a JSON-based molecular representation designed for LLM
/// (large language model) compatibility.  Returns a JS error on invalid input.
#[wasm_bindgen]
pub fn mol_from_moljson(json: &str) -> Result<MolHandle, JsValue> {
    if json.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("MolJSON input too large"));
    }
    let mol = chematic_mol::parse_moljson(json).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {} atoms)",
            WASM_MAX_ATOMS
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Serialise a `MolHandle` to a MolJSON string (pretty-printed).
///
/// Atom IDs are assigned as `"a1"`, `"a2"`, … in molecule atom order.
/// The `hydrogens` field reflects computed implicit H count.
#[wasm_bindgen]
pub fn to_moljson(mol: &MolHandle) -> String {
    chematic_mol::write_moljson(&mol.inner)
}

/// Only the first molecular fragment in the document is returned.
/// Returns a JS error if the document cannot be parsed.
#[wasm_bindgen]
pub fn mol_from_cdxml(cdxml: &str) -> Result<MolHandle, JsValue> {
    if cdxml.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("CDXML input too large"));
    }
    let (mol, _coords) =
        chematic_mol::parse_cdxml(cdxml).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {} atoms)",
            WASM_MAX_ATOMS
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

/// Parse all molecular fragments from a CDXML string.
///
/// Returns a JSON array of SMILES strings, one per fragment:
/// `["CC","c1ccccc1"]`
///
/// Stereochemistry (wedge/dash bonds) is read from the `Display` attribute
/// of bond elements.
#[wasm_bindgen]
pub fn cdxml_to_smiles_json(cdxml: &str) -> Result<String, JsValue> {
    if cdxml.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("CDXML input too large"));
    }
    let fragments =
        chematic_mol::parse_cdxml_all(cdxml).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let parts: Vec<String> = fragments
        .iter()
        .take(WASM_MAX_BATCH_ITEMS)
        .map(|(mol, _)| {
            format!(
                "\"{}\"",
                escape_json_string(&chematic_smiles::canonical_smiles(mol))
            )
        })
        .collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Parse a MOL V2000 string and return 2D coordinates as a JSON array.
///
/// Returns `[[x0,y0],[x1,y1],...]` in atom-insertion order.
/// Coordinates are in Ångström as stored in the MOL file.
#[wasm_bindgen]
pub fn mol_block_coords_json(mol_block: &str) -> Result<String, JsValue> {
    let (_mol, _meta, coords) = chematic_mol::parse_mol_with_coords(mol_block)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let parts: Vec<String> = coords
        .iter()
        .map(|(x, y)| {
            let xf = if x.is_finite() {
                format!("{x:.4}")
            } else {
                "0".to_string()
            };
            let yf = if y.is_finite() {
                format!("{y:.4}")
            } else {
                "0".to_string()
            };
            format!("[{xf},{yf}]")
        })
        .collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Serialize multiple molecules with properties to an SDF string.
///
/// # Arguments
/// * `smiles_json` — JSON array of SMILES strings, e.g. `["CC(=O)O","c1ccccc1"]`
/// * `names_json`  — JSON array of molecule names (same length as `smiles_json`)
/// * `props_json`  — JSON array where each element encodes one molecule's SD data fields
///   as `"key1\tvalue1\nkey2\tvalue2"` (tab-separated key/value, `\n`-separated pairs;
///   pass `""` for a molecule with no properties)
///
/// Returns the SDF string, or a JS error if any SMILES fails to parse or the
/// arrays have mismatched lengths.
///
/// The `\n` and `\t` sequences in `props_json` are JSON-escaped — they are
/// decoded to the actual characters before SDF formatting.
#[wasm_bindgen]
pub fn sdf_from_records_json(
    smiles_json: &str,
    names_json: &str,
    props_json: &str,
) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let names_list = parse_wasm_string_json_array(names_json, "names_json")?;
    let props_list = parse_wasm_string_json_array(props_json, "props_json")?;

    let n = smiles_list.len();
    if names_list.len() != n || props_list.len() != n {
        return Err(JsValue::from_str(
            "sdf_from_records_json: smiles_json, names_json, and props_json must have the same length",
        ));
    }

    let mut sdf = String::new();
    for i in 0..n {
        let mol = chematic_smiles::parse(&smiles_list[i])
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        enforce_wasm_molecule_size(&mol)?;

        // Build 2D coordinates using the same pipeline as to_mol_block.
        let layout = chematic_depict::compute_layout(&mol);
        const SCALE: f64 = 1.5 / 40.0;
        let coords: Vec<(f64, f64)> = (0..mol.atom_count())
            .map(|j| {
                let p = layout.get(chematic_core::AtomIdx(j as u32));
                (p.x * SCALE, -p.y * SCALE)
            })
            .collect();

        let meta = chematic_mol::MolMetadata {
            name: names_list[i].clone(),
            comment: String::new(),
        };
        sdf.push_str(&chematic_mol::write_mol_with_coords(&mol, &meta, &coords));

        // SD data fields from "key\tvalue\nkey\tvalue" format.
        for pair in props_list[i].lines() {
            if let Some((key, value)) = pair.split_once('\t') {
                sdf.push_str(&format!("> <{key}>\n{value}\n\n"));
            }
        }

        sdf.push_str("$$$$\n");
    }

    Ok(sdf)
}

// ---------------------------------------------------------------------------
// Sprint Y: XYZ/PDB I/O, per-atom descriptors, SSSR rings,
//           custom ECFP, stereo isomer enumeration
// ---------------------------------------------------------------------------

/// Parse an XYZ file and return a `MolHandle` (topology only; coordinates are discarded).
///
/// Returns a JS error on parse failure.
#[wasm_bindgen]
pub fn mol_from_xyz(xyz: &str) -> Result<MolHandle, JsValue> {
    if xyz.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str("XYZ input too large"));
    }
    let (mol, _coords) =
        chematic_3d::parse_xyz(xyz).map_err(|e| JsValue::from_str(&format!("{e:?}")))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "molecule too large (max {} atoms)",
            WASM_MAX_ATOMS
        )));
    }
    Ok(MolHandle {
        inner: std::rc::Rc::new(mol),
    })
}

// determine_bonds_from_xyz_json removed to reduce WASM bundle size.

/// Serialize a molecule to XYZ format.
///
/// 3D coordinates are generated via distance-geometry placement.
#[wasm_bindgen]
pub fn to_xyz(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    chematic_3d::write_xyz(&mol.inner, &coords, "")
}

/// Why a PDB block was rejected by [`parse_pdb_molecule_and_coords`].
pub(crate) enum PdbInputError {
    TooLarge,
    TooManyAtoms,
}

impl std::fmt::Display for PdbInputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooLarge => write!(
                f,
                "PDB input too large (max {WASM_MAX_JSON_STRING_BYTES} bytes)"
            ),
            Self::TooManyAtoms => write!(f, "molecule too large (max {WASM_MAX_ATOMS} atoms)"),
        }
    }
}

/// Whether every component of every coordinate is finite (not NaN or
/// Infinity). Shared by `pdb_coords_json` and
/// `mmff94_energy_breakdown_from_coords_json` -- `f64::from_str` parses the
/// literal text "nan"/"inf"/"infinity" (case-insensitively) into their
/// special values rather than erroring, so a PDB coordinate FIELD
/// containing that text (not just truly malformed text, which falls back to
/// 0.0 via `unwrap_or`) passes `parse_pdb_atoms` silently. Without this
/// check, `pdb_coords_json` would serialize a bare `NaN`/`inf` token, which
/// is not valid JSON syntax, into its output.
pub(crate) fn coords_all_finite(coords: &[[f64; 3]]) -> bool {
    coords.iter().flatten().all(|v: &f64| v.is_finite())
}

/// Parse a PDB block into topology and coordinates from a SINGLE pass over
/// the input, shared by `mol_from_pdb` and `pdb_coords_json` -- both must
/// read the exact same atom list in the exact same order, or a caller
/// combining their two outputs would silently pair a topology with the
/// wrong atom's coordinates (issue #90's regression, one level up: the
/// original bug was mol/coords never travelling together at all).
pub(crate) fn parse_pdb_molecule_and_coords(
    pdb: &str,
) -> Result<(chematic_core::Molecule, Vec<[f64; 3]>), PdbInputError> {
    if pdb.len() > WASM_MAX_JSON_STRING_BYTES {
        return Err(PdbInputError::TooLarge);
    }
    let atoms = chematic_3d::parse_pdb_atoms(pdb);
    if atoms.len() > WASM_MAX_ATOMS {
        return Err(PdbInputError::TooManyAtoms);
    }
    let (mol, coords) = chematic_3d::pdb_to_molecule(&atoms);
    let flat: Vec<[f64; 3]> = coords.points.iter().map(|p| [p.x, p.y, p.z]).collect();
    Ok((mol, flat))
}

/// Parse a PDB file and return a `MolHandle` (topology only; coordinates are
/// discarded -- use [`pdb_coords_json`] to recover them in the SAME atom
/// order, and [`mmff94_energy_breakdown_from_coords_json`] to score them
/// without chematic regenerating a fresh conformer).
///
/// Uses CONECT records for connectivity if present; otherwise infers bonds from
/// atom distances (the same heuristic as the internal `pdb_to_molecule` function).
#[wasm_bindgen]
pub fn mol_from_pdb(pdb: &str) -> MolHandle {
    match parse_pdb_molecule_and_coords(pdb) {
        Ok((mol, _coords)) => MolHandle {
            inner: std::rc::Rc::new(mol),
        },
        Err(_) => MolHandle {
            inner: std::rc::Rc::new(chematic_core::MoleculeBuilder::new().build()),
        },
    }
}

/// Extract the atomic coordinates from a PDB block, in the SAME atom order
/// `mol_from_pdb` returns topology for (both read the identical underlying
/// parse via `parse_pdb_molecule_and_coords`, so atom-index correspondence
/// between the two calls is structural, not just conventional -- issue #90).
///
/// Returns JSON `[[x,y,z],...]` (full `f64` precision, not rounded -- for
/// oracle comparison against the Python binding) or `{"error":"<msg>"}` --
/// including when a coordinate field parsed to a non-finite value (see
/// `coords_all_finite`'s doc comment), rather than emitting invalid JSON.
#[wasm_bindgen]
pub fn pdb_coords_json(pdb: &str) -> String {
    match parse_pdb_molecule_and_coords(pdb) {
        Ok((_mol, coords)) => {
            if !coords_all_finite(&coords) {
                return r#"{"error":"PDB contains a non-finite coordinate (NaN or Infinity)"}"#
                    .to_string();
            }
            let parts: Vec<String> = coords
                .iter()
                .map(|c| format!("[{},{},{}]", c[0], c[1], c[2]))
                .collect();
            format!("[{}]", parts.join(","))
        }
        Err(e) => format!(r#"{{"error":"{e}"}}"#),
    }
}

/// Parse a Tripos MOL2 string and return SMILES.
///
/// Returns `"error:<msg>"` on failure.
#[wasm_bindgen]
pub fn mol2_to_smiles(mol2_str: &str) -> String {
    match chematic_mol::parse_mol2(mol2_str) {
        Ok((mol, _)) => chematic_smiles::write(&mol),
        Err(e) => format!("error:{e}"),
    }
}

/// Convert a SMILES to a minimal Tripos MOL2 string (no 3D coordinates).
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn smiles_to_mol2(smiles: &str) -> String {
    match chematic_smiles::parse(smiles) {
        Ok(mol) => chematic_mol::write_mol2(&mol, &[]),
        Err(e) => format!("error:{e}"),
    }
}

/// Write a molecule to AutoDock PDBQT format.
///
/// `coords_json` — JSON array of `[x,y,z]` arrays (Å). Pass `"[]"` for zero coords.
/// `charges_json` — JSON array of partial charges. Pass `"[]"` to write zeros.
/// `name` — ligand name for the REMARK header.
///
/// Returns the PDBQT string, or `"error:<msg>"` on failure.
#[wasm_bindgen]
pub fn smiles_to_pdbqt(smiles: &str, coords_json: &str, charges_json: &str, name: &str) -> String {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:{e}"),
    };
    let coords: Vec<(f64, f64, f64)> = serde_json::from_str::<Vec<[f64; 3]>>(coords_json)
        .unwrap_or_default()
        .into_iter()
        .map(|c| (c[0], c[1], c[2]))
        .collect();
    let charges: Vec<f64> = serde_json::from_str(charges_json).unwrap_or_default();
    chematic_mol::write_pdbqt(&mol, &coords, &charges, name)
}

/// Generate InChI string from SMILES.
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn inchi_from_smiles(smiles: &str) -> String {
    match chematic_smiles::parse(smiles) {
        Ok(mol) => chematic_inchi::inchi(&mol),
        Err(e) => format!("error:{e}"),
    }
}

/// Generate InChIKey from SMILES (27-character identifier).
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn inchikey_from_smiles(smiles: &str) -> String {
    match chematic_smiles::parse(smiles) {
        Ok(mol) => {
            let inchi_str = chematic_inchi::inchi(&mol);
            chematic_inchi::inchi_key(&inchi_str)
        }
        Err(e) => format!("error:{e}"),
    }
}
