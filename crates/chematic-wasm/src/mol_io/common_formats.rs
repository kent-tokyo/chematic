use super::*;

/// Replace only ar bond tokens in a bounded neutral aromatic MOL2 record.
/// Retains atom types and every other source byte; errors throw, no partial text.
#[wasm_bindgen]
pub fn kekulize_mol2_record(text: &str) -> Result<String, JsValue> {
    chematic_mol::kekulize_mol2_record(text).map_err(|error| JsValue::from_str(&error.to_string()))
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
/// Optional `atom_bond_adjacent=true` emits ATOM → BOND → UNITY. Some
/// readers (including Open Babel 3.2.1) ignore UNITY charges in this order.
/// Omit it to retain existing bytes. Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn smiles_to_mol2(smiles: &str, atom_bond_adjacent: Option<bool>) -> String {
    let order = if atom_bond_adjacent.unwrap_or(false) {
        chematic_mol::Mol2SectionOrder::AtomBondAdjacent
    } else {
        chematic_mol::Mol2SectionOrder::default()
    };
    match chematic_smiles::parse(smiles) {
        Ok(mol) => chematic_mol::write_mol2_with_section_order(&mol, &[], order),
        Err(e) => format!("error:{e}"),
    }
}

fn bounded_mol2_record(text: &str) -> Result<chematic_mol::Mol2SourceRecord, String> {
    chematic_mol::parse_mol2_source_record_with_limits(
        text,
        &chematic_mol::Mol2ParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_atoms: WASM_MAX_ATOMS,
            ..Default::default()
        },
    )
    .map_err(|error| error.to_string())
}

fn bounded_cjson_record(text: &str) -> Result<chematic_mol::CjsonRecord, String> {
    chematic_mol::parse_cjson_record_with_limits(
        text,
        &chematic_mol::CjsonParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_atoms: WASM_MAX_ATOMS,
            ..Default::default()
        },
    )
    .map_err(|error| error.to_string())
}

/// Inspect strictly validated CJSON with exact source and opaque field paths.
#[wasm_bindgen]
pub fn parse_cjson_record_json(text: &str) -> Result<String, JsValue> {
    let record = bounded_cjson_record(text).map_err(|error| JsValue::from_str(&error))?;
    let json = record.to_json().to_string();
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "CJSON record JSON exceeds maximum output size",
        ));
    }
    Ok(json)
}

/// Validate the local CJSON profile and retain the exact UTF-8 source.
#[wasm_bindgen]
pub fn roundtrip_cjson_record(text: &str) -> Result<String, JsValue> {
    bounded_cjson_record(text)
        .map(|record| chematic_mol::write_cjson_record(&record))
        .map_err(|error| JsValue::from_str(&error))
}

/// Checked MOL2/CJSON conversion using the shared report schema.
#[wasm_bindgen]
pub fn convert_with_report_json(
    text: &str,
    input_format: &str,
    output_format: &str,
    policy: &str,
) -> Result<String, JsValue> {
    let limits = chematic_mol::ConversionLimits {
        mol2: chematic_mol::Mol2ParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_atoms: WASM_MAX_ATOMS,
            ..Default::default()
        },
        max_output_bytes: WASM_MAX_OUTPUT_BYTES,
    };
    chematic_mol::convert_with_report_json(text, input_format, output_format, policy, &limits)
        .map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Inspect partial charges, atom types and source retention for MOL2.
#[wasm_bindgen]
pub fn parse_mol2_record_json(text: &str) -> Result<String, JsValue> {
    let record = bounded_mol2_record(text).map_err(|error| JsValue::from_str(&error))?;
    let json = record.to_json().to_string();
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "MOL2 record JSON exceeds maximum output size",
        ));
    }
    Ok(json)
}

/// Validate one MOL2 record and preserve its exact source, including metadata.
#[wasm_bindgen]
pub fn roundtrip_mol2_record(text: &str) -> Result<String, JsValue> {
    bounded_mol2_record(text)
        .map(|record| chematic_mol::write_mol2_source_record(&record))
        .map_err(|error| JsValue::from_str(&error))
}

pub(crate) fn common_format_name(format: &str) -> Option<String> {
    match format
        .to_ascii_lowercase()
        .trim_start_matches('.')
        .replace('-', "_")
        .as_str()
    {
        "smi" | "smiles" => Some("smiles".to_string()),
        "mol" | "sdf" => Some("mol".to_string()),
        "v3000" | "mol_v3000" => Some("mol_v3000".to_string()),
        "mol2" => Some("mol2".to_string()),
        "cml" => Some("cml".to_string()),
        "cjson" => Some("cjson".to_string()),
        "moljson" => Some("moljson".to_string()),
        "cdxml" => Some("cdxml".to_string()),
        _ => None,
    }
}

fn parse_common_format(text: &str, format: &str) -> Result<chematic_core::Molecule, String> {
    match format {
        "smiles" => chematic_smiles::parse(text).map_err(|e| e.to_string()),
        "mol" => chematic_mol::parse_mol(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string()),
        "mol_v3000" => chematic_mol::parse_mol_v3000(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string()),
        "mol2" => chematic_mol::parse_mol2(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string()),
        "cml" => chematic_mol::parse_cml(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string()),
        "cjson" => chematic_mol::parse_cjson(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string()),
        "moljson" => chematic_mol::parse_moljson(text).map_err(|e| e.to_string()),
        "cdxml" => chematic_mol::parse_cdxml(text)
            .map(|(mol, _)| mol)
            .map_err(|e| e.to_string()),
        _ => Err(format!("unsupported input format: {format}")),
    }
}

/// Convert between common topology-bearing molecular formats.
///
/// This is the WASM counterpart of Python's `convert_format`. It supports
/// SMILES, MOL/SDF, MOL V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML.
/// MOL2-to-MOL2 validates and retains the exact source. Other conversions use
/// the topology-only bridge and do not carry format-specific metadata.
/// CML-to-CJSON carries complete source Cartesian coordinates, leaving absent
/// 3D geometry absent and refusing partial/nonfinite data and multiple molecules.
#[wasm_bindgen]
pub fn convert_common_format(
    text: &str,
    input_format: &str,
    output_format: &str,
) -> Result<String, JsValue> {
    if text.len() > WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "format input exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            text.len()
        )));
    }
    let input = common_format_name(input_format).ok_or_else(|| {
        JsValue::from_str(&format!("unsupported molecular format: {input_format}"))
    })?;
    let output = common_format_name(output_format).ok_or_else(|| {
        JsValue::from_str(&format!("unsupported molecular format: {output_format}"))
    })?;
    if input == "mol2" && output == "mol2" {
        return roundtrip_mol2_record(text);
    }
    if input == "cml" && output == "cjson" {
        let limits = chematic_mol::CmlParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_atoms: WASM_MAX_ATOMS,
            ..Default::default()
        };
        let (mol, coords) = chematic_mol::parse_cml_with_3d_with_limits(text, &limits)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        let points = coords
            .unwrap_or_default()
            .into_iter()
            .map(|[x, y, z]| (x, y, z))
            .collect::<Vec<_>>();
        return Ok(chematic_mol::write_cjson(&mol, &points));
    }
    let mol = parse_common_format(text, &input).map_err(|e| JsValue::from_str(&e))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "format input exceeds maximum atom count ({WASM_MAX_ATOMS})"
        )));
    }
    let result = match output.as_str() {
        "smiles" => chematic_smiles::canonical_smiles(&mol),
        "mol" => chematic_mol::write_mol(&mol, &chematic_mol::MolMetadata::default()),
        "mol_v3000" => {
            chematic_mol::write_mol_v3000(&mol, &chematic_mol::MolMetadata::default(), &[])
        }
        "mol2" => chematic_mol::write_mol2(&mol, &[]),
        "cml" => chematic_mol::write_cml(&mol, None),
        "cjson" => chematic_mol::write_cjson(&mol, &[]),
        "moljson" => chematic_mol::write_moljson(&mol),
        "cdxml" => chematic_mol::write_cdxml(&mol, &[]),
        _ => {
            return Err(JsValue::from_str(&format!(
                "unsupported output format: {output}"
            )));
        }
    };
    Ok(result)
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
    if smiles.len() > WASM_MAX_INPUT_BYTES
        || coords_json.len() > WASM_MAX_JSON_STRING_BYTES
        || charges_json.len() > WASM_MAX_JSON_STRING_BYTES
        || name.len() > WASM_MAX_JSON_STRING_BYTES
    {
        return "error:input exceeds WASM size limits".to_string();
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:{e}"),
    };
    if mol.atom_count() > WASM_MAX_ATOMS {
        return format!("error:molecule too large (max {WASM_MAX_ATOMS} atoms)");
    }
    let coords: Vec<(f64, f64, f64)> = match serde_json::from_str::<Vec<[f64; 3]>>(coords_json) {
        Ok(values) => values.into_iter().map(|c| (c[0], c[1], c[2])).collect(),
        Err(e) => return format!("error:coords JSON parse error: {e}"),
    };
    let charges: Vec<f64> = match serde_json::from_str(charges_json) {
        Ok(values) => values,
        Err(e) => return format!("error:charges JSON parse error: {e}"),
    };
    chematic_mol::write_pdbqt(&mol, &coords, &charges, name)
}

/// Generate InChI string from SMILES.
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn inchi_from_smiles(smiles: &str) -> String {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return format!(
            "error:SMILES exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            smiles.len()
        );
    }
    match chematic_smiles::parse(smiles) {
        Ok(mol) if mol.atom_count() <= WASM_MAX_ATOMS => chematic_inchi::inchi(&mol),
        Ok(_) => format!("error:SMILES exceeds maximum atom count ({WASM_MAX_ATOMS})"),
        Err(e) => format!("error:{e}"),
    }
}

/// Generate InChIKey from SMILES (27-character identifier).
///
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn inchikey_from_smiles(smiles: &str) -> String {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return format!(
            "error:SMILES exceeds maximum input size ({} > {WASM_MAX_INPUT_BYTES} bytes)",
            smiles.len()
        );
    }
    match chematic_smiles::parse(smiles) {
        Ok(mol) if mol.atom_count() <= WASM_MAX_ATOMS => {
            let inchi_str = chematic_inchi::inchi(&mol);
            chematic_inchi::inchi_key(&inchi_str)
        }
        Ok(_) => format!("error:SMILES exceeds maximum atom count ({WASM_MAX_ATOMS})"),
        Err(e) => format!("error:{e}"),
    }
}
