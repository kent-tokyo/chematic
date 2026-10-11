use super::*;

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
/// Coordinates and format-specific metadata are intentionally not carried
/// across this topology-only bridge; use the format-specific coordinate APIs
/// when those fields must be preserved.
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
