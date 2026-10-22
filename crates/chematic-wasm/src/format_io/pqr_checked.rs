use super::*;

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Record {
    group_pdb: String,
    serial: i64,
    atom_name: String,
    res_name: String,
    chain_id: Option<String>,
    res_seq: i64,
    icode: Option<char>,
    x: f64,
    y: f64,
    z: f64,
    charge: f64,
    radius: f64,
    element: Option<String>,
}
impl Record {
    fn into_atom(self) -> Result<chematic_mol::PqrAtomRecord, String> {
        let element = match self.element {
            Some(symbol) => element_from_symbol(&symbol)?,
            None => chematic_mol::infer_element(&self.group_pdb, &self.res_name, &self.atom_name)
                .ok_or("cannot infer PQR element")?,
        };
        Ok(chematic_mol::PqrAtomRecord {
            group_pdb: self.group_pdb,
            serial: self.serial,
            atom_name: self.atom_name,
            res_name: self.res_name,
            chain_id: self.chain_id,
            res_seq: self.res_seq,
            icode: self.icode,
            x: self.x,
            y: self.y,
            z: self.z,
            charge: self.charge,
            radius: self.radius,
            element,
        })
    }
}

/// Write explicit PQR atom records without numeric rounding or bond inference.
/// Strict JSON array: no duplicate/unknown fields, missing required values or
/// incompatible elements. Non-atom source metadata is not represented.
#[wasm_bindgen]
pub fn write_pqr_json_checked(records_json: &str) -> Result<String, JsValue> {
    check_input_len("PQR records JSON", records_json)?;
    let records: Vec<Record> = serde_json::from_str(records_json)
        .map_err(|error| JsValue::from_str(&format!("PQR records JSON: {error}")))?;
    if records.len() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str("PQR output atom limit exceeded"));
    }
    let atoms = records
        .into_iter()
        .map(Record::into_atom)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| JsValue::from_str(&error))?;
    chematic_mol::write_pqr_checked_with_limits(
        &atoms,
        &chematic_mol::PqrWriteLimits {
            max_atoms: WASM_MAX_ATOMS,
            max_output_bytes: WASM_MAX_OUTPUT_BYTES,
        },
    )
    .map_err(|error| JsValue::from_str(&error.to_string()))
}
