use crate::{MolHandle, WASM_MAX_INPUT_BYTES, WASM_MAX_OUTPUT_BYTES};
use wasm_bindgen::prelude::*;

fn read(
    data: &[u8],
    format: &str,
    allow_metadata_loss: bool,
) -> Result<chematic_mol::ChemDrawBlock, JsValue> {
    let format = chematic_mol::ChemDrawFormat::from_name(format)
        .ok_or_else(|| JsValue::from_str("format must be AUTO, CDX or CDXML"))?;
    let limits = chematic_mol::ChemDrawLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        ..Default::default()
    };
    chematic_mol::parse_chemdraw_block_with_limits(
        data,
        format,
        &chematic_mol::ChemDrawOptions {
            allow_metadata_loss,
        },
        &limits,
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))
}
/// Concrete CDX/CDXML fragment handles, in source fragment order.
#[wasm_bindgen]
pub fn mols_from_chemdraw_block(
    data: &[u8],
    format: &str,
    allow_metadata_loss: bool,
) -> Result<Vec<MolHandle>, JsValue> {
    Ok(read(data, format, allow_metadata_loss)?
        .molecules
        .into_iter()
        .map(|(m, _)| MolHandle {
            inner: std::rc::Rc::new(m),
        })
        .collect())
}
/// Exact source bytes, projected XML, fragment coordinates and CDX omission tags.
/// CDXML uses the existing graph projection; omissions lists CDX tags only.
#[wasm_bindgen]
pub fn chemdraw_block_report_json(
    data: &[u8],
    format: &str,
    allow_metadata_loss: bool,
) -> Result<String, JsValue> {
    let b = read(data, format, allow_metadata_loss)?;
    let coordinates: Vec<_> = b
        .molecules
        .iter()
        .map(|(_, c)| c.iter().map(|&(x, y)| [x, y]).collect::<Vec<_>>())
        .collect();
    let omissions: Vec<_> = b
        .omissions
        .iter()
        .map(|o| serde_json::json!({"offset":o.offset,"tag":o.tag}))
        .collect();
    let result=serde_json::to_string(&serde_json::json!({"format":b.format.name(),"source_bytes":b.source,"projected_cdxml":b.projected_cdxml,"molecule_count":b.molecules.len(),"coordinates":coordinates,"omissions":omissions}))
        .map_err(|e|JsValue::from_str(&e.to_string()))?;
    if result.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("ChemDraw JSON output limit exceeded"));
    }
    Ok(result)
}
