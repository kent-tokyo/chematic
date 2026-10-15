use super::{CcdAtom, CcdBond, CcdError, Row, bond_order, invalid};
use chematic_core::Element;

pub(super) fn take(row: &mut Row, category: &str, field: &str) -> Option<String> {
    row.remove(&format!("_{category}.{field}"))
        .filter(|s| s != "." && s != "?")
}
pub(super) fn required(row: &mut Row, category: &str, field: &str) -> Result<String, CcdError> {
    take(row, category, field)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| invalid(format!("missing {category}.{field}")))
}
fn flag(row: &mut Row, category: &str) -> Result<bool, CcdError> {
    match take(row, category, "pdbx_aromatic_flag").as_deref() {
        None | Some("N") => Ok(false),
        Some("Y") => Ok(true),
        _ => Err(invalid("invalid aromatic flag")),
    }
}
fn stereo(row: &mut Row, category: &str, allowed: &[&str]) -> Result<Option<String>, CcdError> {
    let value = take(row, category, "pdbx_stereo_config");
    if value.as_deref().is_some_and(|s| !allowed.contains(&s)) {
        return Err(invalid("invalid stereo configuration"));
    }
    Ok(value)
}
fn coordinates(row: &mut Row, prefix: &str, suffix: &str) -> Result<Option<[f64; 3]>, CcdError> {
    let mut values = Vec::new();
    for axis in ["x", "y", "z"] {
        if let Some(value) = take(row, "chem_comp_atom", &format!("{prefix}{axis}{suffix}")) {
            let number = value
                .parse::<f64>()
                .ok()
                .filter(|n| n.is_finite())
                .ok_or_else(|| invalid("non-finite or invalid coordinate"))?;
            values.push(number);
        }
    }
    match values.len() {
        0 => Ok(None),
        3 => Ok(Some([values[0], values[1], values[2]])),
        _ => Err(invalid("incomplete coordinate triple")),
    }
}
pub(super) fn atom(mut row: Row, component_id: &str) -> Result<CcdAtom, CcdError> {
    let category = "chem_comp_atom";
    if required(&mut row, category, "comp_id")? != component_id {
        return Err(invalid("atom component ID mismatch"));
    }
    let atom_id = required(&mut row, category, "atom_id")?;
    let symbol = required(&mut row, category, "type_symbol")?;
    let (element, isotope) = element(&symbol)?;
    let charge = take(&mut row, category, "charge")
        .map(|s| s.parse::<i8>().map_err(|_| invalid("invalid atom charge")))
        .transpose()?;
    Ok(CcdAtom {
        atom_id,
        element: element.symbol().into(),
        isotope,
        charge,
        aromatic: flag(&mut row, category)?,
        stereo_config: stereo(&mut row, category, &["N", "R", "S"])?,
        model_coordinates: coordinates(&mut row, "model_cartn_", "")?,
        ideal_coordinates: coordinates(&mut row, "pdbx_model_cartn_", "_ideal")?,
        extra_fields: row,
    })
}
pub(super) fn bond(mut row: Row, component_id: &str) -> Result<CcdBond, CcdError> {
    let category = "chem_comp_bond";
    if required(&mut row, category, "comp_id")? != component_id {
        return Err(invalid("bond component ID mismatch"));
    }
    let atom_id_1 = required(&mut row, category, "atom_id_1")?;
    let atom_id_2 = required(&mut row, category, "atom_id_2")?;
    let value_order = required(&mut row, category, "value_order")?.to_ascii_uppercase();
    bond_order(&value_order)?;
    let aromatic = flag(&mut row, category)? || value_order == "AROM";
    Ok(CcdBond {
        atom_id_1,
        atom_id_2,
        value_order,
        aromatic,
        stereo_config: stereo(&mut row, category, &["N", "E", "Z"])?,
        extra_fields: row,
    })
}

fn element(symbol: &str) -> Result<(Element, Option<u16>), CcdError> {
    let upper = symbol.to_ascii_uppercase();
    match upper.as_str() {
        "D" => return Ok((Element::H, Some(2))),
        "T" => return Ok((Element::H, Some(3))),
        _ => {}
    }
    if !(1..=2).contains(&symbol.len()) || !symbol.bytes().all(|b| b.is_ascii_alphabetic()) {
        return Err(invalid(format!("unknown element {symbol}")));
    }
    let mut normalized = symbol.to_ascii_lowercase();
    normalized[..1].make_ascii_uppercase();
    Element::from_symbol(&normalized)
        .map(|e| (e, None))
        .ok_or_else(|| invalid(format!("unknown element {symbol}")))
}
