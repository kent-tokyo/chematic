//! Known graph fields must not be silently coerced, narrowed or defaulted.
use super::MolJsonError;
use chematic_core::BondOrder;
use serde_json::Value;

pub(super) fn invalid(message: impl Into<String>) -> MolJsonError {
    MolJsonError::InvalidJson(message.into())
}

pub(super) fn flag(value: &Value, field: &'static str) -> Result<Option<bool>, MolJsonError> {
    value
        .get("aromatic")
        .map(|v| {
            v.as_bool()
                .ok_or_else(|| invalid(format!("{field} must be a boolean")))
        })
        .transpose()
}

pub(super) fn charge(value: &Value) -> Result<i8, MolJsonError> {
    value.get("charge").map_or(Ok(0), |v| {
        v.as_i64()
            .and_then(|n| i8::try_from(n).ok())
            .ok_or_else(|| invalid("atoms[].charge must be an integer in -128..=127"))
    })
}

pub(super) fn isotope(value: &Value) -> Result<Option<u16>, MolJsonError> {
    match value.get("isotope") {
        None | Some(Value::Null) => Ok(None),
        Some(v) => v
            .as_u64()
            .and_then(|n| u16::try_from(n).ok())
            .map(Some)
            .ok_or_else(|| invalid("atoms[].isotope must be null or an integer in 0..=65535")),
    }
}

pub(super) fn id<'a>(
    value: &'a Value,
    key: &'static str,
    field: &'static str,
) -> Result<&'a str, MolJsonError> {
    value
        .get(key)
        .ok_or(MolJsonError::MissingField(field))?
        .as_str()
        .filter(|s| !s.is_empty())
        .ok_or_else(|| invalid(format!("{field} must be a nonempty string")))
}

pub(super) fn bond_order(value: &Value) -> Result<BondOrder, MolJsonError> {
    let aromatic = flag(value, "bonds[].aromatic")?;
    let Some(raw) = value.get("order") else {
        return Ok(if aromatic == Some(true) {
            BondOrder::Aromatic
        } else {
            BondOrder::Single
        });
    };
    let order = match raw.as_f64() {
        Some(1.0) => BondOrder::Single,
        Some(1.5) => BondOrder::Aromatic,
        Some(2.0) => BondOrder::Double,
        Some(3.0) => BondOrder::Triple,
        Some(4.0) => BondOrder::Quadruple,
        _ => return Err(invalid("bonds[].order must be exactly 1, 1.5, 2, 3 or 4")),
    };
    if aromatic.is_some_and(|a| a != (order == BondOrder::Aromatic)) {
        return Err(invalid("bonds[].order conflicts with bonds[].aromatic"));
    }
    Ok(order)
}
