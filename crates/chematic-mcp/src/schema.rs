//! Minimal, hand-rolled JSON Schema (2020-12 subset) runtime validator.
//!
//! This is deliberately small: it implements exactly the keywords the 20
//! tool `inputSchema`/`outputSchema` documents in `tools.rs` use (`type`,
//! `properties`, `required`, `additionalProperties`, `items`, `minLength`,
//! `maxLength`, `minItems`, `maxItems`, `minimum`, `maximum`, `enum`,
//! `const`, `oneOf`), not the full 2020-12 vocabulary (no `$ref`, no
//! `$defs`, no `if`/`then`/`else`, no `patternProperties`, no `pattern`).
//! External `$ref` is never resolved — this validator only ever looks at
//! the schema `Value` it was given, so there is no way for it to make a
//! network request even by accident.
//!
//! Used at runtime (before any tool dispatch, modern era only) to map
//! argument-shape problems to `-32602 Invalid Params` per section 9. A
//! *real* general-purpose 2020-12 validator (the `jsonschema` crate) is
//! used as a dev-dependency in tests to independently confirm every
//! `inputSchema`/`outputSchema` in `tools.rs` is itself valid 2020-12, and
//! that sample `structuredContent` values conform to their tool's
//! `outputSchema` — see `tests/schema_conformance.rs`.
//!
//! Validation here is a single recursive walk bounded by the same
//! `MAX_JSON_DEPTH` used elsewhere (`crate::protocol`), so a schema/data
//! pair engineered to recurse forever cannot turn this into an unbounded
//! operation.

use serde_json::Value;

use crate::protocol::MAX_JSON_DEPTH;

/// Validate `data` against `schema`. Returns `Err(message)` describing the
/// first violation found (not an exhaustive list — sufficient for mapping
/// to a single `-32602` message).
pub fn validate(schema: &Value, data: &Value) -> Result<(), String> {
    walk(schema, data, 0)
}

fn walk(schema: &Value, data: &Value, depth: usize) -> Result<(), String> {
    if depth > MAX_JSON_DEPTH {
        return Err(format!(
            "schema validation exceeded maximum depth of {MAX_JSON_DEPTH}"
        ));
    }

    if let Some(one_of) = schema.get("oneOf").and_then(|v| v.as_array()) {
        let mut errors = Vec::new();
        for branch in one_of {
            match walk(branch, data, depth + 1) {
                Ok(()) => return Ok(()),
                Err(e) => errors.push(e),
            }
        }
        return Err(format!(
            "value did not match any oneOf branch: [{}]",
            errors.join(" | ")
        ));
    }

    if let Some(const_val) = schema.get("const")
        && const_val != data
    {
        return Err(format!("expected constant value {const_val}, got {data}"));
    }

    if let Some(enum_vals) = schema.get("enum").and_then(|v| v.as_array())
        && !enum_vals.contains(data)
    {
        return Err(format!(
            "value {data} is not one of the declared enum values"
        ));
    }

    if let Some(ty) = schema.get("type") {
        check_type(ty, data)?;
    }

    match data {
        Value::String(value) => check_string(schema, value)?,
        Value::Array(items) => check_array(schema, items, depth)?,
        Value::Number(number) => check_number(schema, number)?,
        Value::Object(map) => check_object(schema, map, depth)?,
        _ => {}
    }

    Ok(())
}

fn check_string(schema: &Value, value: &str) -> Result<(), String> {
    let length = value.chars().count() as u64;
    if let Some(minimum) = schema.get("minLength").and_then(Value::as_u64)
        && length < minimum
    {
        return Err(format!("string shorter than minLength {minimum}"));
    }
    if let Some(maximum) = schema.get("maxLength").and_then(Value::as_u64)
        && length > maximum
    {
        return Err(format!("string longer than maxLength {maximum}"));
    }
    Ok(())
}

fn check_array(schema: &Value, items: &[Value], depth: usize) -> Result<(), String> {
    let length = items.len() as u64;
    if let Some(minimum) = schema.get("minItems").and_then(Value::as_u64)
        && length < minimum
    {
        return Err(format!("array shorter than minItems {minimum}"));
    }
    if let Some(maximum) = schema.get("maxItems").and_then(Value::as_u64)
        && length > maximum
    {
        return Err(format!("array longer than maxItems {maximum}"));
    }
    if let Some(item_schema) = schema.get("items") {
        for item in items {
            walk(item_schema, item, depth + 1)?;
        }
    }
    Ok(())
}

fn check_number(schema: &Value, number: &serde_json::Number) -> Result<(), String> {
    let value = number.as_f64().unwrap_or(f64::NAN);
    if let Some(minimum) = schema.get("minimum").and_then(Value::as_f64)
        && value < minimum
    {
        return Err(format!("number less than minimum {minimum}"));
    }
    if let Some(maximum) = schema.get("maximum").and_then(Value::as_f64)
        && value > maximum
    {
        return Err(format!("number greater than maximum {maximum}"));
    }
    Ok(())
}

fn check_object(
    schema: &Value,
    map: &serde_json::Map<String, Value>,
    depth: usize,
) -> Result<(), String> {
    if let Some(required) = schema.get("required").and_then(Value::as_array) {
        for property in required {
            let key = property.as_str().unwrap_or_default();
            if !map.contains_key(key) {
                return Err(format!("missing required property: {key}"));
            }
        }
    }
    let properties = schema.get("properties").and_then(Value::as_object);
    let additional_allowed =
        !matches!(schema.get("additionalProperties"), Some(Value::Bool(false)));
    for (key, value) in map {
        match properties.and_then(|known| known.get(key)) {
            Some(sub_schema) => walk(sub_schema, value, depth + 1)?,
            None if !additional_allowed => return Err(format!("unexpected property: {key}")),
            None => {}
        }
    }
    Ok(())
}

fn check_type(ty: &Value, data: &Value) -> Result<(), String> {
    fn matches_one(t: &str, data: &Value) -> bool {
        match t {
            "object" => data.is_object(),
            "array" => data.is_array(),
            "string" => data.is_string(),
            "boolean" => data.is_boolean(),
            "null" => data.is_null(),
            "number" => data.is_number(),
            "integer" => {
                data.is_i64()
                    || data.is_u64()
                    || data
                        .as_f64()
                        .map(|f| f.fract() == 0.0 && f.is_finite())
                        .unwrap_or(false)
            }
            _ => true, // unknown type keyword: don't reject, this is a lenient subset validator
        }
    }

    let ok = match ty {
        Value::String(t) => matches_one(t, data),
        Value::Array(types) => types
            .iter()
            .any(|t| t.as_str().map(|s| matches_one(s, data)).unwrap_or(false)),
        _ => true,
    };

    if ok {
        Ok(())
    } else {
        Err(format!("type mismatch: schema requires {ty}, got {data}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn accepts_valid_object() {
        let schema = json!({
            "type": "object",
            "properties": { "smiles": { "type": "string", "minLength": 1 } },
            "required": ["smiles"],
            "additionalProperties": false
        });
        assert!(validate(&schema, &json!({ "smiles": "CCO" })).is_ok());
    }

    #[test]
    fn rejects_missing_required() {
        let schema = json!({
            "type": "object",
            "properties": { "smiles": { "type": "string" } },
            "required": ["smiles"]
        });
        assert!(validate(&schema, &json!({})).is_err());
    }

    #[test]
    fn rejects_wrong_type() {
        let schema = json!({
            "type": "object",
            "properties": { "smiles": { "type": "string" } },
            "required": ["smiles"]
        });
        assert!(validate(&schema, &json!({ "smiles": 5 })).is_err());
    }

    #[test]
    fn rejects_additional_properties() {
        let schema = json!({
            "type": "object",
            "properties": { "smiles": { "type": "string" } },
            "required": ["smiles"],
            "additionalProperties": false
        });
        assert!(validate(&schema, &json!({ "smiles": "CCO", "extra": 1 })).is_err());
    }

    #[test]
    fn enforces_min_max_items() {
        let schema = json!({
            "type": "array",
            "items": { "type": "string" },
            "minItems": 2,
            "maxItems": 3
        });
        assert!(validate(&schema, &json!(["a"])).is_err());
        assert!(validate(&schema, &json!(["a", "b"])).is_ok());
        assert!(validate(&schema, &json!(["a", "b", "c", "d"])).is_err());
    }

    #[test]
    fn one_of_matches_any_branch() {
        let schema = json!({
            "oneOf": [
                { "type": "string" },
                { "type": "integer" }
            ]
        });
        assert!(validate(&schema, &json!("x")).is_ok());
        assert!(validate(&schema, &json!(5)).is_ok());
        assert!(validate(&schema, &json!(true)).is_err());
    }

    #[test]
    fn nullable_type_union() {
        let schema = json!({ "type": ["string", "null"] });
        assert!(validate(&schema, &json!("x")).is_ok());
        assert!(validate(&schema, &Value::Null).is_ok());
        assert!(validate(&schema, &json!(5)).is_err());
    }
}
