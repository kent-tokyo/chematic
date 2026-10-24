//! Decode explicit ExtXYZ kinds without coercing integers, strings or booleans.
use super::*;
use chematic_mol::{XyzProperty, XyzPropertyKind as Kind, XyzValue as Value};
use pyo3::types::{PyBool, PyInt};

pub(super) fn kind_tag(kind: Kind) -> &'static str {
    match kind {
        Kind::Real => "R",
        Kind::Integer => "I",
        Kind::String => "S",
        Kind::Logical => "L",
    }
}

fn schema_entry(entry: &Bound<'_, PyAny>) -> PyResult<(Kind, usize)> {
    let (tag, count): (String, Bound<'_, PyAny>) = entry
        .extract()
        .map_err(|_| PyValueError::new_err("schema entry must be (kind, component_count)"))?;
    let kind = match tag.as_str() {
        "R" => Kind::Real,
        "I" => Kind::Integer,
        "S" => Kind::String,
        "L" => Kind::Logical,
        _ => return Err(PyValueError::new_err("schema kind must be R, I, S or L")),
    };
    if !count.is_instance_of::<PyInt>() || count.is_instance_of::<PyBool>() {
        return Err(PyValueError::new_err(
            "component_count must be a positive integer",
        ));
    }
    let count: usize = count
        .extract()
        .map_err(|_| PyValueError::new_err("component_count must fit usize"))?;
    if count == 0 {
        return Err(PyValueError::new_err("component_count must be positive"));
    }
    Ok((kind, count))
}

fn scalar(kind: Kind, value: &Bound<'_, PyAny>) -> PyResult<Value> {
    let invalid = || PyValueError::new_err(format!("expected {} property value", kind_tag(kind)));
    match kind {
        Kind::Real => {
            if value.is_instance_of::<PyBool>() {
                return Err(invalid());
            }
            let real: f64 = value.extract().map_err(|_| invalid())?;
            if !real.is_finite() {
                return Err(invalid());
            }
            if value.hasattr("__index__")? {
                let integer = value.call_method0("__index__").map_err(|_| invalid())?;
                if !integer.is_instance_of::<PyInt>() || integer.is_instance_of::<PyBool>() {
                    return Err(invalid());
                }
                if !integer.eq(real)? {
                    return Err(PyValueError::new_err(
                        "integer is not exactly representable as R; declare an I column",
                    ));
                }
            }
            Ok(Value::Real(real))
        }
        Kind::Integer => {
            if !value.is_instance_of::<PyInt>() || value.is_instance_of::<PyBool>() {
                return Err(invalid());
            }
            value.extract().map(Value::Integer).map_err(|_| invalid())
        }
        Kind::String => value.extract().map(Value::Str).map_err(|_| invalid()),
        Kind::Logical => {
            if !value.is_instance_of::<PyBool>() {
                return Err(invalid());
            }
            value.extract().map(Value::Logical).map_err(|_| invalid())
        }
    }
}

fn column(
    name: String,
    rows: &Bound<'_, PyAny>,
    schema: Option<&Bound<'_, PyAny>>,
    atoms: usize,
) -> PyResult<XyzProperty> {
    let rows: Vec<Vec<Bound<'_, PyAny>>> = rows
        .extract()
        .map_err(|_| PyValueError::new_err(format!("properties[{name:?}] must contain rows")))?;
    if rows.len() != atoms {
        return Err(PyValueError::new_err(format!(
            "properties[{name:?}] has {} rows; expected {atoms}",
            rows.len()
        )));
    }
    let (kind, count) = match schema {
        Some(entry) => schema_entry(entry)?,
        None => (Kind::Real, rows.first().map_or(0, Vec::len)),
    };
    let mut values = Vec::with_capacity(rows.len());
    for row in rows {
        if row.len() != count {
            return Err(PyValueError::new_err(format!(
                "properties[{name:?}] row width differs from {count}"
            )));
        }
        values.push(
            row.iter()
                .map(|v| scalar(kind, v))
                .collect::<PyResult<_>>()?,
        );
    }
    Ok(XyzProperty {
        name,
        kind,
        count,
        values,
    })
}

pub(super) fn decode(
    properties: Option<&Bound<'_, PyDict>>,
    schema: Option<&Bound<'_, PyDict>>,
    atoms: usize,
) -> PyResult<Vec<XyzProperty>> {
    if let Some(schema) = schema {
        let mismatch =
            || PyValueError::new_err("property_schema names must exactly match properties");
        if schema.len() != properties.map_or(0, |p| p.len()) {
            return Err(mismatch());
        }
        for key in schema.keys() {
            let name: String = key
                .extract()
                .map_err(|_| PyValueError::new_err("property names must be strings"))?;
            if let Some(properties) = properties
                && !properties.contains(&name)?
            {
                return Err(mismatch());
            }
        }
    }
    let mut result = Vec::new();
    if let Some(properties) = properties {
        for (key, rows) in properties.iter() {
            let name: String = key
                .extract()
                .map_err(|_| PyValueError::new_err("property names must be strings"))?;
            let entry = schema.map(|s| s.get_item(&key)).transpose()?.flatten();
            result.push(column(name, &rows, entry.as_ref(), atoms)?);
        }
    }
    Ok(result)
}
