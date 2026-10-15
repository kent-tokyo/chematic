use super::*;
use pyo3::types::PyBytes;

fn property_key_error(key: &str) -> PyErr {
    pyo3::exceptions::PyKeyError::new_err(format!("'{key}' not found"))
}

fn property_value_error(key: &str, expected: &str, value: &str) -> PyErr {
    PyValueError::new_err(format!(
        "property {key:?} exists but {value:?} cannot be converted to {expected}"
    ))
}

fn parse_get_prop_options<'py>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<(bool, Option<Bound<'py, PyAny>>)> {
    if args.len() > 2 {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "GetProp() takes at most 3 arguments ({} given)",
            args.len() + 1
        )));
    }

    let mut auto_convert = false;
    let mut default = None;
    let positional_auto_convert = !args.is_empty();
    let positional_default = args.len() == 2;

    if positional_auto_convert {
        auto_convert = args.get_item(0)?.extract::<bool>()?;
    }
    if positional_default {
        default = Some(args.get_item(1)?);
    }

    if let Some(kwargs) = kwargs {
        for (raw_key, value) in kwargs.iter() {
            let key = raw_key.extract::<String>()?;
            match key.as_str() {
                "autoConvert" => {
                    if positional_auto_convert {
                        return Err(pyo3::exceptions::PyTypeError::new_err(
                            "GetProp() got multiple values for argument 'autoConvert'",
                        ));
                    }
                    auto_convert = value.extract::<bool>()?;
                }
                "default" => {
                    if positional_default {
                        return Err(pyo3::exceptions::PyTypeError::new_err(
                            "GetProp() got multiple values for argument 'default'",
                        ));
                    }
                    default = Some(value);
                }
                _ => {
                    return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                        "GetProp() got an unexpected keyword argument {key:?}"
                    )));
                }
            }
        }
    }

    Ok((auto_convert, default))
}

fn parse_typed_default<'py>(
    method: &str,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    if args.len() > 1 {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "{method}() takes at most 2 arguments ({} given)",
            args.len() + 1
        )));
    }
    let mut default = if args.is_empty() {
        None
    } else {
        Some(args.get_item(0)?)
    };
    if let Some(kwargs) = kwargs {
        for (raw_key, value) in kwargs.iter() {
            let key = raw_key.extract::<String>()?;
            if key != "default" {
                return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                    "{method}() got an unexpected keyword argument {key:?}"
                )));
            }
            if default.is_some() {
                return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                    "{method}() got multiple values for argument 'default'"
                )));
            }
            default = Some(value);
        }
    }
    Ok(default)
}

fn string_property_to_python(
    py: Python<'_>,
    value: &str,
    property_type: crate::PropertyType,
    auto_convert: bool,
) -> PyResult<Py<PyAny>> {
    if auto_convert {
        match property_type {
            crate::PropertyType::Bool => {
                let value = value == "1";
                return Ok(value
                    .into_pyobject(py)
                    .map(|bound| bound.to_owned().into_any().unbind())?);
            }
            crate::PropertyType::Double => {
                if let Ok(value) = value.parse::<f64>() {
                    return Ok(value
                        .into_pyobject(py)
                        .map(|bound| bound.into_any().unbind())?);
                }
            }
            crate::PropertyType::Int => {
                if let Ok(value) = value.parse::<i64>() {
                    return Ok(value
                        .into_pyobject(py)
                        .map(|bound| bound.into_any().unbind())?);
                }
            }
            crate::PropertyType::Unsigned => {
                if let Ok(value) = value.parse::<u64>() {
                    return Ok(value
                        .into_pyobject(py)
                        .map(|bound| bound.into_any().unbind())?);
                }
            }
            crate::PropertyType::String => {
                if let Ok(value) = value.parse::<i64>() {
                    return Ok(value
                        .into_pyobject(py)
                        .map(|bound| bound.into_any().unbind())?);
                }
                if let Ok(value) = value.parse::<f64>() {
                    return Ok(value
                        .into_pyobject(py)
                        .map(|bound| bound.into_any().unbind())?);
                }
            }
        }
    }
    Ok(value
        .into_pyobject(py)
        .map(|bound| bound.into_any().unbind())?)
}

#[pymethods]
impl Mol {
    // -----------------------------------------------------------------------
    // SD properties (RDKit-compatible)
    // -----------------------------------------------------------------------

    /// Get an SD property by name.
    ///
    /// This follows RDKit 2026.09.1's ``GetProp(key, autoConvert=False,
    /// default=...)`` contract.  When ``default`` is omitted a missing key
    /// still raises ``KeyError``; an explicitly supplied fallback, including
    /// ``None``, is returned unchanged.
    #[pyo3(name = "GetProp")]
    #[pyo3(signature = (key, *args, **kwargs))]
    #[pyo3(text_signature = "($self, key, autoConvert=False, default=...)")]
    fn get_prop<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        let (auto_convert, default) = parse_get_prop_options(args, kwargs)?;
        if let Some(value) = self.props.get(key) {
            string_property_to_python(py, value, self.props.property_type(key), auto_convert)
        } else if let Some(default) = default {
            Ok(default.unbind())
        } else {
            Err(property_key_error(key))
        }
    }

    /// Set an SD property.
    #[pyo3(name = "SetProp")]
    fn set_prop(&mut self, key: String, value: String) {
        self.props
            .insert_typed(key, value, crate::PropertyType::String);
    }

    /// Return ``True`` if the property exists.
    #[pyo3(name = "HasProp")]
    fn has_prop(&self, key: &str) -> bool {
        self.props.contains_key(key)
    }

    /// Return all SD properties as a Python dict.
    #[pyo3(name = "GetPropsAsDict")]
    fn get_props_as_dict<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        let d = PyDict::new(py);
        for (k, v) in &self.props {
            d.set_item(k, v).ok();
        }
        d
    }

    /// Return a list of all property names.
    #[pyo3(name = "GetPropNames")]
    fn get_prop_names(&self) -> Vec<String> {
        self.props.keys().cloned().collect()
    }

    /// Remove a property by name (no-op if not present).
    #[pyo3(name = "ClearProp")]
    fn clear_prop(&mut self, key: &str) {
        self.props.clear(key);
    }

    /// Return the molecule title stored in the conventional ``_Name`` property.
    ///
    /// RDKit 2026.09.1 returns an empty string when no name has been set instead
    /// of raising ``KeyError``.  ``GetName`` and the generic property API share
    /// the same storage, so updates through either surface remain visible.
    #[pyo3(name = "GetName")]
    fn get_name(&self) -> String {
        self.props.get("_Name").cloned().unwrap_or_default()
    }

    /// Set the molecule title in the conventional ``_Name`` property.
    #[pyo3(name = "SetName")]
    fn set_name(&mut self, name: String) {
        self.props
            .insert_typed("_Name".to_string(), name, crate::PropertyType::String);
    }

    /// Pickle support with a versioned, bounded molecule-state payload.
    ///
    /// Python's pickle format is unsafe for untrusted data. The chematic
    /// payload decoder nevertheless validates lengths, indices, enum values,
    /// UTF-8, and trailing bytes before rebuilding a molecule.
    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Py<PyAny>, (Py<PyBytes>,))> {
        let restore = py
            .import("chematic")?
            .getattr("_restore_mol_pickle")?
            .unbind();
        let payload = crate::pickle::encode(self).map_err(PyValueError::new_err)?;
        Ok((restore, (PyBytes::new(py, &payload).unbind(),)))
    }

    /// Set an integer property (stored as its string representation).
    #[pyo3(name = "SetIntProp")]
    fn set_int_prop(&mut self, key: String, val: i64) {
        self.props
            .insert_typed(key, val.to_string(), crate::PropertyType::Int);
    }

    /// Get an integer property, or ``default`` when it is absent.
    #[pyo3(name = "GetIntProp", signature = (key, *args, **kwargs))]
    #[pyo3(text_signature = "($self, key, default=...)")]
    fn get_int_prop<'py>(
        &self,
        key: &str,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<i64> {
        let default = parse_typed_default("GetIntProp", args, kwargs)?;
        match self.props.get(key) {
            Some(value)
                if matches!(
                    self.props.property_type(key),
                    crate::PropertyType::String
                        | crate::PropertyType::Int
                        | crate::PropertyType::Unsigned
                ) =>
            {
                value
                    .parse::<i64>()
                    .map_err(|_| property_value_error(key, "an integer", value))
            }
            Some(value) => Err(property_value_error(key, "an integer", value)),
            None => match default {
                Some(default) => default.extract::<i64>(),
                None => Err(property_key_error(key)),
            },
        }
    }

    /// Set an unsigned integer property (stored as its string representation).
    #[pyo3(name = "SetUnsignedProp")]
    fn set_unsigned_prop(&mut self, key: String, val: u64) {
        self.props
            .insert_typed(key, val.to_string(), crate::PropertyType::Unsigned);
    }

    /// Get an unsigned integer property, or ``default`` when it is absent.
    #[pyo3(name = "GetUnsignedProp", signature = (key, *args, **kwargs))]
    #[pyo3(text_signature = "($self, key, default=...)")]
    fn get_unsigned_prop<'py>(
        &self,
        key: &str,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<u64> {
        let default = parse_typed_default("GetUnsignedProp", args, kwargs)?;
        match self.props.get(key) {
            Some(value)
                if matches!(
                    self.props.property_type(key),
                    crate::PropertyType::String
                        | crate::PropertyType::Int
                        | crate::PropertyType::Unsigned
                ) =>
            {
                value
                    .parse::<u64>()
                    .map_err(|_| property_value_error(key, "an unsigned integer", value))
            }
            Some(value) => Err(property_value_error(key, "an unsigned integer", value)),
            None => match default {
                Some(default) => default.extract::<u64>(),
                None => Err(property_key_error(key)),
            },
        }
    }

    /// Set a float property (stored as its string representation).
    #[pyo3(name = "SetDoubleProp")]
    fn set_double_prop(&mut self, key: String, val: f64) {
        self.props
            .insert_typed(key, val.to_string(), crate::PropertyType::Double);
    }

    /// Get a floating-point property, or ``default`` when it is absent.
    #[pyo3(name = "GetDoubleProp", signature = (key, *args, **kwargs))]
    #[pyo3(text_signature = "($self, key, default=...)")]
    fn get_double_prop<'py>(
        &self,
        key: &str,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<f64> {
        let default = parse_typed_default("GetDoubleProp", args, kwargs)?;
        match self.props.get(key) {
            Some(value)
                if matches!(
                    self.props.property_type(key),
                    crate::PropertyType::String | crate::PropertyType::Double
                ) =>
            {
                value
                    .parse::<f64>()
                    .map_err(|_| property_value_error(key, "a float", value))
            }
            Some(value) => Err(property_value_error(key, "a float", value)),
            None => match default {
                Some(default) => default.extract::<f64>(),
                None => Err(property_key_error(key)),
            },
        }
    }

    /// Set a boolean property (stored as ``"1"`` / ``"0"``).
    #[pyo3(name = "SetBoolProp")]
    fn set_bool_prop(&mut self, key: String, val: bool) {
        self.props.insert_typed(
            key,
            if val { "1" } else { "0" }.to_string(),
            crate::PropertyType::Bool,
        );
    }

    /// Get a boolean property, or ``default`` when it is absent.
    #[pyo3(name = "GetBoolProp", signature = (key, *args, **kwargs))]
    #[pyo3(text_signature = "($self, key, default=...)")]
    fn get_bool_prop<'py>(
        &self,
        key: &str,
        args: &Bound<'py, PyTuple>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<bool> {
        let default = parse_typed_default("GetBoolProp", args, kwargs)?;
        match (
            self.props.get(key).map(String::as_str),
            self.props.property_type(key),
        ) {
            (Some("1"), crate::PropertyType::String | crate::PropertyType::Bool) => Ok(true),
            (Some("0"), crate::PropertyType::String | crate::PropertyType::Bool) => Ok(false),
            (Some(value), _) => Err(property_value_error(key, "a boolean", value)),
            (None, _) => match default {
                Some(default) if default.is_none() => Ok(false),
                Some(default) => default.extract::<bool>(),
                None => Err(property_key_error(key)),
            },
        }
    }

    // -----------------------------------------------------------------------
    // Dunder methods
    // -----------------------------------------------------------------------

    fn __repr__(&self) -> String {
        format!("Mol('{}')", chematic_smiles::canonical_smiles(&self.inner))
    }

    fn __str__(&self) -> String {
        chematic_smiles::canonical_smiles(&self.inner)
    }
}
