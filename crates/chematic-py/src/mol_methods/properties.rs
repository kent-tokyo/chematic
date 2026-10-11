use super::*;

#[pymethods]
impl Mol {
    // -----------------------------------------------------------------------
    // SD properties (RDKit-compatible)
    // -----------------------------------------------------------------------

    /// Get an SD property by name. Raises ``KeyError`` if not present.
    #[pyo3(name = "GetProp")]
    fn get_prop(&self, key: &str) -> PyResult<String> {
        self.props
            .get(key)
            .cloned()
            .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("'{key}' not found")))
    }

    /// Set an SD property.
    #[pyo3(name = "SetProp")]
    fn set_prop(&mut self, key: String, value: String) {
        self.props.insert(key, value);
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
        self.props.remove(key);
    }

    /// Set an integer property (stored as its string representation).
    #[pyo3(name = "SetIntProp")]
    fn set_int_prop(&mut self, key: String, val: i64) {
        self.props.insert(key, val.to_string());
    }

    /// Set a float property (stored as its string representation).
    #[pyo3(name = "SetDoubleProp")]
    fn set_double_prop(&mut self, key: String, val: f64) {
        self.props.insert(key, val.to_string());
    }

    /// Set a boolean property (stored as ``"1"`` / ``"0"``).
    #[pyo3(name = "SetBoolProp")]
    fn set_bool_prop(&mut self, key: String, val: bool) {
        self.props
            .insert(key, if val { "1" } else { "0" }.to_string());
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
