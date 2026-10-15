use super::*;
use pyo3::types::{PyBool, PyBytes, PyInt};
fn integer(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "expected a nonnegative integer, excluding bool",
        ));
    }
    value
        .extract()
        .map_err(|_| PyValueError::new_err("integer is negative or does not fit usize"))
}
fn u32_integer(value: &Bound<'_, PyAny>) -> PyResult<u32> {
    integer(value)?
        .try_into()
        .map_err(|_| PyValueError::new_err("integer must fit u32"))
}
fn boolean(value: &Bound<'_, PyAny>) -> PyResult<bool> {
    if !value.is_instance_of::<PyBool>() {
        return Err(PyValueError::new_err("branched_paths must be bool"));
    }
    value.extract()
}
fn sequence(value: Option<&Bound<'_, PyAny>>, cap: usize) -> PyResult<Option<Vec<u32>>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let mut result = Vec::new();
    for item in value.try_iter()? {
        if result.len() >= cap {
            return Err(PyValueError::new_err(
                "integer sequence exceeds its size limit",
            ));
        }
        result.push(u32_integer(&item?)?);
    }
    Ok(Some(result))
}
fn mask(value: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<u8>>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let value = value
        .cast::<PyBytes>()
        .map_err(|_| PyValueError::new_err("set_only_bits must be packed bytes"))?;
    if value.as_bytes().len() > 8192 {
        return Err(PyValueError::new_err("mask exceeds maximum 8192 bytes"));
    }
    Ok(Some(value.as_bytes().to_vec()))
}
fn layered_error(e: chematic_fp::RdkitLayeredError) -> PyErr {
    let category = if matches!(
        e.code,
        "visit_limit" | "path_limit" | "frontier_limit" | "resource_limit" | "graph_limit"
    ) {
        "resource_limit"
    } else if matches!(
        e.code,
        "unsupported_graph"
            | "unsupported_bond"
            | "index_changing_view"
            | "unrooted_linear_unsupported"
    ) {
        "unsupported"
    } else {
        "malformed"
    };
    crate::errors::input_error(category, e.code, "rdkit_layered_config", e.detail)
}
#[pymethods]
impl Mol {
    /// Configurable experimental RDKit 2026.09.1 Layered bits and atom counts.
    /// Counts are copied and incremented (u32 wrapping), never mutated in place.
    #[pyo3(signature=(*,layer_flags=4294967295,min_path=1,max_path=7,nbits=2048,branched_paths=true,from_atoms=None,set_only_bits=None,atom_counts=None,max_visits=1000000,max_paths=100000))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_layered_config<'py>(
        &self,
        py: Python<'py>,
        #[pyo3(from_py_with=u32_integer)] layer_flags: u32,
        #[pyo3(from_py_with=integer)] min_path: usize,
        #[pyo3(from_py_with=integer)] max_path: usize,
        #[pyo3(from_py_with=integer)] nbits: usize,
        #[pyo3(from_py_with=boolean)] branched_paths: bool,
        from_atoms: Option<&Bound<'py, PyAny>>,
        set_only_bits: Option<&Bound<'py, PyAny>>,
        atom_counts: Option<&Bound<'py, PyAny>>,
        #[pyo3(from_py_with=integer)] max_visits: usize,
        #[pyo3(from_py_with=integer)] max_paths: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let config = chematic_fp::RdkitLayeredConfig {
            layer_flags,
            min_path,
            max_path,
            fp_size: nbits,
            branched_paths,
            from_atoms: sequence(from_atoms, 4096)?
                .unwrap_or_default()
                .into_iter()
                .map(chematic_core::AtomIdx)
                .collect(),
            set_only_bits: mask(set_only_bits)?,
            atom_counts: sequence(atom_counts, 65536)?,
            max_visits: max_visits as u64,
            max_paths,
        };
        let r = chematic_fp::rdkit_layered_config(&self.inner, &config).map_err(layered_error)?;
        let d = PyDict::new(py);
        d.set_item("bits", PyBytes::new(py, &r.bits))?;
        d.set_item("nbits", r.fp_size)?;
        d.set_item("atom_counts", r.atom_counts)?;
        d.set_item("selected_paths", r.selected_paths)?;
        d.set_item("search_visits", r.search_visits)?;
        Ok(d)
    }
}
