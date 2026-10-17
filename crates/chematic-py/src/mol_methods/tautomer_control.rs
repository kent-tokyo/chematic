//! Python callback exceptions leave the enumerator through its normal cancel path.
use super::*;
use pyo3::exceptions::PyTypeError;
use pyo3::types::PyBool;

type Output = (Vec<String>, Vec<String>, String);

#[pymethods]
impl Mol {
    /// Internal bridge for the typed `enumerate_tautomers_checked` API.
    fn _rdkit_tautomers_controlled(
        &self,
        py: Python<'_>,
        max_tautomers: usize,
        max_transforms: usize,
        callback: Bound<'_, PyAny>,
    ) -> PyResult<Output> {
        if !callback.is_callable() {
            return Err(PyTypeError::new_err("callback must be callable"));
        }
        let mut callback_error = None;
        let mut checkpoint = |progress: &chematic_smiles::RdkitTautomerProgress| {
            let call = (|| {
                let keys = PyTuple::new(py, &progress.keys)?;
                let answer = callback.call1((progress.transforms, keys))?;
                if !answer.is_instance_of::<PyBool>() {
                    return Err(PyTypeError::new_err("callback must return bool"));
                }
                answer.extract::<bool>()
            })();
            match call {
                Ok(keep_going) => keep_going,
                Err(error) => {
                    callback_error = Some(error);
                    false
                }
            }
        };
        let result = chematic_smiles::rdkit_enumerate_tautomers_with_callback(
            &self.inner,
            max_tautomers,
            max_transforms,
            &mut checkpoint,
        );
        if let Some(error) = callback_error {
            return Err(error);
        }
        let result = result.map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok((result.smiles, result.keys, format!("{:?}", result.status)))
    }
}
