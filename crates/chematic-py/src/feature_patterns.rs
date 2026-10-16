use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyStringMethods};

/// Owned custom SMARTS feature definitions (up to 32), in bit order.
/// Each matching atom receives its pattern's bit. Empty sets mean zero masks.
/// Fingerprints use native FNV circular hashing, not RDKit Morgan bit IDs.
#[pyclass(name = "FeaturePatternSet", frozen)]
pub(crate) struct PyFeaturePatternSet {
    inner: chematic_fp::feature_patterns::FeaturePatternSet,
}

#[pymethods]
impl PyFeaturePatternSet {
    #[new]
    fn new(patterns: &Bound<'_, PyAny>) -> PyResult<Self> {
        if patterns.is_instance_of::<pyo3::types::PyString>()
            || patterns.is_instance_of::<pyo3::types::PyBytes>()
        {
            return Err(PyTypeError::new_err(
                "patterns must be an iterable of SMARTS strings",
            ));
        }
        let mut owned = Vec::new();
        for item in patterns.try_iter()? {
            let item = item?;
            if owned.len() == 32 {
                return Err(PyValueError::new_err(
                    "at most 32 feature patterns are supported",
                ));
            }
            let text = item.cast::<pyo3::types::PyString>()?.to_str()?;
            if text.len() > 4096 {
                return Err(PyValueError::new_err("feature pattern exceeds 4096 bytes"));
            }
            owned.push(text.to_owned());
        }
        Ok(Self {
            inner: chematic_fp::feature_patterns::FeaturePatternSet::new(&owned)
                .map_err(|error| PyValueError::new_err(error.to_string()))?,
        })
    }

    #[getter]
    fn patterns(&self) -> Vec<String> {
        self.inner.patterns().to_vec()
    }

    /// Per-source-atom u32 masks. No partial result on matching budget exhaustion.
    #[pyo3(signature = (mol, *, max_visits=100_000))]
    fn invariants(&self, mol: &crate::Mol, max_visits: u64) -> PyResult<Vec<u32>> {
        self.inner
            .invariants(&mol.inner, max_visits)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Packed LSB-first native circular feature bits. This is not RDKit-bit-exact.
    #[pyo3(signature = (mol, *, radius=2, nbits=2048, max_visits=100_000))]
    fn fingerprint<'py>(
        &self,
        py: Python<'py>,
        mol: &crate::Mol,
        radius: u32,
        nbits: usize,
        max_visits: u64,
    ) -> PyResult<Bound<'py, PyBytes>> {
        self.inner
            .fingerprint(&mol.inner, radius, nbits, max_visits)
            .map(|bits| PyBytes::new(py, &bits.to_le_bytes()))
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}
