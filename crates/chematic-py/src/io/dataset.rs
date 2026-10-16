use super::*;

/// Indexed random access to a seekable SDF file.
///
/// Construction scans only record boundaries and title lines. Molecular
/// graphs and SD fields are parsed lazily when a record is requested, so a
/// large file does not need to be materialized as Python objects.
#[pyclass(name = "SdfDataset", frozen)]
pub struct PySdfDataset {
    inner: Arc<chematic_mol::SdfDataset>,
}

#[pymethods]
impl PySdfDataset {
    #[new]
    fn new(
        #[pyo3(from_py_with=crate::filesystem::path)] path: std::path::PathBuf,
    ) -> PyResult<Self> {
        Self::open(path)
    }

    #[staticmethod]
    fn open(
        #[pyo3(from_py_with=crate::filesystem::path)] path: std::path::PathBuf,
    ) -> PyResult<Self> {
        chematic_mol::SdfDataset::open(path)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(sdf_dataset_error)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        key: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if key.is_instance_of::<pyo3::types::PyInt>()
            && !key.is_instance_of::<pyo3::types::PyBool>()
        {
            let index = super::dataset_selector::scalar_index(key, self.inner.len())?;
            return Py::new(py, self.get(index as isize)?).map(|v| v.into_bound(py).into_any());
        }
        Py::new(py, self.select(py, key, "raise", 64 << 20)?).map(|v| v.into_bound(py).into_any())
    }

    /// Ordered selection; slices, integer lists/tuples or full-length bool masks.
    /// Malformed records retain their positions with errors="keep".
    #[pyo3(signature=(key,*,errors="keep",max_bytes=67108864))]
    fn select(
        &self,
        py: Python<'_>,
        key: &Bound<'_, PyAny>,
        errors: &str,
        #[pyo3(from_py_with=super::dataset_selector::positive_integer)] max_bytes: usize,
    ) -> PyResult<super::dataset_batch::PySdfSelectionBatch> {
        let raise = super::dataset_selector::error_policy(errors)?;
        let indices = super::dataset_selector::indices(key, self.inner.len(), 10_000)?;
        super::dataset_batch::read_batch(py, &self.inner, &indices, raise, max_bytes)
    }

    /// Pull source-ordered batches while keeping rejected record positions.
    /// The iterator owns the indexed file handle and survives this wrapper.
    #[pyo3(signature=(size=1024,*,errors="keep",max_bytes=67108864))]
    fn batches(
        &self,
        #[pyo3(from_py_with=super::dataset_selector::positive_integer)] size: usize,
        errors: &str,
        #[pyo3(from_py_with=super::dataset_selector::positive_integer)] max_bytes: usize,
    ) -> PyResult<super::dataset_batch::PySdfDatasetBatches> {
        super::dataset_batch::PySdfDatasetBatches::new(
            Arc::clone(&self.inner),
            size,
            errors,
            max_bytes,
        )
    }

    /// Parse one record. Negative indices follow normal Python semantics.
    fn get(&self, index: isize) -> PyResult<PySdfRecord> {
        let index = dataset_index(index, self.inner.len())?;
        self.inner
            .get(index)
            .map_err(sdf_dataset_error)?
            .map(py_sdf_record)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("SDF record index out of range"))
    }

    /// Read one raw record without the terminating ``$$$$`` line.
    fn raw_record(&self, index: isize) -> PyResult<String> {
        let index = dataset_index(index, self.inner.len())?;
        self.inner
            .record_text(index)
            .map_err(sdf_dataset_error)?
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("SDF record index out of range"))
    }

    /// Return byte-range metadata without parsing the molecular graph.
    fn metadata(&self, index: isize) -> PyResult<PySdfRecordMetadata> {
        let index = dataset_index(index, self.inner.len())?;
        self.inner
            .metadata(index)
            .map(PySdfRecordMetadata::from)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("SDF record index out of range"))
    }

    /// Parse a half-open record range. The stop endpoint is clamped to
    /// ``len(dataset)``; malformed records raise ``ChematicInputError``.
    #[pyo3(signature = (start, stop=None))]
    fn read_range(&self, start: usize, stop: Option<usize>) -> PyResult<Vec<PySdfRecord>> {
        let stop = stop.unwrap_or(self.inner.len());
        self.inner
            .read_range(start, stop)
            .map_err(sdf_dataset_error)?
            .into_iter()
            .map(|record| record.map(py_sdf_record).map_err(sdf_dataset_error))
            .collect()
    }
}
