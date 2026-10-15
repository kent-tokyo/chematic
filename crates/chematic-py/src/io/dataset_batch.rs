use super::*;

struct Failure {
    position: usize,
    index: usize,
    code: &'static str,
    message: String,
}

/// Selected source positions, including placeholders and diagnostics for failures.
#[pyclass(name = "SdfSelectionBatch", frozen)]
pub(super) struct PySdfSelectionBatch {
    indices: Vec<usize>,
    records: Vec<Option<Py<PySdfRecord>>>,
    failures: Vec<Failure>,
}
#[pymethods]
impl PySdfSelectionBatch {
    fn __len__(&self) -> usize {
        self.records.len()
    }
    fn source_indices(&self) -> Vec<usize> {
        self.indices.clone()
    }
    fn valid_mask(&self) -> Vec<bool> {
        self.records.iter().map(Option::is_some).collect()
    }
    fn records(&self, py: Python<'_>) -> Vec<Option<Py<PySdfRecord>>> {
        self.records
            .iter()
            .map(|r| r.as_ref().map(|r| r.clone_ref(py)))
            .collect()
    }
    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Option<Py<PySdfRecord>>> {
        let index = dataset_index(index, self.records.len())?;
        Ok(self.records[index].as_ref().map(|r| r.clone_ref(py)))
    }
    fn to_smiles_list(&self, py: Python<'_>) -> Vec<Option<String>> {
        self.records
            .iter()
            .map(|r| r.as_ref().map(|r| r.borrow(py).smiles()))
            .collect()
    }
    fn errors<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        self.failures
            .iter()
            .map(|e| {
                let d = PyDict::new(py);
                d.set_item("position", e.position)?;
                d.set_item("source_index", e.index)?;
                d.set_item("code", e.code)?;
                d.set_item("message", &e.message)?;
                Ok(d)
            })
            .collect()
    }
}

fn selection_error(error: chematic_mol::SdfSelectionError) -> PyErr {
    use chematic_mol::SdfSelectionError;
    match error {
        SdfSelectionError::Index { .. } => {
            pyo3::exceptions::PyIndexError::new_err(error.to_string())
        }
        SdfSelectionError::ResourceLimit {
            resource,
            actual,
            limit,
        } => sdf_dataset_error(chematic_mol::MolParseError::ResourceLimit {
            resource,
            actual,
            limit,
        }),
        SdfSelectionError::Read(error) => sdf_dataset_error(error),
    }
}
pub(super) fn read_batch(
    py: Python<'_>,
    dataset: &chematic_mol::SdfDataset,
    indices: &[usize],
    raise: bool,
    max_bytes: usize,
) -> PyResult<PySdfSelectionBatch> {
    let limits = chematic_mol::SdfSelectionLimits {
        max_records: 10_000,
        max_total_record_bytes: max_bytes,
    };
    let batch = dataset
        .read_indices_with_limits(indices, &limits)
        .map_err(selection_error)?;
    let mut records = Vec::with_capacity(batch.len());
    let mut failures = Vec::new();
    for (position, record) in batch.records.into_iter().enumerate() {
        match record {
            Ok(record) => records.push(Some(Py::new(py, py_sdf_record(record))?)),
            Err(error) => {
                if raise {
                    return Err(sdf_dataset_error(error));
                }
                let code = if matches!(error, chematic_mol::MolParseError::ResourceLimit { .. }) {
                    "resource_limit"
                } else {
                    "malformed"
                };
                failures.push(Failure {
                    position,
                    index: indices[position],
                    code,
                    message: error.to_string(),
                });
                records.push(None);
            }
        }
    }
    Ok(PySdfSelectionBatch {
        indices: batch.indices,
        records,
        failures,
    })
}

/// Lazy indexed-file batches; source indices are created only for the next chunk.
#[pyclass(name = "SdfDatasetBatches")]
pub(super) struct PySdfDatasetBatches {
    dataset: Arc<chematic_mol::SdfDataset>,
    next: usize,
    size: usize,
    raise: bool,
    max_bytes: usize,
}
impl PySdfDatasetBatches {
    pub(super) fn new(
        dataset: Arc<chematic_mol::SdfDataset>,
        size: usize,
        errors: &str,
        max_bytes: usize,
    ) -> PyResult<Self> {
        if size == 0 || size > 10_000 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "size must be 1..=10000",
            ));
        }
        let raise = super::dataset_selector::error_policy(errors)?;
        Ok(Self {
            dataset,
            next: 0,
            size,
            raise,
            max_bytes,
        })
    }
}
#[pymethods]
impl PySdfDatasetBatches {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __len__(&self) -> usize {
        (self.dataset.len() - self.next).div_ceil(self.size)
    }
    fn cancel(&mut self) {
        self.next = self.dataset.len();
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PySdfSelectionBatch>> {
        if self.next >= self.dataset.len() {
            return Ok(None);
        }
        let end = self.next.saturating_add(self.size).min(self.dataset.len());
        let indices: Vec<_> = (self.next..end).collect();
        self.next = end;
        match read_batch(py, &self.dataset, &indices, self.raise, self.max_bytes) {
            Ok(batch) => Ok(Some(batch)),
            Err(error) => {
                self.cancel();
                Err(error)
            }
        }
    }
}
pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySdfSelectionBatch>()?;
    m.add_class::<PySdfDatasetBatches>()?;
    Ok(())
}
