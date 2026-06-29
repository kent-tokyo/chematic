use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

use crate::Mol;

// ---------------------------------------------------------------------------
// PySdfRecord — one entry from an SDF file
// ---------------------------------------------------------------------------

/// One record from an SDF file: a molecule, its name, and SD data properties.
///
///     for rec in chematic.iter_sdf("compounds.sdf"):
///         print(rec.smiles, rec.name)
///         print(rec.properties())          # dict of SD fields
///         activity = rec.get("Activity")   # None if not present
#[pyclass(name = "SdfRecord")]
pub struct PySdfRecord {
    #[pyo3(get)]
    pub(crate) mol: Mol,
    #[pyo3(get)]
    pub name: String,
    pub props: std::collections::HashMap<String, String>,
}

#[pymethods]
impl PySdfRecord {
    /// Canonical SMILES of the molecule.
    #[getter]
    fn smiles(&self) -> String {
        chematic_smiles::canonical_smiles(&self.mol.inner)
    }

    /// Return all SD data fields as a Python dict.
    fn properties<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        let d = PyDict::new(py);
        for (k, v) in &self.props {
            d.set_item(k, v).ok();
        }
        d
    }

    /// Get one SD property by name, returning ``None`` if not present.
    fn get(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(|s| s.as_str())
    }

    fn __repr__(&self) -> String {
        format!(
            "SdfRecord(smiles='{}', name={:?}, n_props={})",
            chematic_smiles::canonical_smiles(&self.mol.inner),
            self.name,
            self.props.len()
        )
    }
}

// ---------------------------------------------------------------------------
// SdfIterStr — lazy Python iterator over pre-collected SDF records
// ---------------------------------------------------------------------------

/// Iterator returned by ``iter_sdf`` and ``iter_sdf_str``.
///
/// Records are parsed eagerly on creation and yielded lazily.
/// Invalid records are silently skipped.
#[pyclass]
pub struct SdfIter {
    records: std::vec::IntoIter<PySdfRecord>,
}

#[pymethods]
impl SdfIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<PySdfRecord> {
        self.records.next()
    }

    fn __len__(&self) -> usize {
        self.records.len()
    }
}

fn parse_to_iter(content: &str) -> SdfIter {
    let records: Vec<PySdfRecord> = chematic_mol::SdfRecordReader::new(content)
        .filter_map(|r| r.ok())
        .map(|rec| PySdfRecord {
            mol: Mol {
                inner: Arc::new(rec.mol),
            },
            name: rec.meta.name,
            props: rec.properties,
        })
        .collect();
    SdfIter {
        records: records.into_iter(),
    }
}

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------

/// Iterate over molecules in an SDF file by file path (streaming).
///
/// Reads one MOL block at a time — suitable for large SDF files where loading
/// the entire content into memory is undesirable.  Invalid records are silently
/// skipped.  IO errors raise ``IOError``.
///
///     for rec in chematic.iter_sdf("chembl.sdf"):
///         print(rec.smiles, rec.get("pChEMBL Value"))
///
///     # Build a Pandas DataFrame:
///     import pandas as pd
///     rows = [{"smiles": r.smiles, **r.properties()} for r in chematic.iter_sdf("data.sdf")]
///     df = pd.DataFrame(rows)
///
/// .. note::
///     Unlike previous versions, this iterator does not support ``len()``
///     because the total record count is not known until the file is fully read.
///     Use ``iter_sdf_str`` if you need ``len()`` on an in-memory string.
#[pyfunction]
pub fn iter_sdf(path: &str) -> PyResult<SdfFileIter> {
    let file = std::fs::File::open(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
    Ok(SdfFileIter {
        inner: chematic_mol::SdfFileReader::new(std::io::BufReader::new(file)),
    })
}

/// Iterate over batches of SDF records from a file (streaming).
///
/// Yields lists of up to ``batch_size`` SDF records.  Suitable for pipelining
/// SDF parsing with bulk computations:
///
///     for batch in chematic.iter_sdf_batched("large.sdf", batch_size=1000):
///         smiles = [r.smiles for r in batch]
///         descs = chematic.bulk.descriptors(smiles)
#[pyfunction]
#[pyo3(signature = (path, batch_size=1000))]
pub fn iter_sdf_batched(path: &str, batch_size: usize) -> PyResult<SdfBatchIter> {
    let file = std::fs::File::open(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
    Ok(SdfBatchIter {
        inner: chematic_mol::SdfFileReader::new(std::io::BufReader::new(file)),
        batch_size,
    })
}

/// Iterate over molecules in an SDF string.
///
///     sdf_text = open("compounds.sdf").read()
///     for rec in chematic.iter_sdf_str(sdf_text):
///         print(rec.smiles, rec.properties())
#[pyfunction]
pub fn iter_sdf_str(content: &str) -> SdfIter {
    parse_to_iter(content)
}

// ---------------------------------------------------------------------------
// SdfFileIter — true streaming iterator backed by a file
// ---------------------------------------------------------------------------

/// Streaming SDF iterator that reads one record at a time from disk.
///
/// Unlike ``SdfIter`` (which loads the entire file into memory), this iterator
/// reads one MOL block at a time and is suitable for large SDF files.
/// Invalid records are silently skipped.  IO errors are raised as Python
/// ``IOError``.
///
/// Does **not** support ``len()`` because the total record count is unknown
/// until the file is fully read.
#[pyclass]
pub struct SdfFileIter {
    inner: chematic_mol::SdfFileReader<std::io::BufReader<std::fs::File>>,
}

#[pymethods]
impl SdfFileIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<PySdfRecord>> {
        loop {
            match self.inner.next() {
                None => return Ok(None),
                Some(Ok(rec)) => {
                    return Ok(Some(PySdfRecord {
                        mol: Mol {
                            inner: Arc::new(rec.mol),
                        },
                        name: rec.meta.name,
                        props: rec.properties,
                    }));
                }
                Some(Err(chematic_mol::MolParseError::Io(msg))) => {
                    return Err(pyo3::exceptions::PyIOError::new_err(msg));
                }
                Some(Err(_)) => {
                    // Parse error → skip malformed record, continue
                    continue;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SdfBatchIter — streaming batch iterator
// ---------------------------------------------------------------------------

/// Streaming SDF iterator that yields batches (lists) of SDF records.
///
/// Useful for pipelining SDF parsing with bulk descriptor computation:
///
///     for batch in chematic.iter_sdf_batched("large.sdf", batch_size=1000):
///         smiles = [r.smiles for r in batch]
///         descs = chematic.bulk.descriptors(smiles)
#[pyclass]
pub struct SdfBatchIter {
    inner: chematic_mol::SdfFileReader<std::io::BufReader<std::fs::File>>,
    batch_size: usize,
}

#[pymethods]
impl SdfBatchIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<Vec<PySdfRecord>>> {
        let mut batch = Vec::with_capacity(self.batch_size);
        loop {
            if batch.len() >= self.batch_size {
                break;
            }
            match self.inner.next() {
                None => break,
                Some(Ok(rec)) => {
                    batch.push(PySdfRecord {
                        mol: Mol {
                            inner: Arc::new(rec.mol),
                        },
                        name: rec.meta.name,
                        props: rec.properties,
                    });
                }
                Some(Err(chematic_mol::MolParseError::Io(msg))) => {
                    return Err(pyo3::exceptions::PyIOError::new_err(msg));
                }
                Some(Err(_)) => continue, // skip malformed
            }
        }
        if batch.is_empty() {
            Ok(None)
        } else {
            Ok(Some(batch))
        }
    }
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySdfRecord>()?;
    m.add_class::<SdfIter>()?;
    m.add_class::<SdfFileIter>()?;
    m.add_class::<SdfBatchIter>()?;
    m.add_function(wrap_pyfunction!(iter_sdf, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_str, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_batched, m)?)?;
    Ok(())
}
