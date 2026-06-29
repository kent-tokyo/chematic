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
                props: Default::default(),
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
                            props: Default::default(),
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
                            props: Default::default(),
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
// SDMolSupplier — RDKit-compatible streaming SDF reader
// ---------------------------------------------------------------------------

/// Streaming SDF reader that yields ``Mol`` objects with SD properties attached.
///
/// RDKit-compatible API:
///
///     sup = chematic.SDMolSupplier("compounds.sdf")
///     for mol in sup:
///         if mol is None:
///             continue          # malformed record (strict_parsing=False only)
///         print(mol.smiles, mol.GetProp("Activity"))
///
/// ``sanitize`` and ``remove_hs`` are accepted for API compatibility but are no-ops:
/// chematic's parser already applies default aromaticity perception and stores
/// only heavy atoms.
#[pyclass(name = "SDMolSupplier")]
pub struct SdMolSupplier {
    inner: chematic_mol::SdfFileReader<std::io::BufReader<std::fs::File>>,
    strict_parsing: bool,
}

#[pymethods]
impl SdMolSupplier {
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature = (path, sanitize=true, removeHs=true, strictParsing=true))]
    fn new(path: &str, sanitize: bool, removeHs: bool, strictParsing: bool) -> PyResult<Self> {
        let _ = (sanitize, removeHs); // accepted, no-op
        let file = std::fs::File::open(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        Ok(SdMolSupplier {
            inner: chematic_mol::SdfFileReader::new(std::io::BufReader::new(file)),
            strict_parsing: strictParsing,
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<pyo3::Py<pyo3::PyAny>>> {
        match self.inner.next() {
            None => Ok(None),
            Some(Ok(rec)) => {
                let mut props = rec.properties;
                if !rec.meta.name.is_empty() {
                    props.insert("_Name".to_string(), rec.meta.name);
                }
                let mol = Mol {
                    inner: std::sync::Arc::new(rec.mol),
                    props,
                };
                Ok(Some(pyo3::Py::new(py, mol)?.into_any()))
            }
            Some(Err(chematic_mol::MolParseError::Io(msg))) => {
                Err(pyo3::exceptions::PyIOError::new_err(msg))
            }
            Some(Err(_)) => {
                if self.strict_parsing {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "malformed SDF record",
                    ));
                }
                Ok(Some(py.None()))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SDWriter — RDKit-compatible SDF writer
// ---------------------------------------------------------------------------

/// Streaming SDF writer that writes ``Mol`` objects with their SD properties.
///
///     with chematic.SDWriter("output.sdf") as w:
///         mol = chematic.from_smiles("c1ccccc1")
///         mol.SetProp("Activity", "7.2")
///         w.write(mol)
///
/// 2D coordinates are computed automatically. The ``_Name`` property, if set,
/// is written into the MOL header line; other ``_``-prefixed properties are
/// omitted from the SD block.
#[pyclass(name = "SDWriter")]
pub struct SdWriter {
    writer: Option<std::io::BufWriter<std::fs::File>>,
}

#[pymethods]
impl SdWriter {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let file = std::fs::File::create(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        Ok(SdWriter {
            writer: Some(std::io::BufWriter::new(file)),
        })
    }

    fn write(&mut self, mol: &Mol) -> PyResult<()> {
        use std::io::Write as _;
        let w = self
            .writer
            .as_mut()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("SDWriter is already closed"))?;
        let layout = chematic_depict::compute_layout(&mol.inner);
        let coords: Vec<(f64, f64)> = layout.coords.iter().map(|p| (p.x, p.y)).collect();
        let name = mol.props.get("_Name").cloned().unwrap_or_default();
        let meta = chematic_mol::MolMetadata {
            name,
            comment: String::new(),
        };
        let record = chematic_mol::write_sdf_record(&mol.inner, &meta, &coords, &mol.props);
        w.write_all(record.as_bytes())
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    fn close(&mut self) {
        self.writer.take();
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[allow(unused_variables)]
    fn __exit__(
        &mut self,
        exc_type: Option<&pyo3::Bound<'_, pyo3::PyAny>>,
        exc_val: Option<&pyo3::Bound<'_, pyo3::PyAny>>,
        exc_tb: Option<&pyo3::Bound<'_, pyo3::PyAny>>,
    ) -> bool {
        self.close();
        false
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
    m.add_class::<SdMolSupplier>()?;
    m.add_class::<SdWriter>()?;
    m.add_function(wrap_pyfunction!(iter_sdf, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_str, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_batched, m)?)?;
    Ok(())
}
