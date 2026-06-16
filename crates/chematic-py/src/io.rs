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
            mol: Mol { inner: Arc::new(rec.mol) },
            name: rec.meta.name,
            props: rec.properties,
        })
        .collect();
    SdfIter { records: records.into_iter() }
}

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------

/// Iterate over molecules in an SDF file by file path.
///
/// Reads the file, parses all valid records, and returns an iterator.
/// Invalid records are silently skipped.
///
///     for rec in chematic.iter_sdf("chembl.sdf"):
///         print(rec.smiles, rec.get("pChEMBL Value"))
///
///     # Build a Pandas DataFrame:
///     import pandas as pd
///     rows = [{"smiles": r.smiles, **r.properties()} for r in chematic.iter_sdf("data.sdf")]
///     df = pd.DataFrame(rows)
#[pyfunction]
pub fn iter_sdf(path: &str) -> PyResult<SdfIter> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("{path}: {e}")
        ))?;
    Ok(parse_to_iter(&content))
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
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySdfRecord>()?;
    m.add_class::<SdfIter>()?;
    m.add_function(wrap_pyfunction!(iter_sdf, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_str, m)?)?;
    Ok(())
}
