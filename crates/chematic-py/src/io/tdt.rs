use pyo3::prelude::*;
use std::sync::Arc;

use crate::Mol;

/// Streaming reader for Daylight TDT (Tagged Data) files.
#[pyclass(name = "TDTMolSupplier")]
pub(super) struct PyTdtMolSupplier {
    inner: chematic_mol::TdtRecordReader<std::io::BufReader<std::fs::File>>,
}

#[pymethods]
impl PyTdtMolSupplier {
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature = (path, nameRecord="NAME".to_string(), confId2D=-1, confId3D=-1, sanitize=true))]
    fn new(
        #[pyo3(from_py_with=crate::filesystem::path)] path: std::path::PathBuf,
        nameRecord: String,
        confId2D: i32,
        confId3D: i32,
        sanitize: bool,
    ) -> PyResult<Self> {
        let _ = sanitize;
        if confId2D >= 0 || confId3D >= 0 {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "chematic's Python Mol wrapper has no coordinate-carrying slot yet; \
                 confId2D/confId3D >= 0 (requesting coordinate data) is not supported here. \
                 Use the Rust API (chematic_mol::tdt::TdtRecordReader with read_2d/read_3d) \
                 for full 2D/3D coordinate support.",
            ));
        }
        let file = std::fs::File::open(&path).map_err(|e| {
            pyo3::exceptions::PyIOError::new_err(format!("{}: {e}", path.display()))
        })?;
        let options = chematic_mol::TdtReaderOptions {
            name_tag: (!nameRecord.is_empty()).then_some(nameRecord),
            read_2d: false,
            read_3d: false,
            strict_parsing: false,
            ..Default::default()
        };
        Ok(Self {
            inner: chematic_mol::TdtRecordReader::new(std::io::BufReader::new(file), options),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<pyo3::Py<pyo3::PyAny>>> {
        match self.inner.next() {
            None => Ok(None),
            Some(Ok(rec)) => {
                let mut props: std::collections::HashMap<String, String> =
                    rec.properties.into_iter().collect();
                if !rec.name.is_empty() {
                    props.insert("_Name".to_string(), rec.name);
                }
                let mol = Mol {
                    inner: Arc::new(rec.mol),
                    props: props.into(),
                };
                Ok(Some(pyo3::Py::new(py, mol)?.into_any()))
            }
            Some(Err(chematic_mol::TdtError::Io(msg))) => {
                Err(pyo3::exceptions::PyIOError::new_err(msg))
            }
            Some(Err(_)) => Ok(Some(py.None())),
        }
    }
}

/// Streaming writer for Daylight TDT files.
#[pyclass(name = "TDTWriter")]
pub(super) struct PyTdtWriter {
    writer: Option<chematic_mol::TdtRecordWriter<std::io::BufWriter<std::fs::File>>>,
    num_mols: usize,
}

#[pymethods]
impl PyTdtWriter {
    #[new]
    fn new(
        #[pyo3(from_py_with=crate::filesystem::path)] path: std::path::PathBuf,
    ) -> PyResult<Self> {
        let file = std::fs::File::create(&path).map_err(|e| {
            pyo3::exceptions::PyIOError::new_err(format!("{}: {e}", path.display()))
        })?;
        Ok(Self {
            writer: Some(chematic_mol::TdtRecordWriter::new(
                std::io::BufWriter::new(file),
                chematic_mol::TdtWriterOptions {
                    write_2d: false,
                    write_3d: false,
                    ..Default::default()
                },
            )),
            num_mols: 0,
        })
    }

    fn write(&mut self, mol: &Mol) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(closed_error)?;
        let name = mol.props.get("_Name").cloned().unwrap_or_default();
        let mut properties: Vec<(String, String)> = mol
            .props
            .iter()
            .filter(|(key, _)| key.as_str() != "_Name")
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect();
        properties.sort_unstable();
        writer
            .write_record(&chematic_mol::MoleculeRecord {
                mol: (*mol.inner).clone(),
                name,
                properties,
                coordinates_2d: None,
                coordinates_3d: None,
            })
            .map_err(io_error)?;
        self.num_mols += 1;
        Ok(())
    }

    #[pyo3(name = "SetProps")]
    fn set_props(&mut self, props: Option<Vec<String>>) -> PyResult<()> {
        self.writer
            .as_mut()
            .ok_or_else(closed_error)?
            .set_properties(props);
        Ok(())
    }

    #[pyo3(name = "SetWriteNames")]
    fn set_write_names(&mut self, value: bool) -> PyResult<()> {
        self.writer
            .as_mut()
            .ok_or_else(closed_error)?
            .set_name_tag(value.then(|| "NAME".to_string()));
        Ok(())
    }

    #[pyo3(name = "SetWrite2D")]
    fn set_write_2d(&mut self, _value: bool) {}

    #[pyo3(name = "SetNumDigits")]
    fn set_num_digits(&mut self, digits: usize) -> PyResult<()> {
        self.writer
            .as_mut()
            .ok_or_else(closed_error)?
            .set_precision(digits);
        Ok(())
    }

    #[pyo3(name = "NumMols")]
    fn num_mols(&self) -> usize {
        self.num_mols
    }

    fn flush(&mut self) -> PyResult<()> {
        if let Some(writer) = self.writer.as_mut() {
            writer.close().map_err(io_error)?;
        }
        Ok(())
    }

    fn close(&mut self) -> PyResult<()> {
        self.flush()?;
        self.writer.take();
        Ok(())
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<&pyo3::Bound<'_, pyo3::PyAny>>,
        _exc_val: Option<&pyo3::Bound<'_, pyo3::PyAny>>,
        _exc_tb: Option<&pyo3::Bound<'_, pyo3::PyAny>>,
    ) -> PyResult<bool> {
        self.close()?;
        Ok(false)
    }
}

fn closed_error() -> PyErr {
    pyo3::exceptions::PyValueError::new_err("TDTWriter is already closed")
}

fn io_error(error: std::io::Error) -> PyErr {
    pyo3::exceptions::PyIOError::new_err(error.to_string())
}
