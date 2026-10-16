use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::gc::{PyTraverseError, PyVisit};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::io::{BufReader, Read};
use std::sync::{Arc, Mutex};

type ReadError = Arc<Mutex<Option<PyErr>>>;

struct PythonRead {
    source: Arc<PythonSource>,
    error: ReadError,
}

struct PythonSource {
    owner: Py<PyAny>,
    read: Py<PyAny>,
}

impl Read for PythonRead {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }
        let result = Python::attach(|py| {
            let result = self.source.read.bind(py).call1((output.len(),))?;
            let bytes = result.cast::<PyBytes>().map_err(|_| {
                PyTypeError::new_err("stream.read(size) must return bytes (open in binary mode)")
            })?;
            let bytes = bytes.as_bytes();
            if bytes.len() > output.len() {
                return Err(PyValueError::new_err(
                    "stream.read(size) returned more than size bytes",
                ));
            }
            output[..bytes.len()].copy_from_slice(bytes);
            Ok(bytes.len())
        });
        result.map_err(|error| {
            // SdfFileReader reduces IO errors to strings; retain the original
            // Python exception separately so callers keep its type and object.
            if let Ok(mut slot) = self.error.lock() {
                *slot = Some(error);
            }
            std::io::Error::other("Python stream read failed")
        })
    }
}

/// Forward-only SDF iterator over a binary object with callable ``read(size)``.
///
/// Construction validates ``read`` without consuming input. Reads request at
/// most 8192 bytes; short reads are supported. The source remains alive while
/// the supplier exists and is never closed by the supplier. UTF-8 titles and
/// SD properties are preserved by the native SDF reader.
///
/// ``sanitize`` and ``removeHs`` are compatibility no-ops, as in SDMolSupplier.
/// Malformed records raise ValueError with strictParsing=True, or yield None
/// otherwise. IO, UTF-8 and resource errors always raise, even in lenient mode.
#[pyclass(name = "ForwardSDMolSupplier")]
pub(super) struct ForwardSdMolSupplier {
    inner: Option<chematic_mol::SdfFileReader<BufReader<PythonRead>>>,
    source: Option<Arc<PythonSource>>,
    read_error: ReadError,
    strict_parsing: bool,
}

#[pymethods]
impl ForwardSdMolSupplier {
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature = (fileobj, sanitize=true, removeHs=true, strictParsing=true))]
    fn new(
        fileobj: &Bound<'_, PyAny>,
        sanitize: bool,
        removeHs: bool,
        strictParsing: bool,
    ) -> PyResult<Self> {
        let _ = (sanitize, removeHs);
        let read = match fileobj.getattr("read") {
            Ok(read) if read.is_callable() => read.unbind(),
            Ok(_) => return Err(PyValueError::new_err("fileobj.read must be callable")),
            Err(error)
                if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(fileobj.py()) =>
            {
                return Err(PyValueError::new_err(
                    "fileobj must have a callable read method",
                ));
            }
            Err(error) => return Err(error),
        };
        let read_error = ReadError::default();
        let source = Arc::new(PythonSource {
            owner: fileobj.clone().unbind(),
            read,
        });
        let reader = PythonRead {
            source: source.clone(),
            error: read_error.clone(),
        };
        Ok(Self {
            inner: Some(chematic_mol::SdfFileReader::new(BufReader::with_capacity(
                8192, reader,
            ))),
            source: Some(source),
            read_error,
            strict_parsing: strictParsing,
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        let record = self.inner.as_mut().and_then(Iterator::next);
        if let Some(error) = self
            .read_error
            .lock()
            .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("stream error lock poisoned"))?
            .take()
        {
            return Err(error);
        }
        match record {
            None => Ok(None),
            Some(Ok(record)) => {
                let mut props = record.properties;
                if !record.meta.name.is_empty() {
                    props.insert("_Name".into(), record.meta.name);
                }
                let mol = crate::Mol {
                    inner: Arc::new(record.mol),
                    props: props.into(),
                };
                Ok(Some(Py::new(py, mol)?.into_any()))
            }
            Some(Err(error @ chematic_mol::MolParseError::Io(_)))
            | Some(Err(error @ chematic_mol::MolParseError::ResourceLimit { .. })) => {
                Err(super::sdf_dataset_error(error))
            }
            Some(Err(error)) if self.strict_parsing => Err(PyValueError::new_err(format!(
                "malformed SDF record: {error}"
            ))),
            Some(Err(_)) => Ok(Some(py.None())),
        }
    }

    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(source) = &self.source {
            visit.call(&source.owner)?;
            visit.call(&source.read)?;
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        self.inner = None;
        self.source = None;
    }
}
