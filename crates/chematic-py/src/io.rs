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
    pub stereo_diagnostics: Vec<chematic_perception::StereoDiagnostic>,
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

    /// Rejected wedge/hash stereocenters for this record.
    ///
    /// A list of ``{"atom_idx": int, "reason": str}`` dicts -- see
    /// :func:`from_mol_block_with_diagnostics` for the reason vocabulary.
    /// Empty unless a wedge/hash bond was present at some center and got
    /// rejected.
    fn stereo_diagnostics<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        crate::formats::stereo_diagnostics_to_py(py, &self.stereo_diagnostics)
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
            stereo_diagnostics: rec.stereo_diagnostics,
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
                        stereo_diagnostics: rec.stereo_diagnostics,
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
                        stereo_diagnostics: rec.stereo_diagnostics,
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
            // SDMolSupplier exposes RDKit-style molecules, not chematic's
            // optional stereo/3D diagnostics. Use the matching lightweight
            // reader so diagnostics do not inflate ordinary supplier reads.
            inner: chematic_mol::SdfFileReader::fast(std::io::BufReader::new(file)),
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
    props_filter: Option<Vec<String>>,
    force_v3000: bool,
    compute_2d: bool,
}

#[pymethods]
impl SdWriter {
    #[new]
    #[pyo3(signature = (path, compute2d=true))]
    fn new(path: &str, compute2d: bool) -> PyResult<Self> {
        let file = std::fs::File::create(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        Ok(SdWriter {
            writer: Some(std::io::BufWriter::new(file)),
            props_filter: None,
            force_v3000: false,
            compute_2d: compute2d,
        })
    }

    fn write(&mut self, mol: &Mol) -> PyResult<()> {
        use std::io::Write as _;
        let w = self
            .writer
            .as_mut()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("SDWriter is already closed"))?;
        let coords: Vec<(f64, f64)> = if self.compute_2d {
            let layout = chematic_depict::compute_layout(&mol.inner);
            layout.coords.iter().map(|p| (p.x, p.y)).collect()
        } else {
            Vec::new()
        };
        let name = mol.props.get("_Name").cloned().unwrap_or_default();
        let meta = chematic_mol::MolMetadata {
            name,
            comment: String::new(),
        };
        let filtered_props = self.props_filter.as_ref().map(|keys| {
            keys.iter()
                .filter_map(|k| mol.props.get(k).map(|v| (k.clone(), v.clone())))
                .collect()
        });
        let props = filtered_props.as_ref().unwrap_or(&mol.props);
        let record = if self.force_v3000 {
            chematic_mol::write_sdf_record_v3000(&mol.inner, &meta, &coords, props)
        } else {
            chematic_mol::write_sdf_record(&mol.inner, &meta, &coords, props)
        };
        w.write_all(record.as_bytes())
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Restrict which SD properties are written. Pass ``None`` to reset (write all).
    #[pyo3(name = "SetProps")]
    fn set_props(&mut self, props: Vec<String>) {
        self.props_filter = Some(props);
    }

    /// No-op: chematic always writes aromatic SMILES internally.
    #[pyo3(name = "SetKekulize")]
    fn set_kekulize(&mut self, _val: bool) {}

    /// Enable or disable automatic 2D coordinate generation.
    ///
    /// The default is `True` for compatibility. Set it to `False` when the
    /// caller only needs graph serialization and does not want depiction
    /// layout work on every record.
    #[pyo3(name = "SetCompute2D")]
    fn set_compute_2d(&mut self, val: bool) {
        self.compute_2d = val;
    }

    /// When ``True``, subsequent ``write()`` calls emit MOL V3000 (Extended
    /// Ctab) blocks instead of V2000 — required for molecules with more
    /// than 999 atoms or bonds.
    #[pyo3(name = "SetForceV3000")]
    fn set_force_v3000(&mut self, val: bool) {
        self.force_v3000 = val;
    }

    /// Flush buffered data to disk.
    fn flush(&mut self) -> PyResult<()> {
        use std::io::Write as _;
        if let Some(w) = self.writer.as_mut() {
            w.flush()
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        Ok(())
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
// SmilesMolSupplier / SmilesWriter — RDKit-compatible SMILES table I/O
// ---------------------------------------------------------------------------

/// Parse an RDKit-style delimiter string into a [`chematic_mol::Delimiter`].
///
/// A string consisting only of spaces/tabs (RDKit's own multi-char
/// delimiter *class* convention, e.g. `" \t"`) maps to
/// [`chematic_mol::Delimiter::Whitespace`] (runs collapsed) — chematic does
/// not offer RDKit's non-collapsing `keep_empty_tokens` behavior for a
/// multi-character class, only for the single-character
/// [`chematic_mol::Delimiter::Tab`]/[`chematic_mol::Delimiter::Custom`] cases.
/// Any other multi-character string is rejected explicitly rather than
/// silently truncated to its first character.
fn parse_delimiter(s: &str) -> PyResult<chematic_mol::Delimiter> {
    if !s.is_empty() && s.chars().all(|c| c == ' ' || c == '\t') {
        return Ok(chematic_mol::Delimiter::Whitespace);
    }
    let mut chars = s.chars();
    match (chars.next(), chars.next()) {
        (Some(','), None) => Ok(chematic_mol::Delimiter::Comma),
        (Some('\t'), None) => Ok(chematic_mol::Delimiter::Tab),
        (Some(c), None) if c.is_ascii() => Ok(chematic_mol::Delimiter::Custom(c as u8)),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unsupported delimiter {s:?}: chematic requires a single ASCII byte, or a \
             space/tab-only string (mapped to whitespace-class splitting)"
        ))),
    }
}

/// Streaming reader for SMILES table files (`.smi`/`.smiles`/`.csv`/`.tsv`/`.txt`).
///
/// RDKit-compatible API:
///
///     sup = chematic.SmilesMolSupplier("compounds.smi")
///     for mol in sup:
///         if mol is None:
///             continue  # malformed row
///         print(mol.smiles, mol.GetProp("_Name"))
///
/// ``sanitize`` is accepted for API compatibility but is a no-op (chematic's
/// SMILES parser already applies default aromaticity perception). Unlike
/// RDKit's own ``SmilesMolSupplier``, this reader is forward-only streaming
/// (no ``len()``/index access) — see ``chematic_mol::smiles_table``'s module
/// docs for the full list of deliberate, documented divergences from RDKit.
#[pyclass(name = "SmilesMolSupplier")]
pub struct PySmilesMolSupplier {
    inner: chematic_mol::SmilesRecordReader<std::io::BufReader<std::fs::File>>,
}

#[pymethods]
impl PySmilesMolSupplier {
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature = (path, delimiter=" ".to_string(), smilesColumn=0, nameColumn=1, titleLine=true, sanitize=true))]
    fn new(
        path: &str,
        delimiter: String,
        smilesColumn: usize,
        nameColumn: i64,
        titleLine: bool,
        sanitize: bool,
    ) -> PyResult<Self> {
        let _ = sanitize; // accepted, no-op (see doc comment)
        let delim = parse_delimiter(&delimiter)?;
        let file = std::fs::File::open(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        let options = chematic_mol::SmilesReaderOptions {
            delimiter: delim,
            smiles_column: smilesColumn,
            name_column: if nameColumn < 0 {
                None
            } else {
                Some(nameColumn as usize)
            },
            title_line: titleLine,
            // RDKit's own SmilesMolSupplier has no strict/lax toggle at all --
            // it unconditionally returns None for a bad row and continues.
            // This constructor matches that behavior; chematic's stricter
            // stop-on-error mode is only reachable via the Rust API.
            strict_parsing: false,
            ..Default::default()
        };
        Ok(PySmilesMolSupplier {
            inner: chematic_mol::SmilesRecordReader::new(std::io::BufReader::new(file), options),
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
                    props,
                };
                Ok(Some(pyo3::Py::new(py, mol)?.into_any()))
            }
            Some(Err(chematic_mol::SmilesTableError::Io(msg))) => {
                Err(pyo3::exceptions::PyIOError::new_err(msg))
            }
            Some(Err(_)) => Ok(Some(py.None())), // malformed row -> None, matches RDKit
        }
    }
}

/// Streaming writer for SMILES table files.
///
///     with chematic.SmilesWriter("output.smi") as w:
///         mol = chematic.from_smiles("c1ccccc1")
///         mol.SetProp("_Name", "benzene")
///         w.write(mol)
///
/// ``isomericSmiles=False`` and ``kekuleSmiles=True`` raise
/// ``NotImplementedError`` -- chematic has no non-isomeric or Kekule-form
/// SMILES writer at present, and this binding does not silently fall back
/// to the isomeric/aromatic form it can actually produce.
#[pyclass(name = "SmilesWriter")]
pub struct PySmilesWriter {
    writer: Option<chematic_mol::SmilesRecordWriter<std::io::BufWriter<std::fs::File>>>,
}

#[pymethods]
impl PySmilesWriter {
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature = (path, delimiter=" ".to_string(), nameHeader="Name".to_string(), includeHeader=true, isomericSmiles=true, kekuleSmiles=false))]
    fn new(
        path: &str,
        delimiter: String,
        nameHeader: String,
        includeHeader: bool,
        isomericSmiles: bool,
        kekuleSmiles: bool,
    ) -> PyResult<Self> {
        if !isomericSmiles {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "chematic has no non-isomeric SMILES writer mode; isomericSmiles=False is not supported",
            ));
        }
        if kekuleSmiles {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "chematic has no Kekule-form SMILES writer mode; kekuleSmiles=True is not supported",
            ));
        }
        let delim = parse_delimiter(&delimiter)?;
        let file = std::fs::File::create(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        let options = chematic_mol::SmilesWriterOptions {
            delimiter: delim,
            name_header: nameHeader,
            include_header: includeHeader,
            properties: Vec::new(),
        };
        Ok(PySmilesWriter {
            writer: Some(chematic_mol::SmilesRecordWriter::new(
                std::io::BufWriter::new(file),
                options,
            )),
        })
    }

    fn write(&mut self, mol: &Mol) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("SmilesWriter is already closed")
        })?;
        let name = mol.props.get("_Name").cloned().unwrap_or_default();
        let mut properties: Vec<(String, String)> = mol
            .props
            .iter()
            .filter(|(k, _)| k.as_str() != "_Name")
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        properties.sort_unstable();
        let record = chematic_mol::MoleculeRecord {
            mol: (*mol.inner).clone(),
            name,
            properties,
            coordinates_2d: None,
            coordinates_3d: None,
        };
        writer
            .write_record(&record)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Restrict (and order) which properties are written, mirroring RDKit's
    /// ``SetProps``. Defaults to writing no extra properties beyond
    /// SMILES + name, matching RDKit's own ``SmilesWriter`` default.
    #[pyo3(name = "SetProps")]
    fn set_props(&mut self, props: Vec<String>) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("SmilesWriter is already closed")
        })?;
        // SmilesRecordWriter doesn't expose a setter for its own options
        // struct post-construction; rebuild is unnecessary since `properties`
        // is the only field this method touches -- reach it directly.
        writer.set_properties(props);
        Ok(())
    }

    fn flush(&mut self) -> PyResult<()> {
        if let Some(w) = self.writer.as_mut() {
            w.flush()
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        Ok(())
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
// TDTMolSupplier / TDTWriter — RDKit-compatible Daylight TDT I/O
// ---------------------------------------------------------------------------

/// Streaming reader for Daylight TDT (Tagged Data) files.
///
/// RDKit-compatible API:
///
///     sup = chematic.TDTMolSupplier("compounds.tdt")
///     for mol in sup:
///         if mol is None:
///             continue  # malformed record
///         print(mol.smiles, mol.GetProp("_Name"))
///
/// ``nameRecord`` defaults to ``"NAME"`` here, **not** RDKit's own default
/// of ``""`` (meaning no tag populates the name) -- a source-confirmed,
/// deliberate divergence: RDKit's own `TDTWriter` always writes a `NAME`
/// tag by default, but its own `TDTMolSupplier` doesn't recognize it back
/// by default either, so a bare RDKit round trip silently loses the name.
/// See `chematic_mol::tdt`'s module doc comment for the full citation.
///
/// ``confId2D``/``confId3D`` (RDKit's opt-in-to-read-coordinates
/// convention: `< 0` disables, `>= 0` enables) are accepted, but a
/// non-negative value raises ``NotImplementedError`` -- chematic's Python
/// ``Mol`` wrapper has no coordinate-carrying slot yet, so 2D/3D
/// coordinates from a TDT file cannot be surfaced through this binding at
/// present (the Rust API, `chematic_mol::tdt::TdtRecordReader`, supports
/// them fully and is tested against a real coordinate-parsing bug found in
/// RDKit itself). ``sanitize`` is accepted for API compatibility but is a
/// no-op.
#[pyclass(name = "TDTMolSupplier")]
pub struct PyTdtMolSupplier {
    inner: chematic_mol::TdtRecordReader<std::io::BufReader<std::fs::File>>,
}

#[pymethods]
impl PyTdtMolSupplier {
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature = (path, nameRecord="NAME".to_string(), confId2D=-1, confId3D=-1, sanitize=true))]
    fn new(
        path: &str,
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
        let file = std::fs::File::open(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        let options = chematic_mol::TdtReaderOptions {
            name_tag: if nameRecord.is_empty() {
                None
            } else {
                Some(nameRecord)
            },
            read_2d: false,
            read_3d: false,
            strict_parsing: false,
            ..Default::default()
        };
        Ok(PyTdtMolSupplier {
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
                    props,
                };
                Ok(Some(pyo3::Py::new(py, mol)?.into_any()))
            }
            Some(Err(chematic_mol::TdtError::Io(msg))) => {
                Err(pyo3::exceptions::PyIOError::new_err(msg))
            }
            Some(Err(_)) => Ok(Some(py.None())), // malformed record -> None, matches RDKit
        }
    }
}

/// Streaming writer for Daylight TDT files.
///
///     with chematic.TDTWriter("output.tdt") as w:
///         mol = chematic.from_smiles("c1ccccc1")
///         mol.SetProp("_Name", "benzene")
///         w.write(mol)
///
/// Only name + properties round-trip through this binding at present --
/// chematic's Python ``Mol`` wrapper carries no 2D/3D coordinate data to
/// write (see ``TDTMolSupplier``'s doc comment); ``SetWrite2D`` is accepted
/// for API compatibility but has no visible effect until `Mol` gains a
/// coordinate-carrying mechanism.
#[pyclass(name = "TDTWriter")]
pub struct PyTdtWriter {
    writer: Option<chematic_mol::TdtRecordWriter<std::io::BufWriter<std::fs::File>>>,
    num_mols: usize,
}

#[pymethods]
impl PyTdtWriter {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let file = std::fs::File::create(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{path}: {e}")))?;
        Ok(PyTdtWriter {
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
        let writer = self.writer.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("TDTWriter is already closed")
        })?;
        let name = mol.props.get("_Name").cloned().unwrap_or_default();
        let mut properties: Vec<(String, String)> = mol
            .props
            .iter()
            .filter(|(k, _)| k.as_str() != "_Name")
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        properties.sort_unstable();
        let record = chematic_mol::MoleculeRecord {
            mol: (*mol.inner).clone(),
            name,
            properties,
            coordinates_2d: None,
            coordinates_3d: None,
        };
        writer
            .write_record(&record)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        self.num_mols += 1;
        Ok(())
    }

    /// Restrict (and order) which properties are written. Pass `None` to
    /// reset to RDKit's own `TDTWriter` default (write all).
    #[pyo3(name = "SetProps")]
    fn set_props(&mut self, props: Option<Vec<String>>) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("TDTWriter is already closed")
        })?;
        writer.set_properties(props);
        Ok(())
    }

    #[pyo3(name = "SetWriteNames")]
    fn set_write_names(&mut self, val: bool) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("TDTWriter is already closed")
        })?;
        writer.set_name_tag(if val { Some("NAME".to_string()) } else { None });
        Ok(())
    }

    /// No visible effect at present -- see the class doc comment.
    #[pyo3(name = "SetWrite2D")]
    fn set_write_2d(&mut self, _val: bool) {}

    #[pyo3(name = "SetNumDigits")]
    fn set_num_digits(&mut self, n: usize) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("TDTWriter is already closed")
        })?;
        writer.set_precision(n);
        Ok(())
    }

    #[pyo3(name = "NumMols")]
    fn num_mols(&self) -> usize {
        self.num_mols
    }

    fn flush(&mut self) -> PyResult<()> {
        if let Some(w) = self.writer.as_mut() {
            w.close()
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        Ok(())
    }

    fn close(&mut self) -> PyResult<()> {
        if let Some(w) = self.writer.as_mut() {
            w.close()
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        self.writer.take();
        Ok(())
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
    ) -> PyResult<bool> {
        self.close()?;
        Ok(false)
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
    m.add_class::<PySmilesMolSupplier>()?;
    m.add_class::<PySmilesWriter>()?;
    m.add_class::<PyTdtMolSupplier>()?;
    m.add_class::<PyTdtWriter>()?;
    m.add_function(wrap_pyfunction!(iter_sdf, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_str, m)?)?;
    m.add_function(wrap_pyfunction!(iter_sdf_batched, m)?)?;
    Ok(())
}
