//! Indexed random access for seekable SDF files.

use std::io::{BufRead, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Mutex;

use crate::error::MolParseError;

use super::{SdfParseLimits, SdfRecord, SdfRecordReader};

/// Lightweight metadata for one indexed SDF record.
///
/// The index stores only byte ranges and the first (title) line. Molecular
/// graphs and SD data fields are parsed lazily when [`SdfDataset::get`] is
/// called.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecordMetadata {
    /// Zero-based record index after empty records have been skipped.
    pub index: usize,
    /// Byte offset of the first byte in the record.
    pub byte_offset: u64,
    /// Number of bytes in the record, excluding the `$$$$` delimiter line.
    pub byte_len: usize,
    /// First line of the MOL block, without its line ending.
    pub title: String,
}

/// Indexed, bounded-memory access to a seekable SDF file.
///
/// Opening the dataset scans delimiter lines once and keeps only one small
/// metadata entry per non-empty record. Calls to [`Self::get`] then seek to
/// and parse just the requested record. This complements [`super::SdfFileReader`]:
/// use the reader for one-pass streams and the dataset for cheap `len()`,
/// metadata inspection, random access, and bounded range reads.
pub struct SdfDataset {
    file: Mutex<std::fs::File>,
    records: Vec<SdfRecordMetadata>,
    limits: SdfParseLimits,
}

impl SdfDataset {
    /// Open and index an SDF file using the default finite limits.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, MolParseError> {
        Self::open_with_limits(path, SdfParseLimits::default())
    }

    /// Open and index an SDF file with explicit input, line, record-size, and
    /// record-count limits.
    pub fn open_with_limits(
        path: impl AsRef<Path>,
        limits: SdfParseLimits,
    ) -> Result<Self, MolParseError> {
        let file = std::fs::File::open(path.as_ref())
            .map_err(|error| MolParseError::Io(error.to_string()))?;
        let records = index_sdf_records(&file, limits)?;
        Ok(Self {
            file: Mutex::new(file),
            records,
            limits,
        })
    }

    /// Number of indexed non-empty records, including malformed records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the dataset contains no non-empty records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Return metadata without parsing the molecular graph.
    pub fn metadata(&self, index: usize) -> Option<&SdfRecordMetadata> {
        self.records.get(index)
    }

    /// Read and parse one record. `Ok(None)` means the index is out of range.
    pub fn get(&self, index: usize) -> Result<Option<SdfRecord>, MolParseError> {
        let Some(text) = self.record_text(index)? else {
            return Ok(None);
        };
        let mut records = SdfRecordReader::with_limits(&text, self.limits);
        records.next().transpose()
    }

    /// Read one raw record without the terminating `$$$$` line.
    pub fn record_text(&self, index: usize) -> Result<Option<String>, MolParseError> {
        let Some(metadata) = self.records.get(index) else {
            return Ok(None);
        };
        self.read_record_text(metadata).map(Some)
    }

    /// Parse a half-open record range while keeping memory bounded by the
    /// returned range. Out-of-range endpoints are clamped to `len()`.
    pub fn read_range(
        &self,
        start: usize,
        end: usize,
    ) -> Result<Vec<Result<SdfRecord, MolParseError>>, MolParseError> {
        let end = end.min(self.len());
        if start >= end {
            return Ok(Vec::new());
        }
        (start..end)
            .map(|index| {
                self.get(index).map(|record| {
                    record.expect("index originates from an in-range dataset interval")
                })
            })
            .map(|result| match result {
                Ok(record) => Ok(Ok(record)),
                Err(MolParseError::Io(error)) => Err(MolParseError::Io(error)),
                Err(error) => Ok(Err(error)),
            })
            .collect()
    }

    fn read_record_text(&self, metadata: &SdfRecordMetadata) -> Result<String, MolParseError> {
        let mut file = self
            .file
            .lock()
            .map_err(|_| MolParseError::Io("SDF dataset file lock is poisoned".to_string()))?;
        file.seek(SeekFrom::Start(metadata.byte_offset))
            .map_err(|error| MolParseError::Io(error.to_string()))?;
        let mut bytes = vec![0; metadata.byte_len];
        file.read_exact(&mut bytes)
            .map_err(|error| MolParseError::Io(error.to_string()))?;
        String::from_utf8(bytes).map_err(|error| MolParseError::Io(error.to_string()))
    }
}

fn index_sdf_records(
    file: &std::fs::File,
    limits: SdfParseLimits,
) -> Result<Vec<SdfRecordMetadata>, MolParseError> {
    let cloned = file
        .try_clone()
        .map_err(|error| MolParseError::Io(error.to_string()))?;
    let mut reader = std::io::BufReader::new(cloned);
    let mut records = Vec::new();
    let mut line = Vec::new();
    let mut absolute_offset = 0u64;
    let mut record_offset = 0u64;
    let mut record_len = 0usize;
    let mut record_has_content = false;
    let mut title: Option<String> = None;

    loop {
        line.clear();
        let bytes = reader
            .read_until(b'\n', &mut line)
            .map_err(|error| MolParseError::Io(error.to_string()))?;
        if bytes == 0 {
            if record_has_content {
                push_sdf_metadata(
                    &mut records,
                    record_offset,
                    record_len,
                    title.take().unwrap_or_default(),
                    limits,
                )?;
            }
            break;
        }

        let total = usize::try_from(absolute_offset)
            .unwrap_or(usize::MAX)
            .saturating_add(bytes);
        if total > limits.max_input_bytes {
            return Err(MolParseError::ResourceLimit {
                resource: "input bytes",
                actual: total,
                limit: limits.max_input_bytes,
            });
        }
        if bytes > limits.max_line_bytes {
            return Err(MolParseError::ResourceLimit {
                resource: "line bytes",
                actual: bytes,
                limit: limits.max_line_bytes,
            });
        }

        let content = line.strip_suffix(b"\n").unwrap_or(&line);
        let content = content.strip_suffix(b"\r").unwrap_or(content);
        if content == b"$$$$" {
            if record_has_content {
                push_sdf_metadata(
                    &mut records,
                    record_offset,
                    record_len,
                    title.take().unwrap_or_default(),
                    limits,
                )?;
            }
            absolute_offset = absolute_offset.saturating_add(bytes as u64);
            record_offset = absolute_offset;
            record_len = 0;
            record_has_content = false;
            title = None;
            continue;
        }

        if title.is_none() {
            title = Some(String::from_utf8_lossy(content).into_owned());
        }
        record_has_content |= !content.iter().all(u8::is_ascii_whitespace);
        record_len = record_len.saturating_add(bytes);
        if record_len > limits.max_record_bytes {
            return Err(MolParseError::ResourceLimit {
                resource: "record bytes",
                actual: record_len,
                limit: limits.max_record_bytes,
            });
        }
        absolute_offset = absolute_offset.saturating_add(bytes as u64);
    }

    Ok(records)
}

fn push_sdf_metadata(
    records: &mut Vec<SdfRecordMetadata>,
    byte_offset: u64,
    byte_len: usize,
    title: String,
    limits: SdfParseLimits,
) -> Result<(), MolParseError> {
    if records.len() >= limits.max_records {
        return Err(MolParseError::ResourceLimit {
            resource: "records",
            actual: records.len().saturating_add(1),
            limit: limits.max_records,
        });
    }
    records.push(SdfRecordMetadata {
        index: records.len(),
        byte_offset,
        byte_len,
        title,
    });
    Ok(())
}
