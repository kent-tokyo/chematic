//! Ordered, bounded selection from an indexed SDF dataset.
use super::{MolParseError, SdfDataset, SdfRecord};

#[derive(Debug, Clone, Copy)]
pub struct SdfSelectionLimits {
    pub max_records: usize,
    /// Sum of indexed source-record bytes, including duplicate selections.
    pub max_total_record_bytes: usize,
}
impl Default for SdfSelectionLimits {
    fn default() -> Self {
        Self {
            max_records: 10_000,
            max_total_record_bytes: 64 << 20,
        }
    }
}

#[derive(Debug)]
pub enum SdfSelectionError {
    Index {
        index: usize,
        len: usize,
    },
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    Read(MolParseError),
}
impl core::fmt::Display for SdfSelectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Index { index, len } => {
                write!(f, "SDF index {index} out of range for {len} records")
            }
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "SDF selection {resource} {actual} exceeds {limit}"),
            Self::Read(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for SdfSelectionError {}

/// Every requested position is retained; duplicate indices are parsed separately.
/// Parse failures stay in place. I/O failures abort the operation.
pub struct SdfIndexedBatch {
    pub indices: Vec<usize>,
    pub records: Vec<Result<SdfRecord, MolParseError>>,
}
impl SdfIndexedBatch {
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    pub fn valid_mask(&self) -> Vec<bool> {
        self.records.iter().map(Result::is_ok).collect()
    }
}
impl SdfDataset {
    /// Select any ordered list of source records with finite default limits.
    /// All indices and cumulative source bytes are checked before any seek/read.
    pub fn read_indices(&self, indices: &[usize]) -> Result<SdfIndexedBatch, SdfSelectionError> {
        self.read_indices_with_limits(indices, &SdfSelectionLimits::default())
    }
    pub fn read_indices_with_limits(
        &self,
        indices: &[usize],
        limits: &SdfSelectionLimits,
    ) -> Result<SdfIndexedBatch, SdfSelectionError> {
        self.check_selection(indices, limits)?;
        let records = indices
            .iter()
            .map(|&index| match self.get(index) {
                Ok(Some(record)) => Ok(Ok(record)),
                Ok(None) => Err(SdfSelectionError::Index {
                    index,
                    len: self.len(),
                }),
                Err(error @ MolParseError::Io(_)) => Err(SdfSelectionError::Read(error)),
                Err(error) => Ok(Err(error)),
            })
            .collect::<Result<_, _>>()?;
        Ok(SdfIndexedBatch {
            indices: indices.to_vec(),
            records,
        })
    }
    fn check_selection(
        &self,
        indices: &[usize],
        limits: &SdfSelectionLimits,
    ) -> Result<(), SdfSelectionError> {
        if indices.len() > limits.max_records {
            return Err(SdfSelectionError::ResourceLimit {
                resource: "records",
                actual: indices.len(),
                limit: limits.max_records,
            });
        }
        let mut total = 0usize;
        for &index in indices {
            let meta = self.metadata(index).ok_or(SdfSelectionError::Index {
                index,
                len: self.len(),
            })?;
            total = total
                .checked_add(meta.byte_len)
                .ok_or(SdfSelectionError::ResourceLimit {
                    resource: "record bytes",
                    actual: usize::MAX,
                    limit: limits.max_total_record_bytes,
                })?;
            if total > limits.max_total_record_bytes {
                return Err(SdfSelectionError::ResourceLimit {
                    resource: "record bytes",
                    actual: total,
                    limit: limits.max_total_record_bytes,
                });
            }
        }
        Ok(())
    }
}
