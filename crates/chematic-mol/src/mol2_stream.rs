//! Bounded, source-preserving iteration over concatenated MOL2 records.
use crate::{Mol2Error, Mol2ParseLimits, Mol2SourceRecord, parse_mol2_source_record_with_limits};

/// Limits for the whole in-memory input and each individual record.
#[derive(Debug, Clone, Copy)]
pub struct Mol2StreamLimits {
    pub max_input_bytes: usize,
    pub max_records: usize,
    pub record: Mol2ParseLimits,
}

impl Default for Mol2StreamLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_records: 100_000,
            record: Mol2ParseLimits::default(),
        }
    }
}

/// Record indices are zero-based; inner parser line numbers are record-local.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mol2StreamError {
    Parse {
        record_index: usize,
        error: Mol2Error,
    },
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl core::fmt::Display for Mol2StreamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Parse {
                record_index,
                error,
            } => write!(f, "MOL2 record {record_index}: {error}"),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "MOL2 stream {resource} limit exceeded: {actual} > {limit}"
            ),
        }
    }
}

impl std::error::Error for Mol2StreamError {}

/// Lazy record parsing over a borrowed, already-materialized UTF-8 input.
/// This is not a file-streaming reader. Each successful item owns its source.
/// Only whole, trimmed, case-insensitive MOLECULE marker lines delimit records.
/// Preamble belongs to the first record; inter-record whitespace/comments to
/// the preceding one. For a nonempty valid record stream, concatenating all
/// successful sources reproduces the input.
/// A malformed record yields an error and iteration may continue at the next
/// known boundary. Global resource errors terminate iteration. Empty/whitespace
/// input yields no records; nonempty input without a marker is a parse error.
pub struct Mol2SourceReader<'a> {
    remaining: &'a str,
    limits: Mol2StreamLimits,
    index: usize,
    done: bool,
}

impl<'a> Mol2SourceReader<'a> {
    pub fn new(input: &'a str) -> Self {
        Self::with_limits(input, Mol2StreamLimits::default())
    }

    pub fn with_limits(input: &'a str, limits: Mol2StreamLimits) -> Self {
        Self {
            remaining: input,
            limits,
            index: 0,
            done: false,
        }
    }

    fn fatal_limit(
        &mut self,
        resource: &'static str,
        actual: usize,
        limit: usize,
    ) -> Option<Result<Mol2SourceRecord, Mol2StreamError>> {
        self.done = true;
        Some(Err(Mol2StreamError::ResourceLimit {
            resource,
            actual,
            limit,
        }))
    }
}

impl Iterator for Mol2SourceReader<'_> {
    type Item = Result<Mol2SourceRecord, Mol2StreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if self.index == 0 && self.remaining.len() > self.limits.max_input_bytes {
            return self.fatal_limit(
                "input bytes",
                self.remaining.len(),
                self.limits.max_input_bytes,
            );
        }
        if self.remaining.trim().is_empty() {
            self.done = true;
            return None;
        }
        if self.index >= self.limits.max_records {
            return self.fatal_limit(
                "records",
                self.index.saturating_add(1),
                self.limits.max_records,
            );
        }
        let mut seen = false;
        let mut end = 0;
        for line in self.remaining.split_inclusive('\n') {
            if line.trim().eq_ignore_ascii_case("@<TRIPOS>MOLECULE") {
                if seen {
                    break;
                }
                seen = true;
            }
            end += line.len();
        }
        let (source, remaining) = self.remaining.split_at(end);
        self.remaining = remaining;
        self.done = remaining.is_empty();
        let record_index = self.index;
        self.index += 1;
        Some(
            parse_mol2_source_record_with_limits(source, &self.limits.record).map_err(|error| {
                Mol2StreamError::Parse {
                    record_index,
                    error,
                }
            }),
        )
    }
}

impl std::iter::FusedIterator for Mol2SourceReader<'_> {}
