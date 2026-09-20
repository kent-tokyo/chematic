//! Deterministic, record-oriented SMILES batch canonicalization.
//!
//! The batch API deliberately keeps parsing and canonicalization policy in
//! `chematic-smiles` so downstream stock/index consumers do not need to
//! duplicate it. It is lazy, preserves input order, and turns a malformed
//! record into a typed rejection without aborting later records.

use crate::canonical_smiles_stable_key;
use crate::{SmilesParseLimits, canonical_smiles, parse_with_limits};
use std::collections::BTreeMap;
use std::io::BufRead;

/// The outcome for one input record in a canonicalization batch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BatchCanonicalization {
    /// The record parsed and was canonicalized successfully.
    Accepted { canonical_smiles: String },
    /// The record was rejected without stopping the batch.
    Rejected { error: String },
}

/// One deterministic, input-order result from [`SmilesBatchCanonicalizer`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchCanonicalRecord {
    /// Zero-based position in the supplied iterator.
    pub input_index: usize,
    /// Original record text, retained for diagnostics and reproducibility.
    pub input: String,
    /// Accepted canonical output or a structured rejection message.
    pub result: BatchCanonicalization,
}

/// Complete, auditable result for an eager canonicalization batch.
///
/// Every supplied input has exactly one record. This operation has no skip or
/// policy-refusal path, so `skipped_count` and `refused_count` are always zero;
/// they remain explicit so callers can use the same accounting equation as
/// other bounded batch operations:
/// `input_count = accepted_count + rejected_count + refused_count + skipped_count`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchCanonicalResult {
    /// One input-order record for every supplied input.
    pub records: Vec<BatchCanonicalRecord>,
    /// Number of supplied inputs.
    pub input_count: usize,
    /// Records canonicalized successfully.
    pub accepted_count: usize,
    /// Records rejected by parsing or a configured resource limit.
    pub rejected_count: usize,
    /// Records refused by an operation policy (always zero for this operation).
    pub refused_count: usize,
    /// Records skipped without an attempted result (always zero for this operation).
    pub skipped_count: usize,
}

impl BatchCanonicalResult {
    fn from_records(records: Vec<BatchCanonicalRecord>) -> Self {
        let input_count = records.len();
        let accepted_count = records
            .iter()
            .filter(|record| matches!(record.result, BatchCanonicalization::Accepted { .. }))
            .count();
        Self {
            records,
            input_count,
            accepted_count,
            rejected_count: input_count - accepted_count,
            refused_count: 0,
            skipped_count: 0,
        }
    }

    /// True only if every input produced an accepted canonical result.
    pub const fn all_succeeded(&self) -> bool {
        self.accepted_count == self.input_count
    }
}

/// Reusable parsing/canonicalization policy for large SMILES batches.
///
/// The context is intentionally lightweight: parser state is created per
/// record, while resource limits are reused. This makes the API safe to reuse
/// across multiple stock files without carrying molecule state between them.
#[derive(Debug, Clone, Copy)]
pub struct SmilesBatchCanonicalizer {
    limits: SmilesParseLimits,
}

impl Default for SmilesBatchCanonicalizer {
    fn default() -> Self {
        Self::new(SmilesParseLimits::default())
    }
}

impl SmilesBatchCanonicalizer {
    /// Create a batch canonicalizer with explicit parser resource limits.
    pub const fn new(limits: SmilesParseLimits) -> Self {
        Self { limits }
    }

    /// Return the parser limits used by this context.
    pub const fn limits(&self) -> SmilesParseLimits {
        self.limits
    }

    /// Create a forward-only reader adapter over newline-delimited SMILES.
    ///
    /// The adapter is lazy and bounded by the caller's `BufRead` buffer
    /// policy. Parse failures are returned as ordinary rejected records;
    /// underlying I/O failures are returned as `Err` items because no reliable
    /// input record can be constructed for them.
    pub fn reader<R>(&self, reader: R) -> SmilesBatchReader<R>
    where
        R: BufRead,
    {
        SmilesBatchReader {
            reader,
            canonicalizer: *self,
            next_index: 0,
            line: String::new(),
            terminated: false,
        }
    }

    /// Lazily canonicalize records in input order.
    ///
    /// Each item is processed independently. A parse or resource-limit error
    /// becomes [`BatchCanonicalization::Rejected`], and later records are
    /// still emitted. Canonical output uses [`canonical_smiles`] directly:
    /// stereochemistry, isotope/charge annotations, explicit hydrogens, and
    /// disconnected components follow the existing parser/writer semantics;
    /// this API does not apply an additional standardization step.
    pub fn iter<I, S>(&self, records: I) -> impl Iterator<Item = BatchCanonicalRecord>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let limits = self.limits;
        records
            .into_iter()
            .enumerate()
            .map(move |(input_index, value)| {
                let input = value.as_ref().to_owned();
                let result = match parse_with_limits(&input, &limits) {
                    Ok(molecule) => BatchCanonicalization::Accepted {
                        canonical_smiles: canonical_smiles(&molecule),
                    },
                    Err(error) => BatchCanonicalization::Rejected {
                        error: error.to_string(),
                    },
                };
                BatchCanonicalRecord {
                    input_index,
                    input,
                    result,
                }
            })
    }

    /// Eagerly collect the same deterministic results as [`Self::iter`].
    pub fn canonicalize<I, S>(&self, records: I) -> Vec<BatchCanonicalRecord>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.iter(records).collect()
    }

    /// Eagerly canonicalize records with complete outcome accounting.
    ///
    /// This is the audited alternative to [`Self::canonicalize`] for callers
    /// that need to retain the input/result conservation invariant alongside
    /// the per-record results.
    pub fn canonicalize_with_result<I, S>(&self, records: I) -> BatchCanonicalResult
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        BatchCanonicalResult::from_records(self.canonicalize(records))
    }

    /// Build a deterministic exact-identity index from input records.
    ///
    /// Only [`canonical_smiles_stable_key`] results are indexed. A parsed
    /// molecule whose canonical representation is not proven stable is
    /// rejected rather than silently becoming an unsafe cache key. Duplicate
    /// identities retain every input position in ascending order.
    pub fn build_identity_index<I, S>(&self, records: I) -> IdentityIndexBuild
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut index = SmilesIdentityIndex::default();
        let mut rejected = Vec::new();
        let mut record_count = 0;

        for (input_index, value) in records.into_iter().enumerate() {
            record_count += 1;
            let input = value.as_ref().to_owned();
            match parse_with_limits(&input, &self.limits) {
                Ok(molecule) => match canonical_smiles_stable_key(&molecule) {
                    Some(key) => index.insert(key, input_index),
                    None => rejected.push(BatchCanonicalRecord {
                        input_index,
                        input,
                        result: BatchCanonicalization::Rejected {
                            error: "canonical identity is not stable for this molecule".to_string(),
                        },
                    }),
                },
                Err(error) => rejected.push(BatchCanonicalRecord {
                    input_index,
                    input,
                    result: BatchCanonicalization::Rejected {
                        error: error.to_string(),
                    },
                }),
            }
        }

        IdentityIndexBuild {
            index,
            record_count,
            rejected,
        }
    }
}

/// Deterministic exact-identity index keyed only by stable canonical SMILES.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SmilesIdentityIndex {
    entries: BTreeMap<String, Vec<usize>>,
}

impl SmilesIdentityIndex {
    fn insert(&mut self, key: String, input_index: usize) {
        self.entries.entry(key).or_default().push(input_index);
    }

    /// Return all input positions sharing an exact stable identity.
    pub fn positions(&self, key: &str) -> Option<&[usize]> {
        self.entries.get(key).map(Vec::as_slice)
    }

    /// Number of distinct stable identities in the index.
    pub fn unique_key_count(&self) -> usize {
        self.entries.len()
    }

    /// Number of indexed input records, including duplicates.
    pub fn record_count(&self) -> usize {
        self.entries.values().map(Vec::len).sum()
    }
}

/// Result of building a [`SmilesIdentityIndex`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityIndexBuild {
    /// Successfully indexed stable identities.
    pub index: SmilesIdentityIndex,
    /// Parse failures and fail-closed unstable identities, in input order.
    pub rejected: Vec<BatchCanonicalRecord>,
    /// Total input records observed, including rejected records.
    pub record_count: usize,
}

/// Lazy newline-delimited SMILES reader using the batch result contract.
pub struct SmilesBatchReader<R> {
    reader: R,
    canonicalizer: SmilesBatchCanonicalizer,
    next_index: usize,
    line: String,
    terminated: bool,
}

/// Terminal I/O failure from [`SmilesBatchReader`].
///
/// Records with indices below [`Self::next_unprocessed_index`] were emitted
/// before the failure. The failing line and every later input are unprocessed:
/// because a stream has no known total length, callers must not turn that
/// unknown suffix into successful or skipped records.
#[derive(Debug)]
pub struct SmilesBatchReaderError {
    /// First input index for which no terminal record was emitted.
    pub next_unprocessed_index: usize,
    /// Underlying I/O failure.
    pub source: std::io::Error,
}

impl std::fmt::Display for SmilesBatchReaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SMILES batch stream failed before input index {}: {}",
            self.next_unprocessed_index, self.source
        )
    }
}

impl std::error::Error for SmilesBatchReaderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

impl<R: BufRead> Iterator for SmilesBatchReader<R> {
    type Item = Result<BatchCanonicalRecord, SmilesBatchReaderError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.terminated {
            return None;
        }
        self.line.clear();
        match self.reader.read_line(&mut self.line) {
            Ok(0) => None,
            Ok(_) => {
                let input = self.line.trim_end_matches(['\n', '\r']).to_owned();
                let input_index = self.next_index;
                self.next_index += 1;
                let result = self
                    .canonicalizer
                    .iter(std::iter::once(input.as_str()))
                    .next()
                    .expect("one input must produce one batch record");
                Some(Ok(BatchCanonicalRecord {
                    input_index,
                    input,
                    result: result.result,
                }))
            }
            Err(source) => {
                self.terminated = true;
                Some(Err(SmilesBatchReaderError {
                    next_unprocessed_index: self.next_index,
                    source,
                }))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{BatchCanonicalization, SmilesBatchCanonicalizer};
    use std::io::{self, BufRead, Read};

    struct FailingReader;

    impl Read for FailingReader {
        fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("synthetic read failure"))
        }
    }

    impl BufRead for FailingReader {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            Err(io::Error::other("synthetic read failure"))
        }

        fn consume(&mut self, _amount: usize) {}
    }

    struct OneLineThenFail {
        data: &'static [u8],
        position: usize,
    }

    impl OneLineThenFail {
        fn new() -> Self {
            Self {
                data: b"CCO\n",
                position: 0,
            }
        }
    }

    impl Read for OneLineThenFail {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let available = self.fill_buf()?;
            let count = available.len().min(buf.len());
            buf[..count].copy_from_slice(&available[..count]);
            self.consume(count);
            Ok(count)
        }
    }

    impl BufRead for OneLineThenFail {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            if self.position == self.data.len() {
                Err(io::Error::other("synthetic read failure"))
            } else {
                Ok(&self.data[self.position..])
            }
        }

        fn consume(&mut self, amount: usize) {
            self.position = (self.position + amount).min(self.data.len());
        }
    }

    #[test]
    fn preserves_input_order_and_partial_errors() {
        let records = SmilesBatchCanonicalizer::default().canonicalize([
            "OCC",
            "C1CC",
            "C[C@H](O)F",
            "[13CH3][O-].[Na+]",
        ]);

        assert_eq!(records.len(), 4);
        assert_eq!(
            records
                .iter()
                .map(|record| record.input_index)
                .collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        assert_eq!(records[0].input, "OCC");
        assert!(matches!(
            records[0].result,
            BatchCanonicalization::Accepted { .. }
        ));
        assert!(matches!(
            records[1].result,
            BatchCanonicalization::Rejected { .. }
        ));
        assert!(matches!(
            records[2].result,
            BatchCanonicalization::Accepted { .. }
        ));
        assert!(matches!(
            records[3].result,
            BatchCanonicalization::Accepted { .. }
        ));
    }

    #[test]
    fn eager_result_has_complete_terminal_outcome_accounting() {
        let result =
            SmilesBatchCanonicalizer::default().canonicalize_with_result(["OCC", "C1CC", "CCN"]);
        assert_eq!(result.input_count, 3);
        assert_eq!(result.records.len(), result.input_count);
        assert_eq!(result.accepted_count, 2);
        assert_eq!(result.rejected_count, 1);
        assert_eq!(result.refused_count, 0);
        assert_eq!(result.skipped_count, 0);
        assert_eq!(
            result.input_count,
            result.accepted_count
                + result.rejected_count
                + result.refused_count
                + result.skipped_count
        );
        assert!(!result.all_succeeded());
    }

    #[test]
    fn applies_reusable_resource_limits_per_record() {
        let limits = crate::SmilesParseLimits {
            max_input_bytes: 3,
            ..Default::default()
        };
        let records = SmilesBatchCanonicalizer::new(limits).canonicalize(["CCO", "CCCC"]);

        assert!(matches!(
            records[0].result,
            BatchCanonicalization::Accepted { .. }
        ));
        assert!(matches!(
            records[1].result,
            BatchCanonicalization::Rejected { .. }
        ));
    }

    #[test]
    fn reader_is_lazy_and_keeps_later_records_after_parse_errors() {
        let input = std::io::Cursor::new("CCO\nC1CC\nCCN\r\n");
        let records: Vec<_> = SmilesBatchCanonicalizer::default()
            .reader(input)
            .map(Result::unwrap)
            .collect();

        assert_eq!(
            records
                .iter()
                .map(|record| record.input_index)
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(records[2].input, "CCN");
        assert!(matches!(
            records[1].result,
            BatchCanonicalization::Rejected { .. }
        ));
        assert!(matches!(
            records[2].result,
            BatchCanonicalization::Accepted { .. }
        ));
    }

    #[test]
    fn reader_io_error_is_terminal_and_exposes_unprocessed_boundary() {
        let mut reader = SmilesBatchCanonicalizer::default().reader(FailingReader);
        let error = reader
            .next()
            .expect("I/O failure must be surfaced")
            .expect_err("synthetic reader must fail");

        assert_eq!(error.next_unprocessed_index, 0);
        assert_eq!(error.source.kind(), io::ErrorKind::Other);
        assert!(
            reader.next().is_none(),
            "stream must fail closed after I/O error"
        );
    }

    #[test]
    fn reader_io_error_keeps_the_next_index_after_completed_records() {
        let mut reader = SmilesBatchCanonicalizer::default().reader(OneLineThenFail::new());
        let record = reader
            .next()
            .expect("first record must be emitted")
            .expect("first record must succeed");
        assert_eq!(record.input_index, 0);
        assert!(matches!(
            record.result,
            BatchCanonicalization::Accepted { .. }
        ));

        let error = reader
            .next()
            .expect("I/O failure must be surfaced")
            .expect_err("reader must fail after the first line");
        assert_eq!(error.next_unprocessed_index, 1);
        assert!(reader.next().is_none(), "stream must remain terminal");
    }

    #[test]
    fn identity_index_is_deterministic_and_fail_closed() {
        let build =
            SmilesBatchCanonicalizer::default().build_identity_index(["CCO", "OCC", "C1CC", "CCN"]);
        let ethanol_key = crate::canonical_smiles_stable_key(&crate::parse("CCO").unwrap())
            .expect("ethanol has a stable identity");

        assert_eq!(build.record_count, 4);
        assert_eq!(build.index.unique_key_count(), 2);
        assert_eq!(build.index.record_count(), 3);
        assert_eq!(build.index.positions(&ethanol_key), Some([0, 1].as_slice()));
        assert_eq!(build.rejected.len(), 1);
        assert_eq!(build.rejected[0].input_index, 2);
        assert!(matches!(
            build.rejected[0].result,
            BatchCanonicalization::Rejected { .. }
        ));
    }
}
