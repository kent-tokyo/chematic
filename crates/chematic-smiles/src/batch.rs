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
}

impl<R: BufRead> Iterator for SmilesBatchReader<R> {
    type Item = std::io::Result<BatchCanonicalRecord>;

    fn next(&mut self) -> Option<Self::Item> {
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
            Err(error) => Some(Err(error)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{BatchCanonicalization, SmilesBatchCanonicalizer};

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
