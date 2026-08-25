//! Canonical round-trip idempotency, full-corpus property test (issue #393).
//!
//! Motivation: #389 (isotope stripping) and #392 (tetrahedral stereo flip on
//! re-canonicalization) were both real, shipped correctness bugs that none of
//! chematic's existing tests caught, despite being simple to state and check:
//! `canon(parse(x)) == canon(parse(canon(parse(x))))`. Both were only found
//! via an external, large-scale (9.47M-compound) real-world corpus scan from
//! a downstream consumer, not from anything in this repo's own test suite.
//! This test runs that exact check against the two real-world-derived
//! corpora already committed to this repo (`scripts/*.smi`, 5,000 lines
//! each) -- no external download needed.
//!
//! Mirrors `chembl_roundtrip.rs`'s own established convention (parse ->
//! canonicalize -> re-parse -> re-canonicalize -> assert stability,
//! collect-all-failures-then-panic), scaled from a 50-molecule hand-picked
//! list to the full 5,000-line corpora.
//!
//! **Known residual, tracked as issue #395, not a regression from this
//! test's own addition.** Running this test against both corpora for the
//! first time found a real, previously-undetected defect, independent of
//! #389/#392/#390: 73/5000 lines in `chembl_accuracy_corpus_4999.smi` and
//! 57/5000 in `descriptor_census_corpus.smi` are not idempotent. Every one
//! of the smallest failing examples shares an explicit-bond-order ring
//! closure (SMILES `-N`/`=N` immediately preceding a ring-closure digit,
//! e.g. `c1-2`) -- an unconfirmed lead, not a diagnosis; see #395. Per this
//! project's own `_known_broken`-fixture convention (e.g. `dg.rs`'s ring-
//! fusion tests) and its "gate fail != regression until redesigned"
//! precedent (issue #70): rather than either (a) leaving this test
//! hard-failing on `main` -- which would make every unrelated future PR's
//! CI run red for a pre-existing, already-tracked defect it didn't cause --
//! or (b) `#[ignore]`-ing it and losing the "visible and tracked, not
//! silently ignored" property issue #393 explicitly asked for, this test
//! pins the CURRENT failure count as a ceiling: it stays green as long as
//! the count doesn't exceed what's measured today, and fails (loudly, with
//! every failing line) the moment it gets worse. Fixing #395 should lower
//! these constants, not raise them.

use chematic_smiles::{canonical_smiles, parse};

const DESCRIPTOR_CENSUS_CORPUS: &str =
    include_str!("../../../scripts/descriptor_census_corpus.smi");
const CHEMBL_ACCURACY_CORPUS: &str =
    include_str!("../../../scripts/chembl_accuracy_corpus_4999.smi");

/// Current known-residual ceiling (issue #395) -- see this file's own doc
/// comment. Must only ever move down, never up.
const DESCRIPTOR_CENSUS_KNOWN_FAILURES: usize = 57;
const CHEMBL_ACCURACY_KNOWN_FAILURES: usize = 73;

/// For every non-empty line in `corpus`: parse, canonicalize once, re-parse,
/// canonicalize twice, check the two canonical strings are identical. Parse
/// failures (of either the original or the once-canonicalized form) count
/// as failures too, not silently skipped -- a corpus line this test can't
/// even parse is itself worth knowing about. Collects every failing line
/// (not just the first) and panics with the full list if the count exceeds
/// `max_known_failures` -- see this file's own doc comment for why this is
/// a ceiling, not a strict `== 0` assertion.
fn assert_corpus_idempotent(label: &str, corpus: &str, max_known_failures: usize) {
    let mut failures: Vec<String> = Vec::new();
    let mut checked = 0usize;

    for (line_no, line) in corpus.lines().enumerate() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        checked += 1;

        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!("line {}: PARSE FAILED '{smi}': {e}", line_no + 1));
                continue;
            }
        };
        let once = canonical_smiles(&mol);
        let reparsed = match parse(&once) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!(
                    "line {}: RE-PARSE FAILED '{smi}' (once='{once}'): {e}",
                    line_no + 1
                ));
                continue;
            }
        };
        let twice = canonical_smiles(&reparsed);
        if once != twice {
            failures.push(format!(
                "line {}: NOT IDEMPOTENT '{smi}': once='{once}' twice='{twice}'",
                line_no + 1
            ));
        }
    }

    assert!(checked > 0, "{label}: corpus was empty, nothing checked");
    assert!(
        failures.len() <= max_known_failures,
        "{label}: {}/{checked} corpus line(s) failed canonical round-trip idempotency \
         (known ceiling: {max_known_failures} -- see issue #395; this test must not regress \
         past that count):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn descriptor_census_corpus_is_canonically_idempotent() {
    assert_corpus_idempotent(
        "descriptor_census_corpus.smi",
        DESCRIPTOR_CENSUS_CORPUS,
        DESCRIPTOR_CENSUS_KNOWN_FAILURES,
    );
}

#[test]
fn chembl_accuracy_corpus_is_canonically_idempotent() {
    assert_corpus_idempotent(
        "chembl_accuracy_corpus_4999.smi",
        CHEMBL_ACCURACY_CORPUS,
        CHEMBL_ACCURACY_KNOWN_FAILURES,
    );
}
