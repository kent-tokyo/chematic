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
//! **Issue #395, fixed.** Running this test against both corpora for the
//! first time found a real, previously-undetected defect, independent of
//! #389/#392/#390: 73/5000 lines in `chembl_accuracy_corpus_4999.smi` and
//! 57/5000 in `descriptor_census_corpus.smi` were not idempotent. Every one
//! of the smallest failing examples shared an explicit-bond-order ring
//! closure (SMILES `-N`/`=N` immediately preceding a ring-closure digit,
//! e.g. `c1-2`).
//!
//! Root cause: `canonical.rs`'s `write_chain` decided whether a ring-closure
//! digit needed an explicit bond-order prefix by checking only the
//! currently-written atom's own aromaticity, never its ring-closure
//! partner's -- unlike the equivalent decision for a tree-edge child, which
//! correctly checks both endpoints (`parent_arom && child_arom`). A bare
//! ring-closure digit between two aromatic atoms is read back by the parser
//! as an *aromatic* bond, so a genuinely `Single`-order ring-closure bond
//! between two atoms that each individually happen to be aromatic (e.g. a
//! non-aromatic fusion bond joining two separately-aromatic ring systems,
//! `c1-2`) silently became an aromatic bond on re-parse whenever the writer
//! omitted its `-` marker. Fixed by checking both endpoints' aromaticity,
//! mirroring the tree-edge `implicit` computation exactly. Both corpora are
//! now **0/5000** failing -- a full fix, not a partial improvement.

use chematic_smiles::{canonical_smiles, parse};

const DESCRIPTOR_CENSUS_CORPUS: &str =
    include_str!("../../../scripts/descriptor_census_corpus.smi");
const CHEMBL_ACCURACY_CORPUS: &str =
    include_str!("../../../scripts/chembl_accuracy_corpus_4999.smi");

/// Known-residual ceiling -- see this file's own doc comment. Issue #395 is
/// now fully fixed (both corpora measured at 0 failures); kept as named
/// constants, not inlined `0`s, so a future regression here fails loudly
/// with the same message shape as when this was a nonzero ceiling. Must
/// only ever move down, never up.
const DESCRIPTOR_CENSUS_KNOWN_FAILURES: usize = 0;
const CHEMBL_ACCURACY_KNOWN_FAILURES: usize = 0;

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
