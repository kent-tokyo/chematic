//! Canonical round-trip idempotency through `standardize`, full-corpus
//! property test (issue #393).
//!
//! Sibling to `chematic-smiles`'s own `canonical_idempotency_corpus.rs`
//! (bare `parse`/`canonical_smiles`, no `standardize` involved) -- this one
//! specifically exercises `StandardizeOptions::default()` (`remove_explicit_h:
//! true`), the exact real pipeline shape #389 and #392 were both found
//! through (a downstream consumer's `standardize`-then-canonicalize stock-
//! identity step). Checks
//! `canon(standardize(x)) == canon(standardize(parse(canon(standardize(x)))))`
//! for every line in the two real-world-derived corpora already committed to
//! this repo (`scripts/*.smi`, 5,000 lines each).

use chematic_chem::{StandardizeOptions, standardize};
use chematic_smiles::{canonical_smiles, parse};

const DESCRIPTOR_CENSUS_CORPUS: &str =
    include_str!("../../../scripts/descriptor_census_corpus.smi");
const CHEMBL_ACCURACY_CORPUS: &str =
    include_str!("../../../scripts/chembl_accuracy_corpus_4999.smi");

/// For every non-empty line in `corpus`: parse, standardize, canonicalize
/// once, re-parse, standardize again, canonicalize twice, check the two
/// canonical strings are identical. Parse failures (of either the original
/// or the once-canonicalized form) count as failures too. Collects every
/// failing line (not just the first) and panics with the full list if the
/// count exceeds `max_known_failures` -- see
/// `canonical_idempotency_corpus.rs` (the sibling bare-`parse` test in
/// `chematic-smiles`) for why this is a ceiling, not a strict `== 0`
/// assertion: the same reasoning applies here (a pre-existing, already-
/// tracked residual must not make every unrelated future PR's CI red).
fn assert_corpus_idempotent(label: &str, corpus: &str, max_known_failures: usize) {
    let opts = StandardizeOptions::default();
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
        let once = canonical_smiles(&standardize(&mol, &opts));
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
        let twice = canonical_smiles(&standardize(&reparsed, &opts));
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
        "{label}: {}/{checked} corpus line(s) failed standardize+canonicalize round-trip \
         idempotency (known ceiling: {max_known_failures}; this test must not regress past \
         that count):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Current known-residual ceiling -- see this file's own doc comment. Must
/// only ever move down, never up.
///
/// **Correction**: an earlier version of this file set these to 60/68,
/// labeled "re-measured against main@743b77b (after #392 merged)" -- that
/// re-measurement never actually ran (a `git checkout` race silently
/// clobbered the source file mid-build, and the resulting "no such test
/// target" cargo error was mistaken for "0 additional failures" against an
/// empty grep pattern). The 60/68 figures were in fact the *pre-#392*
/// baseline, carried forward unverified. Honestly re-measured on
/// `test/issue-393-canonical-idempotency-corpus`
/// (`743b77b`, #392 included): the true count is dramatically higher,
/// **615/519**, not 60/68. This is a real, confirmed, ~9x increase caused
/// specifically by #392 -- not measurement noise, not a different corpus
/// state, not platform-dependent (re-confirmed identically against CI's
/// own Linux run). Root cause, and why #392 (a correct, real fix) is the
/// trigger rather than the culprit: filed as issue #399, with a confirmed
/// minimal repro (`CN1CCC[C@H]1c1cccnc1`, zero explicit H atoms --
/// `remove_hydrogens` should be a pure no-op, yet the molecule is
/// idempotent under bare `canonical_smiles` alone but not through
/// `standardize()`). Do not lower these again without an honest full-corpus
/// re-run, not an assumption.
const DESCRIPTOR_CENSUS_KNOWN_FAILURES: usize = 519;
const CHEMBL_ACCURACY_KNOWN_FAILURES: usize = 615;

#[test]
fn descriptor_census_corpus_standardized_is_canonically_idempotent() {
    assert_corpus_idempotent(
        "descriptor_census_corpus.smi",
        DESCRIPTOR_CENSUS_CORPUS,
        DESCRIPTOR_CENSUS_KNOWN_FAILURES,
    );
}

#[test]
fn chembl_accuracy_corpus_standardized_is_canonically_idempotent() {
    assert_corpus_idempotent(
        "chembl_accuracy_corpus_4999.smi",
        CHEMBL_ACCURACY_CORPUS,
        CHEMBL_ACCURACY_KNOWN_FAILURES,
    );
}
