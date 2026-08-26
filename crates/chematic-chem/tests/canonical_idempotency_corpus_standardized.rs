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
/// **History**: an earlier version of this file set these to 60/68,
/// labeled "re-measured against main@743b77b (after #392 merged)" -- that
/// re-measurement never actually ran (a `git checkout` race silently
/// clobbered the source file mid-build). The 60/68 figures were in fact the
/// *pre-#392* baseline, carried forward unverified. Honestly re-measured on
/// `test/issue-393-canonical-idempotency-corpus` (`743b77b`, #392 included):
/// the true count was **615/519** -- a real, confirmed ~9x increase caused
/// by #392 exposing a pre-existing bug in `standardize.rs` (issue #399):
/// `neutralize_charges`, `normalize_zwitterion`, `normalize_groups`,
/// `remove_isotopes`, `reionize`, `uncharge`, `prefer_organic` and
/// `disconnect_metals` each rebuilt the molecule via a bare
/// `MoleculeBuilder` without carrying `stereo_neighbor_order`/
/// `bond_directions`/`stereo_groups` forward, so any stereocenter surviving
/// one of those stages lost its declared order and fell back to
/// `remove_hydrogens`'s adjacency-based reconstruction -- correct for
/// ring-closing stereocenters, transposed for ring-opening ones.
///
/// **#399 fix (this measurement)**: all 8 functions above now carry the
/// three stereo side tables forward (a simple bulk copy for the 7 that
/// preserve every atom/bond 1:1; `prefer_organic` now delegates to the
/// already-correct `extract_fragment`; `disconnect_metals` remaps
/// `bond_directions` bond-by-bond since it drops metal bonds). Re-measured:
/// **68/60** -- back down to the exact pre-#392 baseline, i.e. the ~547/459
/// failures #392 newly introduced are gone. The remaining 68/60 are
/// confirmed, by direct trace of representative cases (stereo, E/Z, and
/// plain-ring-fusion residuals), to be **unrelated to #399's root cause**:
/// `stereo_neighbor_order`/`bond_directions` survive `standardize()` fully
/// intact for these molecules, and disabling `remove_explicit_h` entirely
/// does not change the outcome -- so neither defect this issue fixed is in
/// play. ~76% of the residual lines share issue #395's exact syntactic
/// signature (an explicit bond symbol directly preceding a ring-closure
/// digit, e.g. `c1-2`); the rest are plain fused-ring systems with no
/// stereocenter at all. Left for #395 (or a new issue if #395 doesn't fully
/// explain them) rather than folded into this fix. Do not lower these again
/// without an honest full-corpus re-run, not an assumption; do not raise
/// them to hide a future regression either.
const DESCRIPTOR_CENSUS_KNOWN_FAILURES: usize = 60;
const CHEMBL_ACCURACY_KNOWN_FAILURES: usize = 68;

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
