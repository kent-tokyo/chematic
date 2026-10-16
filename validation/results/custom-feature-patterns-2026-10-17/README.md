# Custom feature-pattern validation

Pinned base: `d03106b59fdd4d06f40ab267bb540e20508bd988`.
RDKit #9601 and proposed PR #9615 describe ignored custom feature patterns and
the 32-bit pattern limit. The official FingerprintUtil source was consulted for
mask semantics. This candidate implements an owned native SMARTS feature set;
its circular fingerprint uses native FNV hashing, not RDKit Morgan identifiers.
Existing FCFP/Morgan functions and compatibility options remain unchanged.

`rdkit-baseline.json` and `comparison.json` retain an original RDKit 2025.09.3
wrapper probe. On ethanol, O and C definitions both yielded default atom masks
0/0/3, instead of requested 0/0/1 and 1/1/0. The installed Boost wrapper accepted
33 patterns while ignoring them. No C++ undefined-shift experiment was run.

The four initial feature-set tests failed against a retained older native
module with AttributeError because FeaturePatternSet was absent. Current base
source independently lacks the API; this is an API-addition baseline, not a
full rebuild of main or a claim that chematic silently ignored an existing
custom-pattern argument. Artifact details are in `before-change.json`.

Local results:

- `cargo test -p chematic-fp --lib --test feature_patterns --offline`: 346
  existing library + six new integration tests passed.
- Python 3.13.6: 211 passed, including 49 new feature-set cases and 162
  existing compatibility cases; zero failures/skips.
- 150 masks match an independent OR of RDKit SMARTS matches: 15 molecules,
  nine prefix pattern sets (including empty), and one 32-pattern set. Multi-atom
  matches, isotope/recursive queries, simple rings and aromatic/Kekulé spelling
  are included. This is a bounded oracle corpus, not universal SMARTS parity.
- Strict Clippy for chematic-fp/chematic-py, all targets: passed.
- Format, source stubs, refactor budgets, binding inventory/manifest, product
  naming, unsafe-surface and strict MkDocs checks: passed.

Tests also verify all 32 bits, refusal of the 33rd item, bounded consumption of
infinite iterables, ownership/reuse, original iteration exceptions, empty input,
pattern-order effects, budget refusal without partial results, target/pattern
limits and all 20 supported radius/width shapes. Changed/upper-bit patterns
affect native fingerprint bits. No RDKit fingerprint-bit parity is asserted.

`summary.json` and `source-sha256.json` pin source/artifact hashes. Native builds
used one job, no incremental compilation/debug info, opt-level 1, strip=none
and macOS dynamic-lookup flags. Compressed logs preserve actual checks. No
published-wheel, completed hosted-CI or main-merge claim is included.

Reproduce with the built native package on PYTHONPATH:

```sh
python3 validation/rdkit_issues/custom-feature-patterns/probe.py --compare
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q \
  crates/chematic-py/tests/test_feature_patterns.py \
  crates/chematic-py/tests/test_rdkit_compat.py
```
