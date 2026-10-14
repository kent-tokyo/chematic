# COSMolKit 0.5 feature-gap follow-up: fingerprint NumPy output

Date: 2026-10-14

Candidate: `codex/cosmolkit-gap-20261014`

Comparator: COSMolKit 0.5.0 public changelog, release notes, Python API, and
MIT-licensed source at tag `v0.5.0`

## Finding

The COSMolKit 0.5.0 review did not identify a missing fingerprint algorithm in
the compared release surface: chematic already provides Avalon, Pattern,
Layered, Morgan, atom-pair, topological-torsion, path, MACCS, and related batch
operations. The actionable Python gap was output consistency. COSMolKit's
fingerprint value exposes `to_numpy()` for every binary fingerprint, whereas
chematic only provided dense NumPy output for selected families.

This candidate adds one-dimensional `uint8` NumPy methods for the existing
ECFP6, FCFP4, atom-pair, topological-torsion, RDK, Pattern, Layered, Avalon,
and configurable RDKit-compatible Morgan paths. The arrays are copies and use
the same LSB-first order as the corresponding packed-byte method. No
fingerprint algorithm, hash, preprocessing rule, or packed output changed.

## Verification

- `cargo check -p chematic-py`: pass
- `cargo clippy -p chematic-py --all-targets -- -D warnings`: pass
- `cargo test -p chematic-py --no-default-features --lib`: 7 pass
- `pytest crates/chematic-py/tests/test_fp.py`: 80 pass
- full Python suite: 1,243 pass; one pre-existing optional-RDKit failure with
  local RDKit 2025.09.3 (`CalcNumHBA` for caffeine: 6 versus the fixture's
  documented Ertl value 3). The changed fingerprint test module is clean.

The tests expand each packed byte string with NumPy's little-bit-order rule and
require element-for-element equality for 20 fingerprint paths. Separate
cases cover 512-bit chiral Morgan output with and without count simulation and
reject non-byte-aligned Avalon widths instead of risking an indexing panic.

## Remaining reviewed gap

COSMolKit 0.5.0 also introduced arbitrary typed atom/bond properties and
opt-in reaction propagation. chematic intentionally exposes a smaller mutation
contract today. Adding that capability would change core public data models and
reaction ownership rules, so it remains a separate design task rather than an
incidental binding patch.

## Sources

- <https://github.com/cosmol-studio/COSMolKit/releases/tag/v0.5.0>
- <https://github.com/cosmol-studio/COSMolKit/blob/v0.5.0/CHANGELOG.md>
- <https://kit.cosmol.org/python>
