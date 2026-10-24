# MCS query alternatives across relaxed atom and bond comparators

Source candidate based on `origin/main` `ccfd912e857c52f885b9c1b2e8ba3d34b69b545d`.
Comparator: RDKit 2026.09.1. This is a fixed, five-case development oracle,
not an independent holdout or a general FMCS accuracy score.
The [RDKit result](../validation/results/rdkit-mcs-query-alternatives-2026-10-25.json)
records the SMARTS, atom and bond counts, source re-match result, RDKit version,
and SHA-256 of the [generator](../scripts/validation/rdkit_mcs_query_alternatives.py).

## Defect and change

When relaxed comparators mapped unlike atoms or bond orders, the Rust MCS
search found a common topology, but its `QueryMolecule` described only the
first input's atom or bond. For `CO`, `NO`, and `SO` with `CompareAny`, the
previous query could not match `NO`. For `C=N`, `CN`, and `C#N` with
`CompareAny` bonds, it could not match all three inputs.

The query now uses the distinct **observed** atomic numbers and bond orders
at each mapped position. It covers every source without accepting an unseen
alternative solely because comparison was relaxed. The five RDKit oracle rows
include `[#6,#7,#16]-[#8]` and `[#6]=,-,#[#7]`; separate two-input rows confirm
that `PO` and `C#N`, respectively, remain excluded. Rust regression tests also
cover `AnyHeavyAtom` and `match_bonds=false`.

## Verification

```sh
/private/tmp/chematic-rdkit-2026091/bin/python \
  scripts/validation/rdkit_mcs_query_alternatives.py \
  --output /tmp/rdkit-mcs-query-alternatives.json
cargo test --locked -p chematic-smarts --lib mcs::comparator_tests
cargo test --locked -p chematic-smarts --lib
CARGO_INCREMENTAL=0 CARGO_PROFILE_TEST_DEBUG=0 cargo test --locked \
  -p chematic-wasm --lib concrete_mcs_projection_uses_first_input_alternative
CARGO_INCREMENTAL=0 CARGO_PROFILE_TEST_DEBUG=0 cargo test --locked \
  -p chematic-py --no-default-features --lib mcs_projection_tests
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo clippy --locked \
  -p chematic-smarts -p chematic-py -p chematic-wasm -- -D warnings
cargo fmt --all -- --check
python3 scripts/checks/check_benchmark_index.py
```

The focused comparator suite passed 15 tests. The full crate library suite
passed 245 tests. The pinned 660 RDKit-heavy and 1,980 native MCS comparator
cases retained their exact atom and bond counts, and each nonempty Rust query
re-matched every input. The cross-element re-match regression failed before
the fix; the new targeted cases pass after it. The Python and WASM Rust
projection tests also pass for reversed input orders and bond orders. Clippy
passes on all three changed Rust crates with warnings denied. The Python
`pytest` binding case is committed for hosted verification but was not
run locally because no source-built extension was installed.

## Boundary

This corrects Rust `QueryMolecule` output for the measured relaxed-comparator
cases. Python and WASM MCS currently return a concrete first-input projection,
so they do not expose these alternatives as generalized SMARTS. Exact SMARTS
string equality, the frozen FMCS holdout, charge/isotope query constraints,
search speed, and broad RDKit superiority were not established by this record.
