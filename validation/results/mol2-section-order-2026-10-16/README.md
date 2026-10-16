# MOL2 section-order evidence — 2026-10-16

This is bounded source evidence, not a release or full Amber validation.
Open Babel issue #2350 was OPEN when checked (see upstream-issue.json).

- `baseline-native.json`: actual default native output at
  `3e7f04b20593af4b07569d886b56afc1184468fe`, before production changes.
- `candidate-native.json`: actual option/default outputs from the new native
  example, using the same fixture JSON. Every default output equals baseline.
- `candidate-python.json`, `candidate-node.json`: actual newly built bindings.
- `cross-read-validation.json`: independent section-content comparison and
  Open Babel 3.2.1 V2000 readback of both profiles. Raw output and warnings are
  retained. Six charged graph controls lose formal charge only in the adjacent
  profile; the default retains the expected charges in all nine controls.
- `local-validation.json` and compressed logs: tested commands, observed
  results, build profile, fixture/source/artifact hashes and local gates.

The initial experiment changed the default and produced the raw failures in
`exploratory-baseline-tests.log.gz`. External charge loss ruled that policy out;
the final implementation is explicitly opt-in. Those failures demonstrate the
old ordering, not a baseline failure of a new API that did not yet exist.

## Reproduce

Run from a complete checkout of the candidate source. Build native Python and
WASM modules using the normal project tooling; set `PYTHONPATH` to the built
Python package and `CHEMATIC_WASM_MODULE` to the generated Node module. The local
run used debug optimized builds with stripping disabled, not published wheels.

```sh
cargo run -p chematic-mol --example mol2_section_order_dump --offline -- \
  validation/openbabel_issues/mol2-section-order/cases.json \
  > validation/results/mol2-section-order-2026-10-16/candidate-native.json
CHEMATIC_MOL2_ORDER_DUMP=validation/results/mol2-section-order-2026-10-16/candidate-node.json \
  node crates/chematic-wasm/tests/mol2_section_order.test.mjs
python3 validation/openbabel_issues/mol2-section-order/validate.py \
  --results validation/results/mol2-section-order-2026-10-16 --obabel obabel
```

For the baseline, use a separate checkout at the pinned commit. Copy
`validation/openbabel_issues/mol2-section-order/baseline_native_dump.rs` into
that checkout's `crates/chematic-mol/examples/` and run that example with the
absolute path to the candidate's `cases.json`. This harness uses only old
public APIs and produces the recorded baseline JSON schema. No third-party
fixtures or upstream code were copied.

Antechamber was not installed or run. Element-only generated atom types and
zero partial charges predate this change. Consumer warnings are retained;
atom counts alone are not evidence of chemical equivalence.
