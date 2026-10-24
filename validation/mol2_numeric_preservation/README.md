# Native MOL2 numeric regeneration diagnostic

`cases.json` contains original MIT synthetic data. Eight cases participate in
same-format native/Open Babel rewrites with RDKit cross-reading; three finite
extremes are native-only. Newlines or aliases are not counted as extra cases.
Expected coordinate/partial-charge strings are parsed independently as binary64.
The indexed element and bond oracle is explicitly authored in the packet.

Build and run from the repository root (Python must have RDKit installed):

```sh
cargo build --locked -p chematic-mol --example mol2_numeric_dump
python scripts/validation/check_mol2_numeric_preservation.py \
  --native-executable target/debug/examples/mol2_numeric_dump \
  --obabel /path/to/obabel \
  --output /tmp/mol2-numeric-candidate.json.gz
cargo test --locked -p chematic-mol --test 'mol2_*'
python -m unittest scripts.tests.test_mol2_numeric_preservation
```

The runner bounds each CLI to 10 seconds and 1 MiB of combined output, saves
commands/stdout/stderr/resource outcomes, checks original source qualification,
and uses the same raw-field and RDKit observers for both engines. Its candidate
gate requires all native coordinate contracts and both editable record charge
contracts, plus indexed elements and bonds. Graph-only writers intentionally
emit zero placeholder charges, so their charge observations are not an
acceptance gate or a preservation claim. RDKit runs without sanitization,
H removal or substructure cleanup; no chemical-validation claim follows.

Saved reports:

- `../results/mol2-numeric-preservation-baseline-2026-10-24.json.gz`
- `../results/mol2-numeric-preservation-candidate-2026-10-24.json.gz`

Baseline native executable: debug build from `aeae3f53deb34579be60c054cdf22cd686d9254d`
with only the diagnostic example added, before production edits. Candidate:
working-tree build on the same parent. Both use the same packet and observers.
`--baseline --source-revision <baseline SHA> --probe-source <saved example.rs>`
records the baseline production hashes and suppresses the candidate-only gate.
The source revision flag is an operator-supplied provenance declaration; it
does not establish arbitrary executable provenance by itself. Here the baseline
was built and copied before edits, and the report records that executable hash.
Reports include executable, probe, packet, harness, production-source, Open Babel
binary and RDKit parser binary hashes and versions. Runtime paths are local
reproduction context, not committed executables. Gzip timestamps are fixed to zero.

Observed on macOS arm64, Rust debug build with one build job and no debug info or
incremental compilation: native editable output 8/8 coordinate and 8/8 partial
charge vectors exact; previous native and Open Babel 3.2.1 each 3/8 and 5/8.
The joint graph/coordinate/charge result is 8/8 versus 2/8 and 2/8.
This is not a timing result, holdout evaluation, format-wide guarantee or
published release result. Native-only extremes are not external wins.
