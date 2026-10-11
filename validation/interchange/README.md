# Paired interchange development checks

This suite measures MOL2 source retention and MOL2 → CJSON against explicit
field expectations, using the source-built CheMatic 1.0.20 CLI and installed
Open Babel 3.2.1 CLI. These are exposed development fixtures. They are not an
independent holdout, a general chemistry oracle, or a superiority gate.

`protocol.json` fixes eight inputs, two routes, per-field tolerances and
resource budgets. Three inputs form the contiguous-ID common lane, two test
opaque extensions/nonsequential IDs separately, and three are negative inputs
under the declared local profile. Every engine/route/input gets a terminal
row, including refusals, parse/write failures, timeouts and invalid-profile
acceptance. A complete run has 32 rows; exit zero means complete measurement,
not that every field survived.

Expected values are authored explicitly in the protocol. Neither engine's
output supplies the oracle. The separate small output readers cover this
fixture subset only. Formal-charge expectations use the declared `N.4` +1 /
other listed atom types neutral projection; generic MOL2 charge interpretation
is outside scope. Negative/stress results can reflect dialect differences;
they do not by themselves establish an Open Babel bug.

Correspondence uses elements and unlabelled connectivity, then separately
checks bond orders, formal/partial charges, coordinates, IDs, names and types.
It never fits a mapping using the fields being measured. Symmetry and search
budget exhaustion remain unverified. Required coordinate/partial-charge
absolute tolerances are 0.00005 for MOL2 and 1e-12 for CJSON. Exact coordinate
f64 bits are a separate diagnostic; ordinary MOL2 decimal rounding within
tolerance does not fail the required coordinate check. Opaque marker retention
checks literal bytes only, not interpreted or editable extension semantics.

The runner bounds input, combined stdout/stderr while reading, wall time and
mapping search. It kills process groups on timeout/overflow. Reports publish
atomically without replacing an existing record. Protocol, runner, version,
format inventory, input bytes and both output pipes are archived beside the
report. Original fixture paths can be reconstructed from protocol entries and
the archived `<case>.input.mol2` files. Runtime artifacts are explicitly hashed;
the complete OS/library dependency closure is not archived.

Example rerun on the measured host (use a new output path):

```sh
.venv-v1gate/bin/python scripts/check_interchange_preservation.py \
  --chematic /private/tmp/chematic-interchange-integer-native/chematic \
  --chematic-evidence validation/results/interchange-final-build-source-2026-10-11.json \
  --obabel /opt/homebrew/bin/obabel \
  --expected-obabel-version 3.2.1 \
  --expected-obabel-sha256 733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5 \
  --obabel-runtime /opt/homebrew/Cellar/open-babel/3.2.1/lib/libopenbabel.8.dylib \
  --obabel-runtime /opt/homebrew/lib/openbabel/3.2.1/mol2format.so \
  --obabel-runtime /opt/homebrew/lib/openbabel/3.2.1/chemicaljsonformat.so \
  --output /private/tmp/interchange-rerun.json

.venv-v1gate/bin/python -m scripts.check_interchange_cross_read \
  --paired-report /private/tmp/interchange-rerun.json \
  --output /private/tmp/interchange-cross-rerun.json
```

The cross-read runner measures the three common inputs in three directions:
CheMatic CJSON → Open Babel MOL2, Open Babel CJSON → Open Babel MOL2 as a
control, and Open Babel CJSON → CheMatic immutable CJSON retention. It retains
all nine outcomes and verifies the selected artifacts before execution.
Source retention in the reverse direction does not interpret opaque metadata.

The first paired run found that CheMatic CJSON omitted partial charges, the
source charge-method label and molecule name. Checked conversion now emits a
complete declared dataset as `partialCharges.<method>` and the parsed name as
`properties.name`. Missing values are not filled and `NO_CHARGES` columns are
not relabeled as calculated charges. Reparsed partial charges require bitwise
f64 equality through verified atom correspondence. The method label is a
source declaration, not independently validated computation.

Cross reading also reproduced Open Babel 3.2.1 converting bond orders `1.0`
and `2.0` into MOL2 order 0. The CJSON writer now emits integral orders as JSON
integers; aromatic 1.5 stays fractional. The self-read control shows that Open
Babel's default CJSON-read/MOL2-write path also loses names and changes partial
charges from its own CJSON output. Those remaining losses are recorded rather
than claimed fixed by the integer change.

Records: [baseline](../results/interchange-paired-baseline-source-2026-10-11.json),
[partial-charge repair](../results/interchange-paired-candidate-source-2026-10-11.json),
[final candidate](../results/interchange-paired-integer-candidate-source-2026-10-11.json),
[cross-read before integers](../results/interchange-cross-read-before-integer-source-2026-10-11.json),
[cross-read after integers](../results/interchange-cross-read-after-integer-source-2026-10-11.json),
[final build verification](../results/interchange-final-build-source-2026-10-11.json).

Earlier integer-order milestone: 838 MOL/CLI tests, 217 source-wheel Python tests, five Node/WASM
suites, 59 shared report/error cases and 33 runner tests passed. Both paired
runs retain 32 terminal rows. The cross-read controls retain nine rows per
candidate: the integer-order change removes the two order-0 parse failures;
name/partial-charge losses on Open Babel readback remain. See the
[summary](../results/interchange-development-summary-2026-10-11.json).


The checked reverse adapter now validates root CJSON `partialCharges` datasets
and `properties.name`. CJSON → MOL2 preserves a single declared dataset,
explicit coordinates and representable names/method labels with reparse checks.
It does not select between multiple datasets, fill missing geometry or sanitize
names into MOL2 syntax. Those valid-but-unrepresentable records refuse even
under `allow_loss`; malformed known metadata is an input error. Other properties
remain opaque, as does `atoms.partialCharges`.

To reproduce the expanded cross-read lane, use the new
[build record](../results/cjson-mol2-build-source-2026-10-11.json) and
`/private/tmp/chematic-cjson-mol2-native/chematic` in the paired command above,
then add `--include-chematic-mol2` to the cross-read command. This adds both
CheMatic CJSON → CheMatic MOL2 and Open Babel CJSON → CheMatic MOL2 to each
common input, for 15 terminal rows. The original nine-row mode remains available.

[Current paired results](../results/interchange-paired-cjson-mol2-source-2026-10-11.json)
retain the same 32-row outcomes as the integer-order milestone. In the
[expanded cross-read results](../results/interchange-cross-read-cjson-mol2-source-2026-10-11.json),
all six CheMatic MOL2 projections preserve partial charges, source method label
and molecule name; all six Open Babel default readbacks lose those fields.
Mapped graph fields, required coordinate tolerance and coordinate-bit diagnostics
pass in the six new rows. The original nine outcomes and field checks remain
unchanged. All twelve MOL2 outputs still have some field loss; generated atom
names/types are an explicit CheMatic limitation. Required coordinate tolerances
remain distinct from the exact-f64 diagnostic.

Current local gates: 848 unique MOL/CLI tests, 230 isolated source-wheel Python
tests, five Node/WASM suites, 79 full-report/error parity cases (57 reports and
22 errors), and 33 paired-runner tests. Scoped Clippy passes with warnings
denied. The [reverse-adapter summary](../results/cjson-mol2-development-summary-2026-10-11.json)
indexes runtime hashes, logs and remaining boundaries. This is an uncommitted
source candidate tested on exposed fixtures, with no independent holdout or
published-release claim.
