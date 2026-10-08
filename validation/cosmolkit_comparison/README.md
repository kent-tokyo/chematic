# chematic / COSMolKit / RDKit comparison

COSMolKit is an explicit comparison target. This directory holds two harnesses:

* **Corpus-scale comparison** (`run_corpus.py`, `compare_corpus.py`,
  `bench_corpus.py`, engine adapters in `corpus_engines.py`). Each engine
  dumps every operation for every row of a SMILES corpus; the dumps are
  scored against a reference engine (RDKit), and the same operations are
  timed in separate interpreters with rotating order.
* **Smoke harness** (below): a ten-row public corpus and a common JSONL
  contract for adapter validation.

## Corpus-scale comparison

Install each engine in a Python that also has RDKit (round-trip operations
read engine output back with RDKit), then from this directory:

```bash
python run_corpus.py --engine rdkit     --corpus ../../scripts/chembl_accuracy_corpus_4999.smi --output rdkit.jsonl
python run_corpus.py --engine cosmolkit --corpus ../../scripts/chembl_accuracy_corpus_4999.smi --output cosmolkit.jsonl
python run_corpus.py --engine chematic  --corpus ../../scripts/chembl_accuracy_corpus_4999.smi --output chematic.jsonl
python compare_corpus.py --reference rdkit=rdkit.jsonl \
    --engine cosmolkit=cosmolkit.jsonl --engine chematic=chematic.jsonl --output summary.json
python bench_corpus.py --corpus ../../scripts/chembl_accuracy_corpus_4999.smi --blocks 5 --output bench.json
```

Rules: a cell is a `match` only when the value equals the reference exactly
(floats bit for bit; `match_1e-9` is reported separately). Typed refusals,
unsupported operations and errors are counted on their own and never as
matches. Round-trip operations (`*_rt`) are judged against the input
molecule, not against the reference engine's own writer, so RDKit's writer
is scored too. `API_NOTES` in `corpus_engines.py` names the API each
operation uses; where chematic has a native and an RDKit-compatible API the
RDKit-compatible one is used. `check_rdkit_smiles.py` and
`check_rdkit_avalon.py` are focused checks for the RDKit-compatible
canonical SMILES writer and Avalon fingerprint.

The October 8, 2026 results are recorded in
[benchmarks/2026-10-08-cosmolkit-parity.md](../../benchmarks/2026-10-08-cosmolkit-parity.md).

## Smoke harness

Run the two currently available engines from the repository root:

```bash
python3 validation/cosmolkit_comparison/run_engine.py --engine rdkit \
  --output /tmp/rdkit.jsonl
python3 validation/cosmolkit_comparison/run_engine.py --engine chematic \
  --output /tmp/chematic.jsonl
python3 validation/cosmolkit_comparison/score.py \
  /tmp/chematic.jsonl /tmp/rdkit.jsonl
```

Build an operation-level scorecard when more than two engines are available:

```bash
python3 validation/cosmolkit_comparison/scorecard.py \
  --result rdkit=/tmp/rdkit.jsonl \
  --result chematic=/tmp/chematic.jsonl \
  --reference rdkit --output /tmp/scorecard.json
```

The scorecard records each engine's version and source commit, then separates
per-operation status counts from reference comparisons. Unsupported or failed
operations are reported as `uncomparable`; they are never counted as matches.
The output is deterministic and is described by `scorecard.schema.json`.

The `rdkit_morgan_bits` operation uses chematic's promoted RDKit-exact Morgan
API when available; older installed chematic versions report it as
`unsupported` instead of silently substituting native ECFP4.

The corpus-scale harness above includes a COSMolKit adapter
(`corpus_engines.py`). Any other external adapter must emit the same schema
and identify unsupported operations as `unsupported`, never as a passing
value.

An external adapter can be plugged in without changing the harness:

```bash
python3 validation/cosmolkit_comparison/run_external.py \
  --engine cosmolkit --adapter 'python3 path/to/cosmolkit_adapter.py' \
  --output /tmp/cosmolkit.jsonl
python3 validation/cosmolkit_comparison/score.py \
  /tmp/chematic.jsonl /tmp/cosmolkit.jsonl
```

The adapter receives `--corpus PATH --engine NAME` and writes only common-schema
JSONL records to stdout. This keeps competitor installation and API decisions
outside the repository while keeping optional comparison runs reproducible.
The runner enforces a 120-second timeout by default (`--timeout-seconds` can
adjust it), verifies that every emitted record names the requested engine, and
publishes the output only after validation succeeds.

The result contract is versioned by `schema_version`. Each record contains the
engine version, source commit when available, corpus hash, input id/SMILES, and
an operation map. The corpus is a smoke test, not a claim of corpus-scale parity.

`corpus_manifest.json` is the authoritative contract for the smoke corpus. The
validator checks its filename, SHA-256, record order, IDs, and exact SMILES
values before accepting any engine output. Updating the corpus therefore
requires an explicit manifest update and makes accidental fixture drift visible
in CI. Run `python3 validate_manifest.py` for the fast, dependency-free gate.
