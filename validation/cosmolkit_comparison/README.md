# chematic / COSMolKit / RDKit comparison

COSMolKit is an explicit comparison target. This directory holds two harnesses:

* **Corpus-scale comparison** (`run_corpus.py`, `compare_corpus.py`,
  `bench_corpus.py`, engine adapters in `corpus_engines.py`). Each engine
  dumps every operation for every row of a SMILES corpus; the dumps are
  scored against a reference engine (RDKit), and the same operations are
  timed in separate interpreters with rotating order.
* **Smoke harness** (below): a ten-row public corpus and a common JSONL
  contract for adapter validation.

The [ETKDGv3/MMFF94 conformer search comparison](3d-mmff-conformer-2026-10-24.md)
records source-built compatibility and quality profiles against COSMolKit
0.5.0, with row-level energy and convergence results, disjoint checks, and
the candidate-budget limitation.
The [bounded torsion-search follow-up](3d-mmff-torsion-search-2026-10-24.md)
adds a disjoint 12-molecule selection, indexed-graph and stereochemistry
checks, and independent RDKit MMFF94 re-scoring of both engines' coordinates.

The focused mixed V2000/V3000 SDF recognition fixture is
`corpora/mixed-v2000-v3000-v1.sdf`. With RDKit 2025.09.3 and COSMolKit 0.5.0
installed, run `python check_sdf_v3000_recognition.py` to verify the two
external readers against fixed title, graph, charge, and SD-field expectations.
The chematic source reader is covered by `cargo test -p chematic-mol --test
sdf_v3000_recognition`. This two-record development fixture establishes a
specific ingestion contract, not general recognition or stereo superiority.

The bond-level E/Z CIP comparison on the exposed 10,000-molecule corpus is
recorded in [ez-cip-2026-10-16.md](ez-cip-2026-10-16.md). It compares the
source-built default and accurate modes plus COSMolKit against RDKit's CIP
labels, keyed by bond endpoint atoms.

The [normalized-input V3000 comparison](v3000-radical-stereo-2026-10-17.md)
scores CheMatic, RDKit and COSMolKit on three fixed development corpora with
the same RDKit-normalized molecular input. It supersedes the earlier
[316-row raw-input report](v3000-three-engine-2026-10-17.md). Its reproduction
script is `scripts/compare_rdkit_cosmolkit_v3000.py`.
The [charged-bracket follow-up](v3000-charged-radicals-2026-10-17.md) fixes
four additional V3000 read-back identities on `stress-v1` with no regressions
on the three exposed corpora.
The [decorated-wildcard follow-up](v3000-decorated-wildcards-2026-10-17.md)
recovers seven more `stress-v1` identities while the browser and focused
stereo sets remain unchanged.
The [aromatic-ring wildcard follow-up](v3000-ring-wildcards-2026-10-17.md)
recovers five more `stress-v1` identities with no changed browser or focused
stereo rows.
The [pyramidal-sulfur stereo follow-up](v3000-sulfoxide-lone-pair-2026-10-17.md)
recovers four more `stress-v1` identities and checks 12 additional sulfoxide
variants against RDKit, COSMolKit, and CheMatic's own read-back.
The [charged-carbon follow-up](v3000-charged-carbon-2026-10-17.md) recovers
the final CheMatic-only `stress-v1` miss; CheMatic's correct rows now include
every RDKit and COSMolKit correct row on all three fixed V3000 corpora.
The [RDKit 2026.09.1 stereochemistry comparison](rdkit-2026-09-1-stereo-2026-10-17.md)
adds ChEMBL V3000 writer read-back, an independent NCI negative control,
corrected atom/bond CIP comparisons, a ChEMBL unmarked 3D reader sample, and
COSMolKit read-back of CheMatic blocks after omitting redundant V3000 `HCOUNT`.
The [neutral phosphorus follow-up](v3000-neutral-phosphorus-hcount-2026-10-18.md)
omits valence-inferred `[PH]` `HCOUNT`, with 15 additional exact COSMolKit
read-backs and no changed V3000 writer identities on two full corpora.
The [four-reader 3D MOL2 comparison](../openbabel_comparison/mol2-four-engine-stereo-2026-10-18.md)
uses identical MOL2 bytes to compare CheMatic, Open Babel, RDKit, and
COSMolKit on three stereo-eligible development selections.
The [two-vacancy non-tetrahedral comparison](two-vacancy-nontetra-2026-10-17.md)
checks SP2, TB3, and OH4 SMILES recognition across all tags, ligand orders, and
root/incoming positions, including randomized write-back.
The [native SMILES aromatic-stereo oracle follow-up](native-smiles-aromatic-stereo-oracle-2026-10-17.md)
checks a fused-arsenic ring with remote tetrahedral stereo under RDKit 2026.09.1
and guards against accepting an inverted center.
The [large O/S ring and explicit-H rotor follow-up](rdkit-2026-09-macrocycle-rotor-2026-10-18.md)
compares source-built CheMatic 1.1.7, RDKit 2026.09.1, and COSMolKit 0.5.0
on a 256-row operation audit and a 4,072-row three-operation check.

The focused [aromatic sulfur stereo read-back
check](aromatic-s-stereo-2026-10-16.md) compares whether each toolkit keeps
RDKit-enumerated `[s@+]` and `[s@@+]` SMILES distinct across five source rows.

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
