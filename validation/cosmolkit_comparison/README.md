# Direct comparison smoke harness

This directory is the first checked-in slice of ROADMAP Phase 0. It defines a
small public smoke corpus and a common JSONL result contract for chematic, RDKit,
and a future COSMolKit adapter. It deliberately reports parse failures and
unsupported operations separately from value mismatches.

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

The COSMolKit runner is intentionally not bundled yet: its installation/access
method is not resolved in the roadmap. An adapter must emit the same schema and
must identify unsupported operations as `unsupported`, never as a passing value.

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
outside the repository while making the eventual COSMolKit run reproducible.

The result contract is versioned by `schema_version`. Each record contains the
engine version, source commit when available, corpus hash, input id/SMILES, and
an operation map. The corpus is a smoke test, not a claim of corpus-scale parity.

`corpus_manifest.json` is the authoritative contract for the smoke corpus. The
validator checks its filename, SHA-256, record order, IDs, and exact SMILES
values before accepting any engine output. Updating the corpus therefore
requires an explicit manifest update and makes accidental fixture drift visible
in CI. Run `python3 validate_manifest.py` for the fast, dependency-free gate.
