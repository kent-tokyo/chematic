# Full-hydrogen MMFF conformer ranking versus RDKit 2026.09.1

This source-candidate experiment closes a major A6 quality gap without the
cost of embedding every hydrogen with distance geometry. CheMatic embeds the
heavy-atom graph in parallel, places implicit hydrogens with RDKit-compatible
geometry, samples terminal heteroatom-H rotamers deterministically, minimizes
the complete molecule with MMFF94, and ranks the resulting conformers by the
full-atom energy.

It is not a claim that CheMatic beats RDKit on every molecule. The aggregate
quality and exploratory throughput are favorable, but maximum residuals,
second-host replication, holdout evaluation, 21-block timing, memory, and
published-artifact checks remain open.

## Provenance

- Candidate base: `7facc71286d518d59152fd99a0f626a27988e54f`
- Candidate wheel SHA-256:
  `1db99d7ea85f1752d4e192da5015a8869e93e4ffd2ea46658295049e892911bc`
- RDKit: `2026.09.1`
- Host: Mac16,12, arm64, macOS 27.0.1, CPython 3.13.6, 10 logical CPUs
- Corpus: frozen A6 tier A+B, 265 molecules
- MMFF budget: 200 iterations during generation/ranking
- Quality scorer: RDKit MMFF94, with up to 20,000 iterations for both outputs

The machine-readable summary is
[`2026-10-12-a6-full-h-mmff-vs-rdkit.json`](2026-10-12-a6-full-h-mmff-vs-rdkit.json).

## Quality

Both engines use the same candidate count in each row.

| Candidates per molecule | Median CheMatic − RDKit best | p90 | Maximum | Within 1 kcal/mol | At/below RDKit best |
|---:|---:|---:|---:|---:|---:|
| 10 | −0.000012 kcal/mol | 1.041 | 9.896 | 237/265 | 216/265 |
| 20 | −0.000007 kcal/mol | 0.963 | 25.186 | 238/265 | 220/265 |

The median and p90 are favorable, and the 10-candidate lane improves the
previous source result from 193 to 237 molecules within 1 kcal/mol. The
maximum is not favorable: large flexible peptides, non-converged candidates,
and some rigid polycycles still select a higher-energy basin. “Complete
superiority” is therefore not established.

## Exploratory same-work speed

For 10 candidates per molecule, parsing excluded, both engines used all cores
within one molecule and the same 200-iteration MMFF budget.

| Engine | Full A6 median, 2 alternating runs | Successful candidates |
|---|---:|---:|
| CheMatic | 40.337 s | 2,649/2,650 |
| RDKit | 47.303 s | 2,629/2,650 |

The observed ratio is 1.17x in CheMatic's favor. Two repeats are deliberately
reported as exploratory; they do not satisfy the repository's 21-block,
paired-interval, two-host, published-artifact speed acceptance gate.

## Reproduce

Build and install the source wheel in an isolated environment containing
RDKit 2026.09.1. Generate and score quality with:

```bash
python scripts/public_package_3d_chematic.py \
  --output /tmp/a6-chematic.jsonl \
  --metadata-output /tmp/a6-chematic.metadata.json \
  --wheel /path/to/source-wheel.whl \
  --best-of-n 10 \
  --arms chematic_pipeline_v2_mmff94_full_h_lowest_of_n \
  --package-kind source_candidate \
  --source-revision REVISION \
  --source-diff-sha256 DIFF_SHA256

python scripts/a6_conformer_quality_vs_rdkit.py \
  --chematic-rows /tmp/a6-chematic.jsonl \
  --arm chematic_pipeline_v2_mmff94_full_h_lowest_of_n \
  --n-conformers 10 \
  --seed 20260801 \
  --jobs 8 \
  --output /tmp/a6-quality.json \
  --rows-output /tmp/a6-quality.jsonl
```

Run the paired throughput harness with:

```bash
python scripts/a6_conformer_speed_vs_rdkit.py \
  --output /tmp/a6-speed.json \
  --conformers 10 \
  --repeats 21
```

The next gate must target the unresolved high-energy rows, rerun on a second
host and an unused holdout, measure peak memory, and repeat from exact public
artifacts before promoting this to a release claim.
