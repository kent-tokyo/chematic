# chematic / RDKit / Open Babel speed comparison — 2026-09-05

This is a local remeasurement, not a new release claim. The chematic numbers
use the existing local `1.0.6` Python wheel in
`/private/tmp/chematic-speed-venv-2`; RDKit is measured in the same Python
process for each run. Open Babel is measured separately through the `obabel`
CLI, including process startup and conversion output.

Environment:

- macOS 26.5.2 arm64, Python 3.13.6
- chematic 1.0.6
- RDKit 2025.09.3
- Open Babel 3.2.1
- SMILES: `scripts/chembl_accuracy_corpus_4999.smi`, first 5,000 records
- SDF: RDKit `Contrib/PBF/testData/egfr.sdf`, 365 records

## In-process chematic / RDKit results

Each operation was run three times. Values below are medians; ranges show all
three measurements. The multiplier is competitor median divided by chematic
median, so values above 1 favour chematic. All runs parsed/canonicalized/read
the expected number of records.

| Operation | chematic median (range) | RDKit median (range) | chematic multiplier |
| --- | ---: | ---: | ---: |
| SMILES parse, 5,000 | 9.61 µs/mol (9.04–35.17) | 250.46 µs/mol (234.33–251.86) | 26.0x |
| canonical SMILES, 5,000 | 61.29 µs/mol (60.47–238.66) | 143.14 µs/mol (125.63–169.13) | 2.3x |
| SDF read, 365 | 35.47 µs/mol (24.99–66.58) | 401.19 µs/mol (317.84–482.38) | 11.3x |
| SDF write, 365 | 50.27 µs/mol (45.66–55.42) | 554.35 µs/mol (292.11–560.11) | 11.0x |

The broad ranges, especially the parse and canonical lanes, show host-load
sensitivity. These are medians of only three process-level samples and should
not be treated as an SLA or a universal cross-platform multiplier.

## Open Babel CLI boundary

Open Babel was measured with three fresh CLI processes per operation. The
reported values include process startup, CLI parsing, and output conversion;
they are not equivalent to an in-process parser or writer API measurement.

| Operation | Median total | Records | CLI time |
| --- | ---: | ---: | ---: |
| SMILES → SMILES | 1,998.40 ms | 5,000 | 399.68 µs/record |
| SMILES → canonical SMILES | 2,708.79 ms | 5,000 | 541.76 µs/record |
| SDF → SDF | 908.15 ms | 365 | 2,488.07 µs/record |

Reproduction:

```text
TMPDIR=/private/tmp /private/tmp/chematic-speed-venv-2/bin/python \
  scripts/bench_smiles_parse.py scripts/chembl_accuracy_corpus_4999.smi \
  --rdkit --n 5000 --json
TMPDIR=/private/tmp /private/tmp/chematic-speed-venv-2/bin/python \
  scripts/bench_canonical.py scripts/chembl_accuracy_corpus_4999.smi \
  --rdkit --n 5000 --json
TMPDIR=/private/tmp /private/tmp/chematic-speed-venv-2/bin/python \
  scripts/bench_sdf.py --rdkit --json
python3 scripts/bench_openbabel_cli.py --repeats 3
```

For a fair Open Babel library-level comparison, a binding or file-backed
benchmark with the same work definition is still required. The CLI results
are retained as a deployment-boundary baseline only.
