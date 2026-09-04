# RDKit / Open Babel speed comparison

Measured 2026-09-04 on Apple arm64, macOS 26.5.2, Python 3.13.6.

Versions:

- chematic: 1.0.3, built from the current workspace wheel
- RDKit: 2025.09.3
- Open Babel: 3.2.1

## Results

The chematic/RDKit rows use in-process Python bindings. Open Babel was measured
through its `obabel` CLI, including process startup and conversion output. The
Open Babel rows are therefore a separate CLI-conversion lane, not a proof of
parser-only or writer-only superiority.

| Operation | Records | chematic | RDKit | Open Babel CLI | Boundary |
|---|---:|---:|---:|---:|---|
| SMILES parse | 1,000 | 4.98 us/mol | 103.78 us/mol | — | in-process parse API |
| SMILES to SMILES conversion | 1,000 | — | — | 360 us/mol | CLI parse + serialization |
| Canonical SMILES | 1,000 | 23.70 us/mol | 26.37 us/mol | 480 us/mol | seven-process median for in-process APIs; Open Babel CLI `--canonical` |
| SDF graph read | 365 | 18.49 us/mol | 193.98 us/mol | — | in-process SDF reader |
| SDF to SDF conversion | 365 | — | — | 1,342 us/mol | CLI parse + serialization |
| SDF serialization-only write | 365 | 21.18 us/mol | 154.35 us/mol | — | in-process writer |

Canonical SMILES is the median of seven independent benchmark processes. The
chematic samples were 29.30, 24.89, 24.00, 23.64, 23.67, 23.64, and 23.70
us/mol; the RDKit samples were 26.57, 26.28, 26.37, 40.89, 26.04, 26.82, and
25.85 us/mol. The medians give chematic a 1.11x lead on this exact 1,000-row
slice and environment. This meets the scoped 1.1x target, but is not a
cross-platform claim. The other chematic/RDKit rows are single timed
invocations from the existing benchmark scripts.

Open Babel values are medians of three `/usr/bin/time -p` runs; the reported
real times were 0.36 s for 1,000 SMILES records, 0.48 s for 1,000
canonical-SMILES records, and 0.49 s for 365 SDF records.

Two allocation-oriented canonical-writer changes were also A/B tested and
rejected: caching sorted neighbor lists was 1.6% slower, and direct aromatic
symbol emission was about 1% slower. Neither change is retained in the source.

## v1.0.4 follow-up optimization

The fixed-version working tree was optimized further without changing the
public 1.0.3 version during measurement; the retained changes are released in
v1.0.4. Canonical automorphism feasibility now maintains an
inverse mapping instead of linearly scanning the forward mapping, and both
maps use compact `u32::MAX` sentinels. A focused 5,000-record engine run moved
from 131.70 ms to 120.26 ms after the inverse-map change (8.7% less total
time; geometric mean 10.44 to 9.34 us/molecule). Under a later, higher-load
sequence, sentinel compaction moved the repeated-run median from 235.76 ms to
231.94 ms (1.6%). These are source-level optimization checks, not replacements
for the v1.0.3-wheel headline above; search nodes, leaves, orbit tests, and
pruned-child counts were identical.

SDF SD-field serialization now reserves the aggregate property payload once
and formats directly into the record buffer. The same-process A/B benchmark
ran the previous and current implementations alternately over RDKit's
365-record `egfr.sdf`, 50 repeats per round. Five independent invocations all
improved, with speedups 1.482x, 1.260x, 1.248x, 1.145x, and 1.151x (median
1.248x). Each record was first checked for byte-for-byte identity.

## Reproduction

```bash
# Verify versions
python3 -c 'import chematic, rdkit; print(chematic.__version__, rdkit.__version__)'
obabel -V

# chematic/RDKit in-process lanes
python3 scripts/bench_smiles_parse.py scripts/descriptor_census_corpus.smi --rdkit --n 1000 --json
python3 scripts/bench_canonical.py scripts/descriptor_census_corpus.smi --rdkit --n 1000 --json
TMPDIR=/private/tmp python3 scripts/bench_sdf.py --rdkit --json

# Same-process SDF writer A/B (old implementation is embedded as the baseline)
cargo run --release -p chematic-mol --example sdf_writer_ab -- /path/to/egfr.sdf

# Open Babel CLI lanes
head -1000 scripts/descriptor_census_corpus.smi | obabel -ismi -osmi -O /dev/null
head -1000 scripts/descriptor_census_corpus.smi | obabel -ismi -osmi --canonical -O /dev/null
obabel -isdf /path/to/egfr.sdf -osdf -O /dev/null
```

For a fair API-level comparison, Open Babel needs a binding-level harness or a
file-backed benchmark with the same work definition. These CLI numbers are
useful deployment/command-line baselines, but should not be used to claim that
one library's parser or writer is intrinsically faster than another's.
