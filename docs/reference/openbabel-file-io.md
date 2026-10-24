# Open Babel file-I/O comparison

chematic does not try to beat Open Babel by counting format names. Open Babel's
official documentation lists 146 formats, including many specialist and
one-directional adapters. The first chematic target is a smaller production
profile where correctness, failure accounting, safety, and speed can all be
verified.

## Evidence map

| Evidence | What it supports | What it does not support |
|---|---|---|
| [v1.0.42 same-process record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-openbabel-file-io-same-process-v1.0.42.md) | Historical bounded parse/write/round-trip latency on four checked fixtures | Current-release or broad-format performance |
| [SDF read follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-sdf-read-followup.md) | Two unfiltered 1,000-record source slices with semantic and timing gates | Published-package, second-host, or all-SDF superiority |
| [V3000 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-v3000-speed.md) | Bounded native-block and controlled-profile source measurements | A rank for diverse rows whose outputs differ |
| [stereo-preservation record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-15-openbabel-stereo-paired.md) | Paired RDKit-written V2000/V3000 stereo round trips | Independent proof of absolute stereo correctness |
| [latest V2000 experiment](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-18-openbabel-v2000-bond-tail.md) | Full retained evidence for a rejected source optimization and fresh Open Babel lanes | Adopted-main or release performance; the source adoption gate failed |
| [mmCIF numeric preservation](mmcif-numeric-preservation.md) | Authored finite binary64 coordinates, occupancy, B factors and cell parameters through Rust/Python/WASM record output | Complete mmCIF metadata/schema/source-byte retention, CLI support or speed |

The overview therefore reports Open Babel separately from RDKit and COSMolKit.
File-I/O results are not used to imply descriptor, fingerprint, reaction,
force-field, or conformer-generation superiority.

## Tier A scope

The initial contract covers SDF, MOL V2000, MOL V3000, MOL2, CML, CDXML, PDB,
and mmCIF. The machine-readable definition, fixture hashes, semantic fields,
and win criteria are in the
[machine-readable contract](https://github.com/kent-tokyo/chematic/blob/main/validation/openbabel_file_io_contract_v1.json).

Every input must finish in exactly one category: success, typed refusal,
unsupported, invalid input, or internal error. A successful round trip must
preserve the semantic fields declared for that format. Missing information may
be returned through a loss report or typed refusal; it must not disappear
silently.

## Gates

1. **Record accounting and provenance.** Run
   `python3 scripts/checks/check_openbabel_file_io_gate.py --output <path>`. This
   verifies identical fixture bytes, accepted-record counts, failure counts,
   tool versions, and execution boundaries. It is not a semantic or speed
   claim. The first current-source result is the
   [v1.0.42 record-accounting record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-openbabel-file-io-record-accounting-v1.0.42.json).
2. **Semantic round trip.** Compare graph, charge, isotope, bond order,
   stereochemistry, coordinates, and format-specific metadata after
   chematic→Open Babel and Open Babel→chematic conversion. Publish every row as
   exact, typed refusal, documented unsupported, or defect. The first bounded
   source-candidate gate observes both outputs through the same Rust reader for
   V3000, MOL2, CML, and CDXML. On the checked-in fixtures chematic preserves
   all observed fields; Open Babel 3.2.1 changes the MOL2 residue label. This is
   one small fixture per format, not broad corpus evidence. See the
   [clean-commit semantic record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-openbabel-file-io-semantics-v1.0.42.json).
3. **Malformed and bounded-input behavior.** Require zero panic/crash/internal
   errors and enforce time, memory, line, record, atom, and bond limits.
4. **Equivalent-work performance.** Measure cold start, parse, write,
   parse+write, and peak RSS separately. Use at least 21 alternating blocks;
   output must already pass the semantic gate. A win requires every block to
   be faster and the paired 95% speedup lower bound to exceed 1.0.

The first performance lane is a fresh-process CLI parse-plus-same-format-write
measurement for those four semantically checked fixtures. It includes process
startup and therefore does not establish parser-only or writer-only
superiority. All four lanes pass the 21-block rule in the
[clean-commit record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-openbabel-file-io-cli-roundtrip-v1.0.42.md).
Large-file throughput, peak RSS, the remaining Tier-A formats, and published
packages stay open.

The same-process comparison separates parse, write, and round-trip hot loops.
Build the Rust probe in release mode, then compile the small Open Babel 3.2.1
library harness with the installed `pkg-config` flags:

```sh
cargo build --release -p chematic-mol --example file_io_semantic_probe
c++ -std=c++17 -O3 $(pkg-config --cflags openbabel-3) \
  tools/openbabel_file_io_harness.cpp $(pkg-config --libs openbabel-3) \
  -o target/release/openbabel_file_io_harness
python3 scripts/benchmarks/bench_openbabel_file_io_same_process.py --blocks 21 --repeats 1000
```

Each process parses or writes the same fixture repeatedly and reports only the
internal timed loop. This removes executable startup and plugin discovery, but
it remains a tiny-fixture latency test rather than a large-file throughput or
memory measurement.

The first clean-commit run used 5,000 repetitions per process. chematic won
all 21 blocks in all 12 lanes, with median paired speedups from 2.59x to 6.38x
and a weakest paired 95% lower bound of 2.56x. See the
[bounded same-process record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-openbabel-file-io-same-process-v1.0.42.md).

The [1,000-row SDF throughput packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-sdf-speed.md)
extends the input beyond tiny fixtures. It compares native V2000 2D SDF batch
read, prepared rewrite and read-plus-rewrite loops, with 21 rotating blocks
and whole-process peak RSS. Both final rewrite outputs preserve every declared
field on the unfiltered corpus. Separate graph observations expose internal
parse-stereo differences, so parser timings are retained without a semantic
speed-win claim. This remains source-built, single-host evidence; other formats,
streaming memory, cold startup and published artifacts need their own gates.

The [SDF reader follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-sdf-read-followup.md)
corrects ambiguous overlapping-neighbor stereo and checks a second unfiltered
1,000-record slice. Candidate internal graphs and rewrites agree on both slices.
On the additional slice, native parse, write and roundtrip beat Open Babel 3.2.1
by 5.40x, 8.22x and 5.92x, respectively, in every one of 21 blocks. The original
Open Babel internal graph differences remain unranked. Incremental allocation
optimization has mixed results, and this packet does not establish general
Open Babel superiority.

A [paired stereo-preservation record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-15-openbabel-stereo-paired.md)
uses the same RDKit-written V2000/V3000 file bytes for both CLIs. On 155
stereo-sensitive corpus rows per profile, chematic preserves all RDKit
source/reopen identities; Open Babel 3.2.1 differs on 12 V2000 and 9 V3000
rows. Those differences represent four and three distinct source identities,
respectively. This is a one-host source comparison under an RDKit-defined
semantic contract, not independent proof of absolute stereo correctness.

The [V3000 native block packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-v3000-speed.md)
extends the comparison to individual MOL APIs, with collection framing outside
the timer. Candidate reads and rewrites agree on two unfiltered 1,000-row slices;
reading improves 1.055x and 1.087x over the corrected baseline in 21/21 blocks.
Open Babel's diverse-slice graph and rewrite failures leave all six comparator
lanes unranked. A separate controlled profile of 64 linear alkanes (16–1024 atoms)
passes all output gates and beats Open Babel 3.2.1 by 86.675x for parsing, 2.522x
for writing and 49.318x for roundtrip in 21/21 blocks. This controlled source
result does not establish general superiority; incremental writing results are
mixed, including a 0.9% regression on the additional diverse slice.

## Loss-aware MOL2 API

`parse_mol2_record` and `write_mol2_record` preserve Tripos atom types, partial
charges, residue identifiers and names, atom/bond status bits, formal charges
from `UNITY_ATOM_ATTR`, and opaque extension sections. The older `parse_mol2`
and `write_mol2` graph-oriented APIs remain available. Partial charges are no
longer rounded into formal charges; these are distinct MOL2 fields.

Until gates 2–4 pass for a named format and exact package artifact, the valid
claim is only that chematic has a bounded comparison program. It is not valid
to claim general Open Babel file-format superiority.

## External references

- [Open Babel supported formats](https://openbabel.org/docs/FileFormats/Overview.html)
- [Open Babel `obabel` CLI](https://openbabel.org/docs/Command-line_tools/babel.html)
- [Open Babel 3.2.1 release](https://github.com/openbabel/openbabel/releases/tag/openbabel-3-2-1)
