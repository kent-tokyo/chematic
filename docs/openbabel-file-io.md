# Open Babel file-I/O comparison

CheMatic does not try to beat Open Babel by counting format names. Open Babel's
official documentation lists 146 formats, including many specialist and
one-directional adapters. The first CheMatic target is a smaller production
profile where correctness, failure accounting, safety, and speed can all be
verified.

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
   `python3 scripts/check_openbabel_file_io_gate.py --output <path>`. This
   verifies identical fixture bytes, accepted-record counts, failure counts,
   tool versions, and execution boundaries. It is not a semantic or speed
   claim. The first current-source result is the
   [v1.0.42 record-accounting record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-openbabel-file-io-record-accounting-v1.0.42.json).
2. **Semantic round trip.** Compare graph, charge, isotope, bond order,
   stereochemistry, coordinates, and format-specific metadata after
   CheMatic→Open Babel and Open Babel→CheMatic conversion. Publish every row as
   exact, typed refusal, documented unsupported, or defect. The first bounded
   source-candidate gate observes both outputs through the same Rust reader for
   V3000, MOL2, CML, and CDXML. On the checked-in fixtures CheMatic preserves
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

## Loss-aware MOL2 API

`parse_mol2_record` and `write_mol2_record` preserve Tripos atom types, partial
charges, residue identifiers and names, atom/bond status bits, formal charges
from `UNITY_ATOM_ATTR`, and opaque extension sections. The older `parse_mol2`
and `write_mol2` graph-oriented APIs remain available. Partial charges are no
longer rounded into formal charges; these are distinct MOL2 fields.

Until gates 2–4 pass for a named format and exact package artifact, the valid
claim is only that CheMatic has a bounded comparison program. It is not valid
to claim general Open Babel file-format superiority.

## External references

- [Open Babel supported formats](https://openbabel.org/docs/FileFormats/Overview.html)
- [Open Babel `obabel` CLI](https://openbabel.org/docs/Command-line_tools/babel.html)
- [Open Babel 3.2.1 release](https://github.com/openbabel/openbabel/releases/tag/openbabel-3-2-1)
