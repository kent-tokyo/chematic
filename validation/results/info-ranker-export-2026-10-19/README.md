# Typed InfoTheory ranking export — 2026-10-19

Implements the previously missing `InfoBitRanker.WriteTopBitsToFile` alongside
the typed InfoTheory work related to open
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601). No upstream algorithm
defect or whole-issue resolution is claimed.

## Evidence

299 Python tests passed: 23 new export tests and 276 existing ranker,
mask/chi-square/bias and compatibility tests. Default mypy
(`--follow-imports=skip`) passed for the changed module; a scoped strict check
also passed with existing root-package stub diagnostics ignored. An unscoped
strict attempt reports 81 pre-existing missing-generic-parameter errors in
`__init__.pyi`; package-wide strict typing is not claimed. Black checks on the
new test/oracle files, Python 3.9 grammar parsing of all three changed/new
Python files, `git diff --check` and the MkDocs site build passed.

`evidence.json` includes all 512 paired runs: 64 seeded datasets, four metrics,
and explicit/sparse RDKit source variants. Both runtimes' output files are
checked against the documented fixed-width/five-significant-digit table built
from their respective ranking arrays. After arrays are edited and extra votes
are accumulated, output must retain the original ranking snapshot.

- 1,024 format checks and 1,024 snapshot checks pass.
- 98 paired output files are byte-identical. Full file parity is not claimed:
  existing tie ordering and tiny score differences affect the other files.
- 508 complete reference permutations have equal per-class counts and scores
  agreeing within `1e-12`; maximum absolute score difference is
  `1.7763568394002505e-14`. Eight formatted score fields differ.
- Four reference results are not complete permutations of bit IDs (duplicates
  and omissions). Their rows are retained as `reference_bit_ids_not_permutation`;
  these are not silently excluded or counted as numeric matches. Export format
  and snapshot checks still run for these cases.
- Candidate results are asserted to contain each configured bit exactly once.

RDKit reference: 2026.09.1. Source inspected:
`Code/ML/InfoTheory/InfoBitRanker.cpp`, Git blob
`81037caa45399213eaf665a405d67268e88c5523` (GitHub master on 2026-10-19).
The writer uses field widths 12/12/10, score precision 5 and class header
suffixes outside the ten-character label field. The candidate keeps its
existing ascending-bit-ID tie policy.

An initial isolated reference probe attempted export on a newly constructed
ranker before `GetTopN` and exited with signal 11 (shell status 139). It was
not repeated and is not a claim about every RDKit build. The committed oracle
only exports valid populated rankings. The candidate explicitly refuses export
before a successful `GetTopN`, with no file creation/truncation; its tests cover
both no-vote and voted-but-unranked states.

## Runtime and reproduction

The changed ranker module, compatibility namespace and other Python modules
were copied from the candidate source to an isolated package. Its temporary
initializer exported the native extension's public symbols. The unchanged
published 1.1.8 CPython 3.13 macOS arm64 extension was reused solely for imports
(SHA-256 `ed5a04679987ee9b0a0a99daadca9f7a907b7af00e9e0c4180eaf49b9dbfcc55`).
No native APIs are invoked by the ranker or export implementation. The full
current package initializer/native build was not tested locally and remains a
CI integration check. This is not a wheel-release validation.

With a candidate package installed and RDKit available for the oracle:

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q \
  crates/chematic-py/tests/test_info_ranker_export.py \
  crates/chematic-py/tests/test_info_bit_ranker.py \
  crates/chematic-py/tests/test_info_ranker_masks.py \
  crates/chematic-py/tests/test_chi_square_bit_ranker.py \
  crates/chematic-py/tests/test_biased_bit_ranker.py \
  crates/chematic-py/tests/test_rdkit_compat.py
python3 scripts/validation/validate_info_ranker_export.py --output /tmp/ranker-export.json
```

Existing ranker limits remain: at most 65,536 features, 256 classes and
1,048,576 count entries. There is one copied numeric snapshot, replaced on
successful ranking. Export formats the complete bounded table before opening
the file; it does not promise atomic filesystem replacement.
