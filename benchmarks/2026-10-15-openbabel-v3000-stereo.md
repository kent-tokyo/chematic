# V3000 stereochemistry: CheMatic and Open Babel default writers

On one macOS arm64 host, both command-line tools received the same 300 SMILES
from the checked-in [Stereo Torture development suite](../validation/stereo_torture_suite_development.jsonl).
Each wrote V3000 MOL. Open Babel 3.2.1 was run both normally and with
`--gen2D`; CheMatic revision `2300923e` used its default `convert` path.
Python 3.13.6 and RDKit 2025.09.3 reread
each MOL block and compared canonical isomeric SMILES with the source. Atom map
numbers were cleared before comparison because map retention is a separate
metadata question. All 300 source SMILES parsed in RDKit; 226 source identities
contained explicit stereo tokens. Each case was run once; no timing claim is
made.

| Writer | All source identities retained | Explicit-stereo source identities retained |
| --- | ---: | ---: |
| CheMatic V3000 | 300 / 300 | 226 / 226 |
| Open Babel V3000 default | 74 / 300 | 0 / 226 |
| Open Babel V3000 with `--gen2D` | 92 / 300 | 19 / 226 |

For example, `F[C@H](Cl)Br` rereads from CheMatic's V3000 MOL with its
stereocenter intact. Open Babel's V3000 MOL for the same input has no bond
`CFG` marker, including with `--gen2D`, and RDKit rereads it as `FC(Cl)Br`.
The existing CheMatic V3000 writer produces this result; it is not an effect
of the 0D V2000 reader change in this branch.

Reproduce the [per-molecule JSON record](2026-10-15-openbabel-v3000-stereo.json),
including output MOL hashes and tool versions, with:

```bash
python3 scripts/compare_openbabel_v3000_stereo.py \
  --cli /path/to/chematic \
  --report /tmp/openbabel-v3000-stereo.json
```

This is a writer-to-RDKit-readback comparison on an exposed development
corpus. The writers use different default coordinate behavior. The added
`--gen2D` run checks Open Babel with generated 2D coordinates, but this record
does not establish an independent holdout result, reader accuracy for arbitrary
V3000 files, or overall superiority over Open Babel.
