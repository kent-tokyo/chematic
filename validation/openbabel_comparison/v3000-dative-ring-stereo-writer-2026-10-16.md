# V3000 writer: dative-ring E/Z preservation

## Scope and result

This is a selected development set of 316 SMILES with demanding stereo syntax,
stored in `v3000-stereo-writer-stress-316.jsonl` (SHA-256
`f6728abce73d744260967520b5bdd07d25469c435fb5a875d7b98aaf2b82004c`).
It is not an independent or representative chemistry sample. Each writer
received the same SMILES; RDKit 2025.09.3 read each resulting V3000 MOL block.
The score compares mapless canonical isomeric SMILES to RDKit's direct read of
the source SMILES, so it requires the graph, isotope, charge, and stereochemistry
to agree. Open Babel 3.2.1 was tested both with its default settings and with
`--gen2D`.

| Writer | Exact identity | Explicit stereo-token sources |
| --- | ---: | ---: |
| CheMatic, before this change | 268/316 | 220/263 |
| CheMatic, after this change | 272/316 | 223/263 |
| Open Babel, default | 46/316 | 0/263 |
| Open Babel, `--gen2D` | 74/316 | 29/263 |

The four repaired source IDs are `stress-line-3212`, `stress-line-3954`,
`stress-line-3995`, and `stress-line-4069`. None of the previously correct
CheMatic cases regressed; there were no cases in this set where Open Babel
preserved exact identity and CheMatic did not. The remaining 44 CheMatic
identity misses are outside this E/Z repair. These counts do not establish
universal superiority over Open Babel.

## Defect and fix

The writer must mark a double bond with unknown E/Z as V3000 `CFG=2` when its
2D coordinates could otherwise be interpreted as a specific geometry. A
shortest-path check incorrectly treated a dative bond as closing a small ring,
so it omitted that mark. On one case, the dative neighbor was also counted as
an E/Z substituent, preventing the mark for a second reason. The ring check
and substituent check now exclude dative bonds. RDKit read the four affected
unmarked blocks with invented E/Z; it reads the corrected blocks with the
source's unspecified geometry.

## Reproduction

From the repository root, build `chematic-cli`, then run:

```sh
python3 scripts/compare_openbabel_v3000_stereo.py \
  --cli target/debug/chematic \
  --corpus validation/openbabel_comparison/v3000-stereo-writer-stress-316.jsonl \
  --report /tmp/chematic-v3000-dative-ring-fixed-316.json
```

The measured CheMatic binary SHA-256 was
`7ed99eb377518bb823a192295ee9184811bb210ad6f6627e0b3c45eeb106d6f9`.
The raw JSON report SHA-256 was
`dbb1a4401e431a780f34c0e90acdb30945c26f924453ca20bc325eb1b1899c54`.
The source base was `2251a9f5d2e706087efaf55cb0f8ab1ff20b8016`, plus
the changes in this branch. The 786 `chematic-mol` library tests passed,
including a regression test requiring `CFG=2` on the affected double bonds.
