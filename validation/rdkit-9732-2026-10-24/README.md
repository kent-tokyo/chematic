# Checked reaction CXSMILES atom-group references

RDKit #9732 was OPEN when checked on 2026-10-24 JST. Its example is retained
as the first authored fixture: `O.C[C@H](O)CC>>CCC(C)=O |&1:2|`. RDKit
2026.09.1 writes `CC[C@H](C)O.O>>CCC(C)=O |&1:3|` in canonical mode. Index 3
then refers to the methyl branch, not the stereocentre at index 2.

The candidate adds explicit enhanced-atom-stereo reaction parsing/writing.
Indices are remapped through each molecule's actual emitted atom permutation
after sorting components within each role. The kind/number namespace is global
across all roles; group partitions can therefore span participants. This is an
explicit convention of the new API, not a blanket claim about reaction execution
or other toolkits' correlation models. Legacy APIs remain unchanged.

## Fixed packet

`cases.json` contains 20 authored positive inputs (two modes = 40 runs), including
plain/empty controls, plus 22 refusal cases. Both Python and WASM produce 40/40
identical output strings and pass the independent RDKit observer on every run.
The observer checks role, map, isotope, rooted atom identity, group partitions,
and CIP configuration. It ignores equivalent numeric group renaming and uniform
AND/OR group inversion. Six controls include harmless renaming/reordering,
permitted relative inversion, prohibited ABS inversion and wrong-atom changes.

| Captured writer | Observer result |
|---|---|
| Candidate Python | 40 matches |
| Candidate WASM | 40 matches |
| RDKit 2026.09.1 | 29 matches, 8 atom-membership differences, 3 cross-participant partition differences |
| Prior local Chematic runtime | New API unavailable; no chemistry-failure count |

The three partition differences are reported separately from the wrong-atom
problem because the candidate declares global group correlations. Raw captures
and before/after group observations are retained. The observer sanitizes and
assigns stereo after parsing; a group pointing to a non-stereogenic atom can
therefore be recorded as a lost group. The raw captured CX string retains the
incorrect reference for inspection. Numeric spelling is not a correctness
oracle, and these controls do not establish general stereo accuracy or speed.

The prior runtime was built from `29c5c192` (main `38b4c92c` plus the unrelated
V3000 either-bond change). It is not a separately rebuilt main/release binary.
Candidate sources, fixtures, observer scripts, captures and runtime module hashes
are pinned in `summary.json` and the gzip packets. No package publication or
main integration is implied.

## Reproduce

Build the candidate Python and Node/WASM artifacts, then:

```sh
python3 scripts/validation/validate_reaction_cxsmiles.py --engine chematic \
  --python-root /path/to/python-root --output /tmp/reaction-cx-python.json.gz
node scripts/validation/validate_reaction_cxsmiles_node.mjs \
  /path/to/node/chematic_wasm.js /tmp/reaction-cx-wasm.json.gz
/path/to/rdkit-2026.09.1/bin/python scripts/validation/validate_reaction_cxsmiles.py \
  --engine rdkit --output /tmp/reaction-cx-rdkit.json.gz
/path/to/rdkit-2026.09.1/bin/python scripts/validation/validate_reaction_cxsmiles.py \
  --observe /tmp/reaction-cx-python.json.gz /tmp/reaction-cx-wasm.json.gz \
  /tmp/reaction-cx-rdkit.json.gz --output /tmp/reaction-cx-observed.json.gz
cargo test --locked -p chematic-rxn --lib --test reaction_cxsmiles
PYTHONPATH=/path/to/python-root PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
  python3 -m pytest crates/chematic-py/tests/test_reaction_cxsmiles.py
CHEMATIC_WASM_MODULE=/path/to/node/chematic_wasm.js \
  node crates/chematic-wasm/tests/reaction_cxsmiles.test.mjs
```

Native tests separately reverse atom numbering and check role/map/group/CIP
identity. The profile refuses unsupported CX channels, invalid or overlapping
indices, repeated group declarations, non-tetrahedral members and malformed
numbers/titles. It does not retain arbitrary source text or full CX metadata.
