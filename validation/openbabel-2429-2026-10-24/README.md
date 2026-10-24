# V3000 either bonds: Open Babel #2429

The upstream issue was OPEN when checked on 2026-10-24 JST. The authored
`either-single.mol` represents butene with a wavy single bond adjacent to the
double bond. Sparse source IDs make the output-ID remapping observable.
Expected alkene stereo is unknown despite non-collinear coordinates.

The old Chematic reader dropped the single-bond `CFG=2` and inferred defined
E/Z. Its graph-only writer could also give a previously unknown alkene defined
E/Z after a coordinate-bearing round trip. `regression-before.json` records the
new Python regression failing against the retained prior runtime.

The fix retains single/double-bond `CFG=2` in normalized `MolMetadata`, suppresses
coordinate-based E/Z on double bonds adjacent to a wavy carrier's endpoints,
and records the resulting explicitly unknown alkene state on the graph. Rust
source-metadata writers retain the original wavy bond and endpoint order,
including 3D output. Python/WASM graph-only writers may express that same state
with double-bond `CFG=2`; they do not retain the original wave's location. An
ordinary graph without an explicit unknown marker may still express supplied
E/Z geometry. Existing exact-source APIs remain available and unchanged.

## Fixed comparison

`cases.json` has 16 wavy variants (left/right, endpoint direction, mirrored
geometry, 2D/3D), one double-either control and one unmarked geometry control.
Every input is four carbons and three bonds. The observer checks the expected
defined/unknown state before and after writing and reparsing. It retains raw
blocks and SMILES and does not compare literal spelling between engines.

| Runtime / pipeline | Read state | Read, write and reread state |
|---|---|---|
| Prior local Chematic runtime | 2/18 | 1/18 |
| Patched Chematic Python | 18/18 | 18/18 |
| RDKit 2026.09.1 | 18/18 | 18/18 |
| Open Babel 3.2.1 mol/can pipeline | 1/18 | 1/18 |

The prior binary was built from `c7e348f6` (main base `38b4c92c` plus the unrelated
3D zero-seed fix), not a separately rebuilt main/release executable. The
candidate starts independently from main `38b4c92c`. Engine/module, source,
fixture and observer hashes are retained in the summary and gzip packets.
Open Babel subprocesses each have a ten-second timeout. No timing or memory
measurement is claimed.

Native controls additionally cover a wavy carrier shared by two double bonds,
a remote wavy bond that must not suppress a different alkene, other bond
attributes, continuations, declared E/Z, source ID normalization, partial
coordinate defaults, 2D and 3D rewriting. The frozen parser differential keeps
unaffected rows; single-CFG=2 rows intentionally differ and have dedicated
semantic tests. The frozen writer profiles adopt only the new either-marker
selection policy; the string-building comparisons remain in place.

This is an alkene/V3000 contract, not SVG wave rendering, general tetrahedral
either perception, V2000 output,
arbitrary malformed CFG validation, or comprehensive stereo preservation.
The fixed packet does not establish toolkit-wide superiority or publication.

## Reproduce

Build the extension for the desired source revision, then run:

```sh
python3 scripts/validation/validate_v3000_either_single.py \
  --engine chematic --python-root /path/to/python-root --label candidate \
  --output /tmp/v3000-either-candidate.json.gz
/path/to/rdkit-2026.09.1/bin/python scripts/validation/validate_v3000_either_single.py \
  --engine rdkit --label rdkit --output /tmp/v3000-either-rdkit.json.gz
python3 scripts/validation/validate_v3000_either_single.py \
  --engine openbabel --label openbabel --output /tmp/v3000-either-ob.json.gz
cargo test --locked -p chematic-mol --lib --test v3000_either_single \
  --test v3000_lookup_contract --test v3000_semantic_interop
PYTHONPATH=/path/to/python-root PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
  python3 -m pytest crates/chematic-py/tests/test_v3000_either_single.py
CHEMATIC_WASM_MODULE=/path/to/node/chematic_wasm.js \
  node crates/chematic-wasm/tests/v3000_either_single.test.mjs
```
