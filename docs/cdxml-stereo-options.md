# CDXML stereo inference options

[RDKit #9359](https://github.com/rdkit/rdkit/issues/9359) asks for an opt-in
way to infer stereo from ChemDraw's non-directional Bold and Hash bonds.
The issue is open as inspected on 2026-10-16. RDKit 2026.09.1 leaves these
unassigned by default; its discussion explicitly favors optional inference.

chematic already had `CdxmlParseOptions::infer_nondirectional_stereo` in
Rust. This change exposes that existing policy to Python and WASM. It also
fixes local reader defects found during independent validation: CDXML's
screen Y-down coordinates must be reflected for Cartesian parity perception;
solid/hashed sense is separate from Begin/End ownership; only the narrow
endpoint of a directional wedge owns its stereo marker. Returned coordinates
remain in their original CDXML convention.

```python
import chematic

mol = chematic.from_cdxml(text, infer_nondirectional_stereo=True)
```

The argument is keyword-only and defaults to `False`. Existing positional
and named `cdxml_str` calls retain their signatures. Only the first molecular
fragment is returned, matching the existing Python API.

```javascript
const mol = wasm.mol_from_cdxml_with_options(text, true);
const fragments = JSON.parse(wasm.cdxml_to_smiles_json_with_options(text, true));
mol.free();
```

Existing `mol_from_cdxml` and `cdxml_to_smiles_json` keep the default policy.
The new WASM functions retain the existing input/atom/fragment limits and
result shapes. No new Python export or public Rust error variant is added.

## Interpretation boundary

- Bold means toward the viewer; Hash and Dash mean away when the caller
  opts in. A unique endpoint with three or four neighbors is the candidate.
  If both endpoints qualify, the non-directional marker is ambiguous and
  does not assign either center. Missing or degenerate geometry also abstains.
- Non-directional B/E reversal does not change the inferred configuration.
  Only single-bond stereo glyphs participate; multiple-bond display styles
  keep their bond orders.
- Per the [ChemDraw bond-display specification](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/properties/Bond_Display.htm),
  `WedgeBegin`/`WedgedHashBegin` have the narrow end at B;
  `WedgeEnd`/`WedgedHashEnd` have it at E. Solid wedges point toward the
  viewer and hashed wedges point away. The reader normalizes End glyphs to
  that reference endpoint. This corrects earlier directional configurations,
  including mirror inversions and tags spuriously assigned at the wide end.
- This is heuristic molecular reading, not a general ChemDraw drawing
  interpreter. Bold may be visual emphasis, which is why inference is optional.
  Document/presentation APIs, binary CDX, reactions and CDXML writer stereo
  fidelity are outside this validation. Source spellings and B/E ordering are
  not preserved by the molecular graph.

## Reproduction and validation

The authored fixtures are in `validation/rdkit_issues/cdxml/`. The Rust suite
`cdxml_stereo_semantics` checks absolute reference CIP, screen-coordinate
preservation, End ownership and a multiple-bond negative control. Python and
real Node/WASM regressions cover default compatibility, opposite configurations,
B/E reversal, ambiguity, missing geometry, first-fragment behavior and errors.
WASM additionally checks batch output and unchanged safety limits.

`scripts/validate_cdxml_stereo_options.py` uses RDKit 2026.09.1 as an independent
reference for 216 drawings: three displays, two B/E orders, six rotations,
two mirror states and three node permutations. Each intended configuration
is expressed as an unambiguous directional CDXML control with the unique
carbon endpoint at the narrow end. The script compares both CIP and the full
isomeric graph after RDKit reparses chematic's SMILES. Default non-directional
controls are checked separately. Another 288 directional drawings check all
four Begin/End solid/hashed glyphs, including markers at a terminal endpoint.
Python and real Node/WASM match all 216 intended configurations, all 288
directional controls and all 216 default non-directional controls; the WASM
batch reader matches the same 216 intended configurations. The retained
pre-fix source runtime lacks the Python option and disagrees with 216 of the
288 directional controls. Local validation passes 70 Rust tests and 20 Python
tests (16 new and four existing CDXML contract tests), plus the Node/WASM
regression script. These are bounded source-candidate checks,
not a released-wheel or general ChemDraw compatibility result.

```sh
PYTHONPATH=/path/to/source/runtime:/path/to/pinned/rdkit \
  python3 scripts/validate_cdxml_stereo_options.py \
  --wasm-module /path/to/node/chematic_wasm.js \
  --baseline-pythonpath /path/to/earlier/source/runtime
```

The summary records source/runtime hashes and the retained pre-fix runtime's
missing Python option. The complete compressed packet retains every drawing,
its directional control and comparison output under
`validation/results/rdkit-cdxml-stereo-options-2026-10-16*`.
