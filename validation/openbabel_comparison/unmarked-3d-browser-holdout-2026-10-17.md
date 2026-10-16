# Open Babel 3.2.1: browser-corpus 3D SDF stereo holdout

## Measured result

This development comparison selected 100 tetrahedral and 100 alkene-stereo
source identities from the independent 10,000-row browser corpus. RDKit
2025.09.3 embedded each source in 3D, removed CTAB stereo markers, and read
the resulting V2000 and V3000 blocks back. Only records whose RDKit read-back
matched the source's canonical isomeric SMILES were retained. CheMatic,
Open Babel 3.2.1, and COSMolKit 0.5.0 then read the **same SDF bytes**.

| Input | CheMatic before | CheMatic after | Open Babel | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| V2000, exact mapless stereo | 196/200 | **197/200** | 197/200 | 200/200 |
| V3000, exact mapless stereo | 196/200 | **197/200** | 197/200 | 200/200 |

Every CheMatic and Open Babel output retained the mapless molecular graph
(200/200). The only Open Babel-only identity before this change was source
row 3944, with three coupled stereocenters in a six-membered ring. CheMatic's
first 3D perception pass assigned the two outer centers. A narrowly gated
second pass now lets those configurations distinguish the remaining ring
center. The new V2000 and V3000 fixtures pin its exact identity.

There are no changes to the other 199 identities in this selection. The three
remaining misses (source rows 2508, 442, and 2718) are shared by CheMatic and
Open Babel: both assign an additional tetrahedral configuration that the source
left unspecified. COSMolKit leaves those centers unspecified. This is a
source-identity comparison, so the extra tags count as mismatches; the 3D
coordinates alone do not prove what the source intended for an unspecified
center.

The earlier, separate 300-record SDF set remains 300/300 for CheMatic on
both V2000 and V3000 after this change. Open Babel remains 290/300 and
288/300 respectively. The new holdout establishes **parity** with Open Babel
on these 200 identities, and the prior set retains a measured CheMatic lead.

## Change and regression guard

CTAB still uses its original first pass. Its second pass considers only an
unbracketed, degree-three ring carbon with one implicit hydrogen and two
unbracketed ring neighbors that both received stereo on the first pass. MOL2
retains its existing two-pass behavior. The 14,072 pinned graph inputs and all
five-vertex graphs retain their frozen assignments; the full `chematic-mol`
test suite passes, including the new V2000/V3000 fixtures.

## Provenance and reproduction

- Source corpus: `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`,
  SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
  Selection seed `20261018`, 100 per category, no excluded source rows.
- [Manifest and per-row reports](results/unmarked-3d-browser-holdout-2026-10-17.json.gz)
  contain the RDKit/COSMolKit source gate, CheMatic before and after, Open
  Babel before and after, the previous 300-row cross-check, versions, and
  input hashes. Compressed SHA-256:
  `ea74cbe8eeca34f6e1d7e037d258a6d610f7822a68f26c6505acc2636643afde`.
- Exact inputs: [V2000 SDF](results/unmarked-3d-browser-holdout-v2000-2026-10-17.sdf.gz)
  (uncompressed SHA-256 `93a5c6c3a09f92d2a25c6877c9203931c341f535ab1ce4e02752757aa9870a4d`)
  and [V3000 SDF](results/unmarked-3d-browser-holdout-v3000-2026-10-17.sdf.gz)
  (uncompressed SHA-256 `fec70b699c64dab5c8bb32f959db390abd77e7b724f9783e0771fa6168b76f9f`).
- CheMatic source-built example SHA-256:
  `5a12fdae0bf18b0a27dd092973b020085f7727cd01fecb8b8a3b51b6ddd1dfb2`.
  Open Babel binary SHA-256:
  `733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
  Generator and Open Babel comparator SHA-256 values are in the report.

From the repository root, with RDKit 2025.09.3 and Open Babel 3.2.1 installed:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline \
  -p chematic-mol --example sdf_v3000_accuracy_dump
gzip -dc validation/openbabel_comparison/results/unmarked-3d-browser-holdout-v2000-2026-10-17.sdf.gz > /tmp/browser-v2000.sdf
gzip -dc validation/openbabel_comparison/results/unmarked-3d-browser-holdout-v3000-2026-10-17.sdf.gz > /tmp/browser-v3000.sdf
gzip -dc validation/openbabel_comparison/results/unmarked-3d-browser-holdout-2026-10-17.json.gz > /tmp/browser-report.json
python - <<'PY'
import json
from pathlib import Path
report = json.loads(Path('/tmp/browser-report.json').read_text())
Path('/tmp/browser-manifest.json').write_text(json.dumps(report['manifest']))
PY
python validation/openbabel_comparison/check_unmarked_3d_stereo.py \
  --case /tmp/browser-manifest.json /tmp/browser-v2000.sdf \
  --case /tmp/browser-manifest.json /tmp/browser-v3000.sdf \
  --chematic target/debug/examples/sdf_v3000_accuracy_dump \
  --output-prefix /tmp/browser-check
```

Alternatively, regenerate the inputs with
`validation/cosmolkit_comparison/check_unmarked_3d_stereo.py` using the corpus,
seed, and category count above.

The RDKit source/reread gate is useful for a controlled reader comparison,
but it is not independent chemical ground truth or evidence of universal
accuracy superiority over Open Babel, RDKit, or COSMolKit.
