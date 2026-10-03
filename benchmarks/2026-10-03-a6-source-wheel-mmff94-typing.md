# A6 MMFF94 atom typing: opt-in source-wheel gate

This checks a **locally built source wheel**, not the published v1.0.31
package. The wheel SHA-256 is
`41a360a1b5eb9404515d6395f2c55dd429961cf7891810aea78fc0490d36324b`.
The input is the pinned 10,000-SMILES corpus (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`),
compared to RDKit **2026.03.6** on macOS arm64/Python 3.13. The source
contains the bounded compact-ring MMFF aromatic-state propagation. Its
native Rust source differential is recorded [separately](2026-10-03-a6-mmff94-aromatic-propagation-source.md).

| Atom-typing outcome | Published v1.0.31 wheel | Local source wheel |
|---|---:|---:|
| Comparable molecules | 9,774 | 9,774 |
| RDKit unsupported / CheMatic errors | 204 / 22 | 204 / 22 |
| Heavy atoms compared / differing | 214,990 / 2,908 | 214,990 / **451** |
| Explicit H compared / differing | 190,737 / 56 | 190,737 / **8** |
| RDKit carbon type 37 / CheMatic type 2 | 1,322 | **16** |

The 16 remaining type-37/type-2 carbons are in two rows that each lack one
RDKit ring in CheMatic's symmetrized ring view. This is a diagnosis, not a
reason to force an aromatic flag on those atoms. The aggregate improvement
does not establish per-row non-regression, parameter parity, force-field
convergence, geometry, stereo safety, clash absence, equal-coordinate energy,
or conformer quality. The [published quality gate](2026-10-03-a6-published-v1031-mmff94-quality.md)
still reports only 100/265 converged at the shipped 200-iteration limit.

The local gated census JSON has SHA-256
`ca6c1d84c76002d8f72aa37668451147d0343e3c93e6e398dc64452b489427a4`.
It is a temporary measurement, not a checked-in package artifact. Recreate it
using a freshly installed source wheel and RDKit 2026.03.6:

```sh
python scripts/mmff94_atom_type_census.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --output validation/results/mmff94-type-source-wheel-10k-local.json \
  --wheel dist/hba-wheel/chematic-*.whl --module-root .venv \
  --type37-context --source-wheel-gate
```

The gate checks import provenance and the corpus hash, records the wheel hash,
and checks all row-status
denominators, the heavy/H difference totals, and the type-37/type-2 bucket.
Post-merge PR [#736](https://github.com/kent-tokyo/chematic/pull/736), run
`37102043846`, built release-profile source wheels on Linux and macOS arm64.
Both retained `rdkit-hba-source-wheel-5k` artifacts were downloaded and read:
each reports 9,774 comparable rows, 204 RDKit unsupported, 22 CheMatic
typing errors, **451** differing heavy types and **8** differing H types.
All 76 CI checks passed (two intentionally skipped). These are source-wheel
results, not published-package quality or proof of per-row non-regression.
