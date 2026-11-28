# Reconstructed explicit hydrogen and near-planar carbon screening

The previous optional screen applied a normalized three-bond-direction volume
threshold of 0.15 to a carbon with three graph neighbors and one implicit H,
but not to the same center after that H was made explicit. A fixed synthetic
fixture in `crates/chematic-py/tests/test_stereo_3d_diagnostics.py` reproduces
the difference: after reconstructing the H along the implicit-H phantom
direction, the implicit and explicit forms used to receive different screen
results. RDKit 2026.09.2 leaves the center unassigned in both forms on that
fixture. This one fixture is not a general RDKit accuracy measurement.

The candidate also applies the existing 0.15 threshold to carbon with exactly
three non-hydrogen neighbors and one explicit H **when the explicit H direction
matches the reconstructed implicit-H direction**. The matching condition uses
a unit-vector dot product of at least `1 - 1e-8`; an independent explicit H
position is not subject to the added volume test. The normal non-screened
stereo verifier remains unchanged. This heuristic does not prove that other
explicit-H directions are safe or that a passing structure is physically valid.

## Exposed development consistency

`development.py` reads the 1,345-row exposed COD development input used in
`../2026-11-27-tetrahedral-geometry-screen/`. It skips 680 SMILES without `@`
and finds no target in 97 others. In the remaining 568 molecules it selects
reported carbon stereo centers with three non-hydrogen graph neighbors and one
implicit H. It generates one chematic ETKDG structure, adds explicit H atoms,
checks preservation of the original atom ordering and H-parent bonds, and
places each target H along the same reconstructed direction used by the
implicit verifier. Other added hydrogen coordinates are placeholders; this is a
local representation-consistency experiment, not a complete explicit-H
conformer or crystal-geometry assessment. There were no processing errors or
skipped target centers.

Across 1,471 target centers, the complete screened status-and-reason entries
matched between implicit and reconstructed-explicit H in 1,451 cases before
and all 1,471 after. The 20 changed explicit-H results were 11 `satisfied` and
9 `violated` decisions that became `unevaluable`, matching the implicit form.
The recorded implicit screened decisions were unchanged. The per-center before and
after outputs are retained. No reserved molecular holdout was used.

```sh
python validation/results/2026-11-27-explicit-h-planarity/development.py > output.json
```

## Synthetic explicit hydrogen comparison

`synthetic.py` generates explicit-H `F[C@]([H])(Cl)Br` geometries with RDKit
used only to reassign the 3D chiral tag. Each scenario contains 30,000 seeded
geometries of one molecule, not independent molecular data. In one scenario H
is placed along the reconstructed implicit direction; in the other its
position is sampled independently. The script saves the complete contingency
counts and a per-case classification-and-screen-status sequence digest.

| Synthetic case | RDKit unassigned but screen decided, before → after | RDKit-agreeing cases decided, before → after |
|---|---:|---:|
| Reconstructed H | 4 → 0 | 13,525 → 13,376 |
| Independent H | 0 → 0 | 2,373 → 2,373 |

The full per-case status sequence digest for independent H was identical before
and after in this sample. For reconstructed H, the screen additionally withheld
149 RDKit-agreeing decisions. In the independent-H sample there were 9,386
raw chematic/RDKit disagreements, all already unevaluable under both versions
of the optional screen. These artificial distributions and counters do not
establish general agreement or physical conformer quality.

```sh
python validation/results/2026-11-27-explicit-h-planarity/synthetic.py > output.json
```

The baseline is commit `9b6438c7`, local debug wheel SHA-256
`dd1fc9a0527dd3ccbe7add431e646c44d9ab9fed143f6dc27b6af55a4c16ec50`.
The candidate is that commit plus the accompanying code changes, local debug
wheel SHA-256 `d0c655bd38fcc6d133a7ffa5c933caaf044421ea588eaa47f67ef2bef0d3aaa0`.
Both are macOS arm64 CPython 3.13.6 source builds, not published wheels;
synthetic comparisons used NumPy 2.5.3 and RDKit 2026.09.2.
