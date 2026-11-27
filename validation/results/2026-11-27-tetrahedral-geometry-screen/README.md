# Tetrahedral geometry screen: exposed development comparison

The input is the 1,345-row exposed COD development partition stored at
`../2026-11-27-selenium-crystal-stereo/source-development-input.json` (SHA-256
`e3c2cf4aa5c7c5740fd3d23c670e26535515b97f16a02ee2bf387661513726a1`). No
reserved holdout was used. The script skips the 680 input SMILES without an
`@` marker. It parses the remaining 665 in chematic and RDKit, verifies matching
atom symbols, atomic numbers, charges and bond-index pairs, and generates one
chematic ETKDG geometry per molecule. No experiment-level parse or mapping
errors occurred.

For every RDKit-recognized declared tetrahedral tag, the script compares the
input's RDKit chiral tag with the tag RDKit reassigns from the generated
coordinates. It compares that agreement with chematic's raw declared-stereo
result and the opt-in geometry-screened result. Direct chiral tags are
compared on the same RDKit molecular graph, avoiding CIP-priority changes caused
by assignments at other centers. Tags that RDKit did not reassign as
`CHI_TETRAHEDRAL_CW` or `CHI_TETRAHEDRAL_CCW`, and chematic's already
unevaluable statuses, are retained in the per-center output but excluded from
the comparable-decision totals. This is an independent implementation check,
not a crystal-coordinate or energy reference.

Among 2,094 declared tetrahedral centers, 2,055 were comparable by that rule;
the other 39 included 2 raw unevaluable centers and 37 for which RDKit's 3D
result was not a tetrahedral tag. In the 2,055 comparable decisions:

| Raw chematic vs RDKit | Screen retained decision | Screen made unevaluable |
|---|---:|---:|
| Agreed | 975 | 890 |
| Disagreed | 0 | 190 |

Thus the screen catches all 190 disagreements in this exposed fixed sample,
while withholding 890 decisions that agreed. The retained 975 decisions
agreed in this sample. These counts are not a general accuracy estimate or a
claim that all distorted or erroneous cases will be detected. The test inputs
are generated geometries, and the heuristic may withhold chemically valid
centers. It does not screen E/Z geometry, unsupported stereo, bond lengths,
energies, or steric clashes.

The candidate was a local macOS arm64 debug wheel from commit `b746d445` plus
the accompanying source changes, built using CPython 3.13.6 and
`maturin build --locked`; SHA-256
`58ff82ddb167bdec2036d8b275099d7a046d58c5ca08fc4b2f1ec9ebbe5163c9`. RDKit was
2026.09.2. This is not a published wheel. Reproduce with both installed:

```sh
python validation/results/2026-11-27-tetrahedral-geometry-screen/run.py > results.json
```

The screen requires every center-in-ligand-tetrahedron barycentric weight to
exceed 1e-6 and all ligand-direction angles to lie in 45–165 degrees, including
reconstructed implicit ligand positions. The normal verifier remains the
unchanged default. A separately retained prototype fixture reproduces one
RDKit disagreement and verifies that the optional screen reports it as
unevaluable.
