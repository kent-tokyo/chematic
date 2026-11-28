# Common-coordinate MMFF development comparison

The 20 ordinary development inputs in `chematic-ordinary20-forcefields-input.json`
were selected independently of force-field results from the reserved COD
development partition. They contain no boron or five-neighbor center and at
most 25 heavy atoms. RDKit 2026.09.2 generated explicit-hydrogen ETKDGv3
coordinates with seed 42 on macOS arm64, Python 3.13.6. Both engines evaluated
those same coordinates. The exact coordinates, per-atom errors and source
SMILES are retained in `validation/results/2026-11-26-3d-wheel-pilot/`.
The chematic optimized source library SHA-256 is
`af8b3a615fc163232d6d277f11bf5829d7334b94eef5adf986dfad5f37af3893`.
The reference wheel SHA-256 is
`3860ba53505f27cb24d05383027f78e3f82aefd6caf193c217a94d876d91b023`.

Both engines parameterized 19 of the 20 molecules for MMFF94 and MMFF94s;
both declined the remaining molecule. In the 19 common cases, all numeric
atom types matched and the largest absolute partial-charge difference was
1.12e-16 elementary charges. For both variants, the maximum total-energy
difference was 8.53e-14 kcal/mol and the largest componentwise gradient
difference was 3.59e-7 kcal/mol/Å. The largest difference occurred in a
torsion term. Each of the seven term classes was compared separately; the
largest per-term energy difference was 5.96e-14 kcal/mol and the largest
non-torsion componentwise gradient difference was 3.13e-13 kcal/mol/Å.

UFF was parameterizable in all 20 cases; its largest total-energy difference
was 2.28e-13 kcal/mol and its largest componentwise gradient difference was
3.46e-7 kcal/mol/Å. Six evenly spaced gradient components per molecule were
also checked by central finite differences with step 1e-5 Å. These small,
previously unscored development inputs do not establish general force-field
parity or coverage of boron and selenium.
