# Checked conformer ensembles

Experimental Python methods for generating and MMFF94-minimizing ETKDG
conformers with additional checks are available in v1.2.4 and later. They are
opt-in and do not establish output equivalence to RDKit.

```python
import chematic

mol = chematic.from_smiles("C[C@H](O)F").add_hydrogens()
conformers = mol.etkdg_embed_multiple_confs_mmff_mixed_checked(
    20, random_seed=42, max_mmff_iters=200, repair_mmff_iters=400
)

for status, energy, coordinates in conformers:
    if status == 0:
        print(energy, coordinates)
```

The returned coordinates are in atom order and include the explicit
hydrogens. Each tuple contains the minimization status, MMFF94 energy in
kcal/mol, and coordinates in angstroms. Status zero requires a finite energy
and coordinates and a final MMFF gradient with no component greater than
0.01 kcal/mol/angstrom in magnitude. A nonzero status means the candidate
was not confirmed converged; it can still be returned. Some attempts may
fail or be discarded, so the returned list can be shorter than the requested
count or empty. Invalid options or force-field setup failures raise
`ValueError` in the Python binding.

The mixed method divides the requested attempts between the compatibility
and corrected ETKDG torsion gradients, using the first half for compatibility
and the remainder for corrected. The two modes share the same attempt seed
schedule and MMFF94 minimization. For a compatibility-gradient ensemble
without that split, use
`mol.etkdg_embed_multiple_confs_mmff_compatible_checked(...)` with the same
arguments. Both methods apply the source's declared-stereo checks. If a
minimized structure fails a check, the implementation can retain the original
embedded structure with a nonzero status, or omit the candidate. These checks
do not establish correctness for every possible stereochemical specification.

`repair_mmff_iters` defaults to zero. A positive value allows that many extra
MMFF iterations when attempting to repair a near-linear three-coordinate
nitrogen; a repair is retained only if its checks pass. The example uses 400
to match the ongoing evaluation, not as a universally optimal setting.
The [evaluation protocol](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-platinum-holdout-protocol.md)
and [macrocycle evaluation](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-macrocycle-pdb-evaluation.md)
record the tested candidates, datasets, screens and results. The 0.2-second
acceptance deadline used in timed evaluations does not interrupt an attempt.
