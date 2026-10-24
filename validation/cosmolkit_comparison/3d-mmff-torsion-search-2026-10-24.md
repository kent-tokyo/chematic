# Bounded MMFF94 torsion search against COSMolKit 0.5.0

The source candidate at `39e6a460` adds an opt-in torsion search around
CheMatic's best converged quality-profile conformer. This follow-up tests it
against the [COSMolKit 0.5.0 package](https://pypi.org/project/cosmolkit/).
The metric is the lowest **converged** MMFF94 energy, with differences of at
most 0.01 kcal/mol classified as equivalent. Lower force-field energy
measures the search of this model, not experimental conformer accuracy.

## Selection and conditions

The selection starts at index 48 of the committed Tier B manifest, beyond the
38 molecules in the [earlier comparison](3d-mmff-conformer-2026-10-24.md).
It keeps single-fragment molecules with 8–25 heavy atoms, 3–8 RDKit
rotatable bonds and MMFF94 parameters, excluding the manifest's macrocycle
category. The 35 eligible rows were ordered by SHA-256 of their source SMILES;
the first 12 were selected before scoring either engine. This is a small,
structure-filtered check, not a representative accuracy sample. The selected
manifest indices are 192, 73, 60, 159, 174, 53, 146, 126, 95, 169, 68, and
66. The [row-level result](results/3d-mmff-torsion-hash12-2026-10-24.json)
records the selection rule and the manifest SHA-256.

All arms use explicit H, ETKDGv3 seed 42, 24 conformers per seed batch,
MMFF94, a 100 Å nonbonded threshold, and at most 1,000 optimizer iterations.
COSMolKit uses one thread and 24 candidates. CheMatic's opt-in quality profile
can draw 48–96 candidates and retry unfinished minima. The torsion arm starts
from that profile and makes at most 24 further MMFF minimization attempts,
while returning the best 24 candidates. This is an accuracy-first workflow
comparison with a larger CheMatic compute budget. Runtime was not compared.

## Results

| CheMatic arm versus COSMolKit | CheMatic lower by >0.01 | COSMolKit lower by >0.01 | Within 0.01 kcal/mol |
| --- | ---: | ---: | ---: |
| Quality profile | 6 | 1 | 5 |
| Quality plus torsion search | **7** | **0** | 5 |

Each engine returned 24 converged conformers on every row. Torsion search
lowered CheMatic's best energy on five rows. On the quality profile's sole
loss, `chembl_tier_b_0066`, CheMatic moved from 110.549510 to 110.375326
kcal/mol; COSMolKit gave 110.375328 kcal/mol. The difference after torsion
search is about 0.0000023 kcal/mol, within the comparison threshold. The
other changed classification is `_0169`, which moves from equivalent to
CheMatic lower.

The comparator validates indexed atom elements and bond topology/order for
both engines against RDKit, and COSMolKit's formal charges and isotopes. It
checks that specified tetrahedral and double-bond stereochemistry survive in
the selected coordinates. RDKit independently re-scores both engines' best
converged coordinates with the same MMFF94 settings. The largest absolute
energy discrepancy from either engine's reported value is below
5×10⁻¹⁴ kcal/mol. The [source result stream](results/3d-mmff-torsion-hash12-source-2026-10-24.jsonl)
contains CheMatic's returned coordinates; its SHA-256 and the executable,
script, torsion-source, and COSMolKit extension hashes are in the row-level
result.

These 12 rows show that this opt-in search closes the observed COSMolKit
energy gap on the selected workload. They do not establish a same-compute
advantage, a lower-energy result for every molecule, or superior physical
accuracy. Coverage, force-field typing, speed, and conformer diversity need
separate comparisons.

## Reproduce

Install COSMolKit 0.5.0 and RDKit 2026.09.1 in the same Python interpreter.
From the repository root, build the source example and run the complete
comparison:

```sh
cargo build -p chematic-3d --example etkdg_torsion_search_probe
python validation/cosmolkit_comparison/compare_3d_torsion_search.py \
  --binary target/debug/examples/etkdg_torsion_search_probe \
  --source-output /tmp/chematic-cosmolkit-torsion-source.jsonl \
  --output /tmp/chematic-cosmolkit-torsion.json
```

The script pins both package versions, applies the selection rule, runs the
CheMatic example and COSMolKit API, validates graph/stereo and coordinates,
and writes row-level counts and SHA-256 hashes. For re-scoring a retained
CheMatic source stream, pass `--source-results PATH`; that mode avoids
rerunning the CheMatic conformer search.
