# ETKDGv3/MMFF94 conformer search against COSMolKit 0.5.0

This is a source-candidate comparison of CheMatic's 3D quality profile with
the [COSMolKit 0.5.0 Python package](https://pypi.org/project/cosmolkit/).
It measures the minimum **converged** MMFF94 energy in each generated ensemble,
not experimental conformer accuracy. The CheMatic candidate starts from
`origin/main` at `99c5c07e` and adds a conditional fourth seed batch: it runs
only when the third batch lowers the best finite energy by more than
0.01 kcal/mol. The existing stereo check on optional MMFF retries is
unchanged. The compatibility embedding API and its seeds are unchanged.

## Development set

The first 24 rows of the committed Tier B manifest were used to find and fix
the residual search gap. Each arm uses explicit hydrogens, ETKDGv3, seed 42,
one thread, 24 conformers per batch, MMFF94, nonbonded threshold 100 Å, and
up to 1,000 optimizer iterations. COSMolKit and CheMatic's compatibility arm
each generate 24 candidates. The quality arm generates 48–96 candidates and
may also retry unfinished minima, so its energy comparison has a larger
compute budget.

| Arm versus COSMolKit | CheMatic lower by >0.01 | COSMolKit lower by >0.01 | Within 0.01 kcal/mol |
| --- | ---: | ---: | ---: |
| CheMatic compatibility, 24 candidates | 0 | 4 | 20 |
| CheMatic quality, 48–96 candidates | **4** | **0** | 20 |

All 24 COSMolKit and CheMatic output ensembles had 24 converged conformers.
The four resolved cases are `chembl_tier_b_0009`, `_0020`, `_0022`, and
`_0023`. The newly added fourth batch changes `_0020` from −19.561605 to
−20.314056 kcal/mol; COSMolKit yields −19.818846 kcal/mol. The other three
cases were resolved by the previous quality profile. The
[row-level record](results/3d-mmff-conformer-dev-2026-10-24.json) includes
both source profiles, COSMolKit values and statuses, all 24 SMILES, package
versions, binary hashes and the input-manifest hash.
The [source result stream](results/3d-mmff-conformer-dev-source-2026-10-24.jsonl)
matches the SHA-256 recorded in the row-level record.

The metric is MMFF94 energy under a fixed graph and force-field model. A
lower minimum can indicate a better search of that model's conformational
space, but it is not evidence of a more realistic structure. The quality
profile's larger candidate budget means its result is an accuracy-first
workflow comparison, not a same-work algorithm win. The compatibility arm
remains worse in four cases under equal candidate and iteration budgets.

## Disjoint checks

Two checks use rows outside the 24 development molecules. The first is the
next six consecutive Tier B rows (manifest indices 24–29). A longer run was
stopped when the following macrocycle became expensive; the six completed
rows are retained as the contiguous check. The second selects the eight
smallest heavy-atom counts from manifest indices 24–47, breaking ties by
name: indices 33, 41, 40, 46, 37, 42, 39, and 43. Selection used molecule
size alone, before inspecting engine outcomes. The two checks are disjoint;
together they cover 14 additional molecules, with a size bias in the second.

| Check and arm | CheMatic lower by >0.01 | COSMolKit lower by >0.01 | Within 0.01 kcal/mol |
| --- | ---: | ---: | ---: |
| Next six, compatibility | 3 | 1 | 2 |
| Next six, quality | 5 | 0 | 1 |
| Eight smallest, compatibility | 6 | 1 | 1 |
| Eight smallest, quality | 7 | 0 | 1 |

Every compared arm returned 24 converged conformers for these rows. The
[next-six comparison](results/3d-mmff-conformer-next-six-2026-10-24.json)
and [source stream](results/3d-mmff-conformer-next-six-source-2026-10-24.jsonl),
and the [size-selected comparison](results/3d-mmff-conformer-size-holdout-2026-10-24.json)
and [source stream](results/3d-mmff-conformer-size-holdout-source-2026-10-24.jsonl)
retain row-level results and hashes. These small checks support the quality
profile on the selected workload. They do not establish general conformer
accuracy or a same-compute advantage.

## Reproduce

Install pinned COSMolKit 0.5.0 and RDKit 2026.09.1 in the Python interpreter,
then run from the repository root:

```sh
cargo run -p chematic-3d --example cosmolkit_quality_probe --offline -- \
  validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json 24 \
  > /tmp/chematic-3d-source.jsonl
python validation/cosmolkit_comparison/compare_3d_conformer_quality.py \
  --manifest validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json \
  --source-results /tmp/chematic-3d-source.jsonl \
  --source-bin target/debug/examples/cosmolkit_quality_probe \
  --output /tmp/chematic-cosmolkit-3d.json
```

The comparison script requires complete source coverage, checks COSMolKit's
graph and explicit-H atom count against RDKit, and rejects nonfinite converged
energies. The score counts only rows where both arms have a converged result.
The `0.01 kcal/mol` threshold excludes tiny floating-point differences.
