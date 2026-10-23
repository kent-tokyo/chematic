# Bounded torsion search after ETKDG/MMFF94

Date: 2026-10-24. This source candidate includes `origin/main` at
`02f5a37e` and adds an opt-in quality profile. The comparator is RDKit
2026.09.1 on macOS arm64 with Python 3.13.6. The fixed workloads and source
binary hashes are in the [Tier A](../validation/results/etkdg-torsion-search-tier-a-2026-10-24.json),
[Tier B](../validation/results/etkdg-torsion-search-tier-b-2026-10-24.json),
and [ChEMBL](../validation/results/etkdg-torsion-search-chembl-2026-10-24.json)
row records.

The existing quality profile remains unchanged. The new
`etkdg_embed_multiple_confs_mmff_torsion_search` API takes its lowest
converged MMFF94 conformer, rotates each of at most twelve strict acyclic
rotatable bonds by ±120°, and reoptimizes each trial for at most 200 MMFF94
iterations in this measurement. A trial is accepted only if it converges,
has finite coordinates, lowers energy by more than 1e-6 kcal/mol, and keeps
declared stereochemistry. At most 24 extra minimizations are made. This
targets local conformational basins while retaining the original result
when no valid improvement exists.

Both engines use explicit hydrogens, ETKDGv3, seed 42, MMFF94, and up to
200 iterations per generated candidate. CheMatic's base quality profile
uses up to 80 candidates in four 20-seed batches and may make three bounded
retries; it returns at most 20. The fourth batch runs only when the third
lowers the best energy by more than 0.01 kcal/mol. RDKit generates 60
candidates and returns all of them. The new torsion arm has a larger
optimization budget. It is a
search-quality comparison, **not** an equal-compute algorithm comparison.
The metric is the lowest **converged** MMFF94 energy; 1e-6 kcal/mol is the
comparison tolerance.

| Input set | Comparable | CheMatic lower | RDKit lower | Within tolerance | No RDKit converged result |
| --- | ---: | ---: | ---: | ---: | ---: |
| Tier A, 9 eligible rows | 8 | 1 | 0 | 7 | 1 |
| Tier B, 21 prior development rows | 20 | 11 | 0 | 9 | 1 |
| ChEMBL, 12 hash-selected rows | 11 | 4 | 0 | 7 | 1 |

The torsion arm lowered its own base quality-profile minimum on 1/9 Tier A,
5/21 Tier B, and 2/12 ChEMBL rows. The ChEMBL rows were selected before
scoring either engine: from the committed 5,000-line corpus, 851 molecules
met the predeclared single-fragment, at most 25-heavy-atom, 3–8-rotatable-bond,
MMFF-parameter, and no-high-order-ETKDG-torsion filters; the 12 smallest
SHA-256 hashes of SMILES were used, excluding exact Tier A/B SMILES. The
Tier A/B manifests and ChEMBL corpus have all been used in the repository;
the ChEMBL sample is disjoint from these two manifest sets, not a sealed
external test set.

The known `chembl_tier_b_0022` gap closes: its lowest converged energy moves
from −52.729760 to −53.719974 kcal/mol. RDKit re-scoring of the resulting
CheMatic coordinates returns the same energy within floating-point error.
The largest absolute same-coordinate re-score difference among the 42 rows
is below 1e-12 kcal/mol. The comparison rejects changes to declared
tetrahedral or double-bond stereo on each returned best conformer. For three
rows, RDKit yielded no converged result at 200 iterations, so no energy rank
is assigned.

## Reproduce

With RDKit 2026.09.1 installed in the Python interpreter:

```sh
cargo build --locked -p chematic-3d --example etkdg_torsion_search_probe
for tier in a b chembl; do
  python3 scripts/benchmarks/etkdg_torsion_search_probe.py \
    --binary target/debug/examples/etkdg_torsion_search_probe \
    --tier "$tier" --output "/tmp/etkdg-torsion-search-$tier.json"
done
```

The scorer checks atom-indexed graph and bond orders, independently scores
both CheMatic minima with RDKit MMFF94, verifies declared stereo, and stores
the exact source binary, script, corpus, and torsion-source hashes.

The experiment does not measure experimental geometry, conformer diversity,
population coverage, wall time, memory, equal work, or published-wheel
behavior. The force-field formula is unchanged, so this does not establish
better physical MMFF accuracy or broad 3D superiority over RDKit.
