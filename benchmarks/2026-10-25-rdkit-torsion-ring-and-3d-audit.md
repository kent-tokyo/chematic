# RDKit torsion ring and ETKDG/MMFF94 audit

Date: 2026-10-25. This source audit started from `origin/main` at
`eb850d4b`. The fixed comparator is RDKit 2026.09.1 on macOS arm64.

## Topological torsion fingerprints

The source previously warned that asymmetrically substituted cyclopropanes
could add more than one three-membered-ring closure. That warning confused
RDKit's extra four-distinct-atom walks with closures. In RDKit's
`FindAllPathsOfLengthN(mol, 4, useBonds=False)`, the tested cyclopropanes
have one closure per ring. CheMatic already generated that closure, so no
fingerprint algorithm change was needed.

The new [ring audit](../validation/results/rdkit-torsion-triangle-parity-2026-10-25.json)
enumerates 7³ substitutions (none, C, N, O, F, Cl, Br at each carbon) and
three rooted SMILES spellings, deduplicated to 343 inputs. All **343/343**
match RDKit's sparse topological-torsion counts and 2048-bit fingerprint
exactly. The result retires that specific documented limitation. It does not
prove universal parity, assess chiral torsion options, or show higher
accuracy than RDKit.

## 3D conformer search

The [12-row ChEMBL report](../validation/results/etkdg-torsion-search-chembl-hash36-47-2026-10-25.json)
uses the next preselected SHA-256 SMILES positions 36–47 among 851 eligible
rows, disjoint from the preceding 0–35 comparison slices. The existing
two-pass CheMatic torsion profile gave a lower converged MMFF94 energy on
5/12 inputs and tied within 1e-6 kcal/mol on 7/12; RDKit gave no lower
energy in this slice. Torsion refinement improved CheMatic's own base
quality-profile minimum on 3/12. Independent RDKit re-scoring of CheMatic's
coordinates differed by at most 4.7e-14 kcal/mol. The scorer also checked
the atom-indexed graph and declared tetrahedral/double-bond stereo.

Both engines used ETKDGv3 seed 42 and MMFF94 with a 200-iteration limit per
minimization. CheMatic's quality profile can search up to four 20-candidate
seed batches and perform bounded retries, then the torsion profile can add
up to 48 minimizations. RDKit generated 60 candidates. Thus this is a
**minimum-energy search comparison with unequal work**, not evidence of
superior physical force-field accuracy, equal-compute performance, conformer
diversity, experimental geometry agreement, or broad RDKit superiority.
These ChEMBL rows are repository development data, not a sealed external
holdout. A trial of an additional 180° torsion move on the preceding 24–35
slice found no lower minima and was dropped because it increased work.

## Reproduce

```sh
cargo build --locked -p chematic-fp --example rdkit_torsion_ring_probe
python3 scripts/validation/rdkit_torsion_triangle_parity.py \
  --binary target/debug/examples/rdkit_torsion_ring_probe \
  --output /tmp/rdkit-torsion-triangle-parity.json

cargo build --locked -p chematic-3d --example etkdg_torsion_search_probe
python3 scripts/benchmarks/etkdg_torsion_search_probe.py \
  --binary target/debug/examples/etkdg_torsion_search_probe \
  --tier chembl --chembl-start 36 \
  --output /tmp/etkdg-torsion-chembl-hash36-47.json
```

Use a Python environment with RDKit 2026.09.1. Both scorers enforce that
version and record source/input hashes. The ring audit follows the official
[RDKit fingerprint generator](https://rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html)
and [RDKit Book](https://github.com/rdkit/rdkit/blob/master/Docs/Book/RDKit_Book.rst)
definitions; the path-enumeration finding is an inference from pinned RDKit
outputs, not a claim about every molecule.
