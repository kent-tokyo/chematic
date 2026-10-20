# Check supplied molecular geometry

[Open Babel #2651](https://github.com/openbabel/openbabel/issues/2651), still open
on October 19, 2026, asks how to inspect generated 3D coordinates and calculate
a batch acceptance rate. Finite coordinates alone cannot detect collapsed
bonds, overlapping atoms, or a contradicted stereocenter.

`chematic-3d::check_geometry` and Python `Mol.check_geometry` inspect **existing
Cartesian coordinates in Angstrom, in the molecule's original atom order**.
They reuse the native embedding bond-length model and declared-stereo verifier.
They do not generate, optimize, reorder, repair or replace coordinates.

```python
import chematic

mol = chematic.from_smiles("CC")
report = mol.check_geometry([[0, 0, 0], [10, 0, 0]])
assert report["bond_violations"][0]["bond"] == 0
assert report["bond_violations"][0]["distance"] == 10.0
```

## Report and thresholds

| Field | Meaning |
| --- | --- |
| `atom_count` | Original molecule and coordinate count |
| `bond_relative_tolerance` | Applied fractional deviation threshold; default 0.5 |
| `nonbonded_min_distance` | Applied absolute clash cutoff in Angstrom; default 1.2 |
| `bonds_checked` | Bonds with a supported reference length |
| `skipped_bonds` | Original indices of bonds without a supported length reference |
| `max_relative_bond_deviation` | Largest checked fractional deviation, or `None` without checked bonds |
| `bond_violations` | Original bond index, endpoint indices, distance, reference length, fractional deviation |
| `nonbonded_clashes` | Original atom-index pairs and measured distances, sorted by pair |
| `stereo` | Per-element tetrahedral/E/Z outcomes and declared/satisfied/violated/unevaluable counts |
| `stereo_groups_not_evaluated` | Enhanced stereo groups whose group semantics were not evaluated |

A bond is flagged when `abs(distance-reference)/reference` is **greater than**
the tolerance. Reference lengths are the existing covalent-radius sum multiplied
by 1.00 (single/directional), 0.87 (double), 0.78 (triple) or 0.93 (aromatic).
These are heuristic screening references, not experimentally fitted equilibrium
lengths. Dative, zero, quadruple and query bonds, and bonds touching wildcard
atoms, are listed as skipped. They do not count as checked or pass silently.

A clash is a graph-unconnected pair **strictly closer than** the cutoff. Every
graph bond excludes its endpoints from this comparison, including skipped bond
types. Pairs separated by two bonds and atoms in different fragments are
included. The default is a coarse absolute-distance screen, not an
element-specific van der Waals overlap calculation.

Each stereo row includes `status` (`satisfied`, `violated`, `unevaluable`),
`rejection_reason`, and its original `atom` or `bond` index. Unevaluable entries
remain unknown even when the violation count is zero. Undeclared stereo is not
inferred. Other stereo classes are unsupported; enhanced group semantics,
valence, angles, torsions, energies and global physical realizability are not
validated. There is deliberately no universal `valid` or `sound` result.
The existing embedding pipeline and its legacy soundness gate are unchanged.

## Batch decisions use an explicit policy

This example measures acceptance under one stated screen, with all input rows
retained. It does not measure the percentage of chemically valid molecules.

```python
packet = [
    ("CC", [[0, 0, 0], [1.54, 0, 0]]),
    ("CC", [[0, 0, 0], [10, 0, 0]]),
    ("CC", [[0, 0, 0], [0, 0, 0]]),
]
results = []
for index, (smiles, coordinates) in enumerate(packet):
    try:
        r = chematic.from_smiles(smiles).check_geometry(coordinates)
    except ValueError as error:
        results.append((index, False, str(error)))
        continue
    accepted = (
        not r["bond_violations"]
        and not r["nonbonded_clashes"]
        and not r["skipped_bonds"]
        and not r["stereo_groups_not_evaluated"]
        and r["stereo"]["n_violations"] == 0
        and r["stereo"]["n_unevaluable"] == 0
    )
    results.append((index, accepted, r))
assert len(results) == len(packet)
assert [row[1] for row in results] == [True, False, False]
acceptance_rate = sum(row[1] for row in results) / len(results)
```

## Rust and input limits

```rust
use chematic_3d::{check_geometry, Coords3D, GeometryCheckConfig, Point3};

let molecule = chematic_smiles::parse("CC").unwrap();
let coordinates = Coords3D {
    points: vec![Point3::zero(), Point3::new(10.0, 0.0, 0.0)],
};
let report = check_geometry(&molecule, &coordinates, &GeometryCheckConfig::default()).unwrap();
assert_eq!(report.bond_violations.len(), 1);
```

Requires 1–2048 atoms and exactly one finite coordinate triple per atom, including
explicit hydrogen atoms. Each component must be within ±1,000,000 Angstrom to
bound stereo arithmetic. The bond tolerance must be finite and in `[0, 1)`;
the clash cutoff must be finite and positive. Rust returns `GeometryCheckError`;
Python raises `ValueError` for invalid values/counts and may raise `TypeError`
for wrong types. Python consumes at most `atom_count + 1` outer iterable items
and propagates source errors. Runtime and worst-case clash output are O(N²).

Validation uses original analytic controls for collapsed/stretched bonds,
nonbonded coincidences, bond-order scales, mirror/planar stereochemistry,
unsupported references, enhanced groups, input boundaries, rigid transforms,
and unchanged molecule/coordinate data. No upstream fixture corpus is imported.
Independent checks compare 360 indexed-distance reports with Python/RDKit graph
references and 192 nondegenerate tetrahedral reports with RDKit 2026.09.1
coordinate stereoperception. These measure diagnostic agreement, not physical
conformer quality.
See [the dated validation record](https://github.com/kent-tokyo/chematic/tree/main/validation/results/geometry-check-2026-10-19).
