#!/usr/bin/env python3
"""Independent indexed-distance and RDKit stereo checks on original controls.

Requires the candidate native extension and RDKit 2026.09.1 on PYTHONPATH.
No external molecule files or generated coordinates are imported.
"""

import itertools
import json
import math
import random

import chematic
from rdkit import Chem, rdBase


def distance_checks() -> int:
    rng = random.Random(2651)
    smiles = [
        "CC",
        "CCO",
        "C1CCCCC1",
        "c1ccccc1",
        "C/C=C/C",
        "C/C=C\\C",
        "FC(Cl)(Br)I",
        "C#N",
        "CC.C",
        "[H]O[H]",
    ]
    params = Chem.SmilesParserParams()
    params.removeHs = False
    count = 0
    orders = {
        "SINGLE": 1.0,
        "UP": 1.0,
        "DOWN": 1.0,
        "DOUBLE": 2.0,
        "TRIPLE": 3.0,
        "AROMATIC": 1.5,
    }
    for text in smiles:
        native = chematic.from_smiles(text)
        reference = Chem.MolFromSmiles(text, params)
        assert reference is not None
        assert [a[1] for a in native.atom_table] == [
            a.GetAtomicNum() for a in reference.GetAtoms()
        ]
        native_edges = {
            (min(a, b), max(a, b)): orders[o] for a, b, o, _ in native.bond_table
        }
        edges = {
            (
                min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            ): b.GetBondTypeAsDouble()
            for b in reference.GetBonds()
        }
        assert native_edges == edges, text
        n = reference.GetNumAtoms()
        for attempt in range(12):
            rows = [[rng.uniform(-2.5, 2.5) for _ in range(3)] for _ in range(n)]
            if attempt == 0:
                rows = [[0.0, 0.0, 0.0] for _ in range(n)]
            for cutoff in (0.8, 1.2, 1.6):
                r = native.check_geometry(rows, nonbonded_min_distance=cutoff)
                expected = {
                    (i, j): math.dist(rows[i], rows[j])
                    for i in range(n)
                    for j in range(i + 1, n)
                    if (i, j) not in edges and math.dist(rows[i], rows[j]) < cutoff
                }
                actual = {
                    tuple(v["atoms"]): v["distance"] for v in r["nonbonded_clashes"]
                }
                assert actual.keys() == expected.keys(), (text, attempt, cutoff)
                assert all(
                    math.isclose(v, expected[p], abs_tol=1e-12)
                    for p, v in actual.items()
                )
                assert not r["skipped_bonds"]
                assert r["bonds_checked"] == len(edges)
                for v in r["bond_violations"]:
                    a, b = v["atoms"]
                    assert math.isclose(
                        v["distance"], math.dist(rows[a], rows[b]), abs_tol=1e-12
                    )
                count += 1
    return count


def stereo_checks() -> int:
    base = [
        [1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ]
    count = 0
    for text in ("F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I", "F[C@H](Cl)Br", "F[C@@H](Cl)Br"):
        reference = Chem.MolFromSmiles(text)
        expected_cip = reference.GetAtomWithIdx(1).GetProp("_CIPCode")
        native = chematic.from_smiles(text)
        for mirror in (False, True):
            for axes in itertools.permutations(range(3)):
                inversions = sum(
                    axes[i] > axes[j] for i in range(3) for j in range(i + 1, 3)
                )
                for signs in itertools.product((-1, 1), repeat=3):
                    if (-1) ** inversions * math.prod(signs) != 1:
                        continue
                    rows = [
                        [(-1 if mirror and a == 0 else 1) * p[a] for a in range(3)]
                        for p in base[: reference.GetNumAtoms()]
                    ]
                    rows = [
                        [signs[a] * p[axes[a]] + (20, -10, 5)[a] for a in range(3)]
                        for p in rows
                    ]
                    candidate = Chem.Mol(reference)
                    Chem.RemoveStereochemistry(candidate)
                    conf = Chem.Conformer(len(rows))
                    for i, row in enumerate(rows):
                        conf.SetAtomPosition(i, row)
                    candidate.AddConformer(conf)
                    Chem.AssignAtomChiralTagsFromStructure(
                        candidate, confId=0, replaceExistingTags=True
                    )
                    Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
                    actual_cip = candidate.GetAtomWithIdx(1).GetProp("_CIPCode")
                    r = native.check_geometry(rows)["stereo"]
                    assert r["n_declared"] == 1 and r["n_unevaluable"] == 0
                    assert r["n_violations"] == int(actual_cip != expected_cip), (
                        text,
                        rows,
                        r,
                    )
                    count += 1
    return count


if __name__ == "__main__":
    assert rdBase.rdkitVersion == "2026.09.1", rdBase.rdkitVersion
    print(
        json.dumps(
            {
                "rdkit": rdBase.rdkitVersion,
                "native_version": chematic.__version__,
                "indexed_distance_reports": distance_checks(),
                "tetrahedral_stereo_reports": stereo_checks(),
            },
            indent=2,
        )
    )
