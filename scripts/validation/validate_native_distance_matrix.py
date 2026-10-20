#!/usr/bin/env python3
"""Pinned, independent heavy-subgraph distance oracle for RDKit #9647 scope."""
from __future__ import annotations

import hashlib
import itertools
import json
import pathlib
import sys

import chematic
from rdkit import Chem, rdBase

ROOT = pathlib.Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "validation/rdkit_issues/native-distance/fixtures.json"
OUTPUT = ROOT / "validation/results/native-distance-matrix-2026-10-15.json"


def digest(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(data).hexdigest()


def packet():
    fixture = json.loads(FIXTURES.read_text())
    for i, smiles in enumerate(fixture["smiles"]):
        yield f"recipe-{i:02d}", smiles
    edges = list(itertools.combinations(range(5), 2))
    for mask in range(1024):
        model = Chem.RWMol()
        for _ in range(5):
            model.AddAtom(Chem.Atom(6))
        for bit, (a, b) in enumerate(edges):
            if mask & (1 << bit):
                model.AddBond(a, b, Chem.BondType.SINGLE)
        smiles = Chem.MolToSmiles(model, canonical=False)
        yield f"graph5-{mask:04d}", smiles
    for size in range(1, 33):
        yield f"chain-{size:02d}", "C" * size
        if size > 2:
            yield f"cycle-{size:02d}", "C1" + "C" * (size - 1) + "1"


def oracle(smiles):
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params)
    assert mol is not None, smiles
    model = Chem.RWMol(mol)
    # Explicit vertex deletion avoids RDKit's isotope/bridging-H removal policy.
    for atom in reversed(list(mol.GetAtoms())):
        if atom.GetAtomicNum() == 1:
            model.RemoveAtom(atom.GetIdx())
    matrix = Chem.GetDistanceMatrix(model, useBO=False, useAtomWts=False)
    return [[4294967295 if value >= 1e8 else int(value) for value in row] for row in matrix]


def main():
    assert rdBase.rdkitVersion == "2026.09.1", rdBase.rdkitVersion
    records = []
    cells = 0
    for identifier, smiles in packet():
        expected = oracle(smiles)
        actual = chematic.from_smiles(smiles).topological_distance_matrix()
        assert actual == expected, (identifier, smiles, actual, expected)
        cells += sum(map(len, actual))
        records.append(dict(id=identifier, smiles=smiles, payload_sha256=digest(smiles),
                            expected_sha256=digest(expected), rows=len(actual)))
    guards = []
    for n in [1025, 46341, 65535, 65536, 65537]:
        mol = chematic.from_smiles(".".join(["C"] * n))
        try:
            mol.topological_distance_matrix(max_atoms=sys.maxsize)
        except ValueError as error:
            assert f"{n} atoms exceeds 1024" in str(error), str(error)
            guards.append(dict(input_atoms=n, exception=type(error).__name__, message=str(error)))
        else:
            raise AssertionError(f"missing cap: {n}")
    result = dict(
        schema_version=1, issue=dict(number=9647, state="OPEN", verified_date="2026-10-15",
            url="https://github.com/rdkit/rdkit/issues/9647", updated_at="2026-10-10T16:48:18Z"),
        scope="Checked native heavy-subgraph matrix API and existing Python method; no descriptor-wide guard claim",
        oracle=dict(rdkit=rdBase.rdkitVersion, version_pin="2026.09.1",
            algorithm="Explicit H-vertex deletion, RDKit unweighted distance matrix, disconnected sentinel normalization"),
        fixtures=dict(path=str(FIXTURES.relative_to(ROOT)), sha256=hashlib.sha256(FIXTURES.read_bytes()).hexdigest(),
                      provenance=json.loads(FIXTURES.read_text())),
        results=dict(exact_matrix_cases=len(records), exact_cells=cells, mismatches=0,
                     oversized_input_guards=guards, upstream_oversized_crash_executed=False),
        records=records,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["results"], indent=2))


if __name__ == "__main__":
    main()
