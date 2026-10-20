#!/usr/bin/env python3
"""Record checked matrix evidence and ordinary RDKit parity on synthetic inputs."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
from pathlib import Path
import random
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def classified(value):
    if math.isnan(value):
        return {"category": "nan"}
    if math.isinf(value):
        return {"category": "infinity", "sign": 1 if value > 0 else -1}
    return {"category": "finite", "value": value}


def rdkit_matrix(Chem, smiles, coords, weights):
    mol = Chem.MolFromSmiles(smiles)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(coords):
        conformer.SetAtomPosition(i, xyz)
    mol.AddConformer(conformer)
    return Chem.Get3DDistanceMatrix(mol, useAtomWts=weights).tolist()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "validation/results/rdkit-checked-3d-matrix-2026-10-16.json")
    args = parser.parse_args()
    import chematic
    from rdkit import Chem, rdBase
    rng = random.Random(9647)
    ordinary = []
    max_difference = 0.0
    legacy_disagreements = 0
    for smiles in ["CCO", "C=CO", "c1ccccc1", "C.O.N", "[H]C([H])([H])[H]", "[2H]C", "[NH4+]", "ClCBr"]:
        n = Chem.MolFromSmiles(smiles).GetNumAtoms()
        native = chematic.from_smiles(smiles)
        for _ in range(8):
            coords = [[rng.randint(-80,80)/4 for _ in range(3)] for _ in range(n)]
            for weights in [False,True]:
                matrix = native.rdkit_distance_matrix_3d_checked(coords,weights)
                legacy = native.rdkit_distance_matrix_3d(coords,weights)
                reference = rdkit_matrix(Chem,smiles,coords,weights)
                legacy_disagreements += int(matrix != legacy)
                difference = max(abs(a-b) for ra,rb in zip(matrix,reference) for a,b in zip(ra,rb))
                max_difference = max(max_difference,difference)
                assert difference < 1e-12
                ordinary.append({"smiles":smiles,"coords":coords,"weights":weights,"checked":matrix,"rdkit":reference})
    assert legacy_disagreements == 0
    fixture = ROOT / "validation/rdkit_issues/checked-3d-matrix/cases.json"
    edges = []
    for case in json.loads(fixture.read_text(encoding="utf-8")):
        coords = [[float(value) for value in row] for row in case["coords"]]
        mol = chematic.from_smiles(case["smiles"])
        legacy = mol.rdkit_distance_matrix_3d(coords,case["weights"])
        reference = rdkit_matrix(Chem,case["smiles"],coords,case["weights"])
        result = {"id":case["id"],"legacy":classified(legacy[0][-1]),"rdkit":classified(reference[0][-1])}
        try:
            checked = mol.rdkit_distance_matrix_3d_checked(coords,case["weights"])
            assert "expected_error" not in case
            assert all(math.isfinite(value) for row in checked for value in row)
            if "expected_distance" in case:
                assert checked[0][1] == case["expected_distance"]
            result["checked"] = classified(checked[0][-1])
        except ValueError as error:
            assert case["expected_error"] in str(error)
            result["checked"] = {"category":"rejected","reason":str(error)}
        edges.append(result)
    packet = args.output.with_suffix(".packet.json")
    packet.write_text(json.dumps({"ordinary":ordinary,"numeric_edges":edges},indent=2,allow_nan=False)+"\n",encoding="utf-8")
    sources = ["crates/chematic-smiles/src/rdkit/matrices/checked_3d.rs", "crates/chematic-smiles/src/rdkit/matrices.rs", "crates/chematic-py/src/mol_methods/transformations.rs", "crates/chematic-wasm/src/checked_matrices.rs"]
    report = {"schema_version":1,"checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "issue":{"url":"https://github.com/rdkit/rdkit/issues/9647","state_checked":"open","relation":"Numeric extension of the existing allocation guard; no new upstream heap-crash claim"},
        "source_base_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "source_candidate":True,"published_artifact":False,"source_sha256":{p:sha(ROOT/p) for p in sources},
        "extension_sha256":sha(importlib.import_module("chematic.chematic").__file__),
        "rdkit":{"version":rdBase.rdkitVersion,"rdBase_sha256":sha(rdBase.__file__)},
        "ordinary_cases":len(ordinary),"max_absolute_rdkit_difference":max_difference,"legacy_disagreements":legacy_disagreements,
        "numeric_edge_cases":len(edges),"failed_cases":0,"oversized_rdkit_crash_reproducer_executed":False,
        "fixture_sha256":sha(fixture),"packet_sha256":sha(packet)}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"ordinary_cases":len(ordinary),"max_rdkit_difference":max_difference,"legacy_disagreements":legacy_disagreements,"numeric_edge_cases":len(edges)}))


if __name__ == "__main__":
    main()
