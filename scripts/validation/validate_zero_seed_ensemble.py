#!/usr/bin/env python3
"""Bounded RDKit #9738 diagnostic: 4 flexible graphs, 3 seeds, 5 conformers.

Run once per engine/runtime; output records coordinates and rigid-motion-free
heavy-atom distance differences. No alignment, energy or coverage claim.
"""
import argparse
import datetime
import gzip
import hashlib
import itertools
import json
import math
import pathlib
import platform
import sys

SMILES = ["CCCCCCO", "CCCCCC", "NCCCCCC", "CC(C)CCCO"]


def digest(path):
    path = pathlib.Path(path)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=["chematic", "rdkit"], required=True)
    parser.add_argument("--python-root")
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if args.python_root:
        sys.path.insert(0, args.python_root)
    if args.engine == "chematic":
        import chematic
        from chematic import chematic as native
        runtime = digest(native.__file__)
        runtime["version"] = getattr(chematic, "__version__", "unknown")
    else:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdDistGeom
        runtime = digest(rdDistGeom.__file__)
        runtime["version"] = rdkit.__version__
    cases = []
    for smiles, seed in itertools.product(SMILES, [0, 1, 42]):
        if args.engine == "chematic":
            graph = chematic.from_smiles(smiles)
            heavy_atoms = graph.atom_count
            molecule = graph.add_hydrogens()
            run = lambda n: molecule.rdkit_embed_multiple_confs(n, random_seed=seed)
        else:
            graph = Chem.MolFromSmiles(smiles)
            heavy_atoms = graph.GetNumAtoms()
            molecule = Chem.AddHs(graph)

            def run(n):
                ids = AllChem.EmbedMultipleConfs(
                    molecule, numConfs=n, randomSeed=seed,
                    pruneRmsThresh=-1.0, numThreads=1,
                )
                return [molecule.GetConformer(i).GetPositions().tolist() for i in ids]

        coords = run(5)
        assert len(coords) == 5
        assert all(math.isfinite(x) for c in coords for p in c for x in p)
        repeat = run(5)
        prefix = run(3)
        pairs = list(itertools.combinations(range(heavy_atoms), 2))
        differences = [
            max(abs(math.dist(a[i], a[j]) - math.dist(b[i], b[j])) for i, j in pairs)
            for a, b in itertools.combinations(coords, 2)
        ]
        cases.append({
            "smiles": smiles, "seed": seed, "heavy_atoms": heavy_atoms,
            "conformers": coords, "repeat_equal": coords == repeat,
            "prefix_equal": coords[:3] == prefix,
            "max_heavy_distance_change_by_conformer_pair": differences,
            "has_nonrigid_variation_over_0_1_angstrom": max(differences) > 0.1,
            "unique_coordinate_arrays": len({json.dumps(c) for c in coords}),
        })
    result = {
        "schema": "chematic.zero-seed-ensemble.v1", "engine": args.engine,
        "label": args.label, "runtime": runtime, "observer": digest(__file__),
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "platform": platform.platform(), "python": platform.python_version(),
        "issue": "https://github.com/rdkit/rdkit/issues/9738",
        "parameters": {"num_confs": 5, "seeds": [0, 1, 42], "pruning": False,
                       "explicit_hydrogens": True, "mmff": False},
        "cases": cases,
    }
    payload = json.dumps(result, indent=2, allow_nan=False).encode() + b"\n"
    if args.output.suffix == ".gz":
        payload = gzip.compress(payload, mtime=0)
    args.output.write_bytes(payload)
    for seed in [0, 1, 42]:
        selected = [case for case in cases if case["seed"] == seed]
        print(json.dumps({"seed": seed, "cases": len(selected),
                          "nonrigid_variation": sum(c["has_nonrigid_variation_over_0_1_angstrom"] for c in selected),
                          "reproducible": sum(c["repeat_equal"] for c in selected),
                          "prefix_stable": sum(c["prefix_equal"] for c in selected)}))


if __name__ == "__main__":
    main()
