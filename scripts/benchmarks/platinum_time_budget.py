#!/usr/bin/env python3
"""Run one-conformer ETKDG/MMFF attempts under a shared per-engine time budget."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import time
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from platinum_development_probe import graph_matches, rms_values, stereo_violations
from platinum_geometry import check_geometry, prepare_geometry_checks


def attempts(engine, mol, other, seconds, max_candidates, repair_iters=0):
    accepted = []
    statuses = []
    attempted = 0
    start = time.perf_counter()
    deadline = start + seconds
    while attempted < max_candidates and time.perf_counter() < deadline:
        seed = 42 * (attempted + 1)
        attempted += 1
        if engine == "chematic":
            result = other.etkdg_embed_multiple_confs_mmff_corrected(
                1, random_seed=seed, max_mmff_iters=200, repair_mmff_iters=repair_iters
            )
            coords = result[0][2] if result else None
            status = result[0][0] if result else 1
        else:
            generated = Chem.Mol(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.numThreads = 1
            if AllChem.EmbedMolecule(generated, params) == 0:
                status = AllChem.MMFFOptimizeMolecule(generated, maxIters=200)
                if status == 0:
                    props = AllChem.MMFFGetMoleculeProperties(generated)
                    field = AllChem.MMFFGetMoleculeForceField(generated, props)
                    gradient = field.CalcGrad()
                    if any(not math.isfinite(g) or abs(g) > 0.01 for g in gradient):
                        status = 1
                coords = [
                    list(generated.GetConformer().GetAtomPosition(i))
                    for i in range(generated.GetNumAtoms())
                ]
            else:
                coords = None
        if coords is not None and time.perf_counter() <= deadline:
            accepted.append(coords)
            statuses.append(status)
    return accepted, {
        "attempted": attempted,
        "accepted": len(accepted),
        "elapsed": time.perf_counter() - start,
        "statuses": statuses,
    }


def score_engines(result, order, mol, other, reference, checks, seconds, max_candidates, repair_iters):
    """Retain and score each engine independently if the other engine fails."""
    ensembles = {}
    engine_errors = {}
    for engine in order:
        try:
            conformers, timing = attempts(
                engine, mol, other, seconds, max_candidates, repair_iters
            )
            ensembles[engine] = conformers
            result[engine] = timing
        except Exception as exc:
            engine_errors[engine] = {"stage": "generation", "error": repr(exc)}
    for engine, conformers in ensembles.items():
        try:
            result[engine]["rmsd"] = rms_values(mol, reference, conformers)
            result[engine]["stereo"] = stereo_violations(mol, conformers)
            if checks is not None:
                result[engine]["geometry"] = [check_geometry(coords, checks) for coords in conformers]
        except Exception as exc:
            engine_errors[engine] = {"stage": "scoring", "error": repr(exc)}
    if engine_errors:
        result["engine_errors"] = engine_errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--seconds", type=float, default=1.0)
    parser.add_argument("--max-candidates", type=int, default=20)
    parser.add_argument("--repair-iters", type=int, default=0)
    parser.add_argument("--geometry", action="store_true")
    args = parser.parse_args()
    if (
        args.limit < 1
        or not math.isfinite(args.seconds)
        or args.seconds <= 0
        or args.max_candidates < 1
        or args.repair_iters < 0
    ):
        parser.error("limits must be positive")
    baseline = json.loads(args.baseline.read_text())
    if hashlib.sha256(args.sdf.read_bytes()).hexdigest() != baseline["source_sha256"]:
        parser.error("SDF hash differs from baseline")
    rows = sorted(
        baseline["results"], key=lambda r: hashlib.sha256(r["name"].encode()).digest()
    )[: args.limit]
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    output = {
        "rdkit": rdBase.rdkitVersion,
        "chematic": importlib.metadata.version("chematic"),
        "platform": platform.platform(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": baseline["source_sha256"],
        "seconds": args.seconds,
        "max_candidates": args.max_candidates,
        "repair_iters": args.repair_iters,
        "geometry_checked": args.geometry,
        "engine_error_policy": "independent generation and scoring per engine",
        "seed_schedule": "42 * (attempt + 1)",
        "timing": "individual embedding and minimization attempts, accept only attempts completed before the deadline; scoring excluded",
        "results": [],
    }
    warm_mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    warm_other = chematic.from_smiles("CCO").add_hydrogens()
    for engine in ("chematic", "rdkit"):
        _, warmup = attempts(engine, warm_mol, warm_other, 5.0, 1)
        if warmup["accepted"] != 1:
            raise RuntimeError(f"{engine} warmup did not complete")
    for index, row in enumerate(rows):
        result = {"name": row["name"], "index": row["index"]}
        try:
            reference = supplier[row["index"]]
            assert reference is not None and reference.GetProp("_Name") == row["name"]
            reference = Chem.RemoveHs(reference)
            mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
            other = chematic.from_smiles(row["smiles"]).add_hydrogens()
            assert graph_matches(other, mol)
            checks = prepare_geometry_checks(mol) if args.geometry else None
            if checks is not None:
                result["geometry_reference"] = {
                    "bonds_supported": len(checks[0]),
                    "bonds_total": mol.GetNumBonds(),
                    "angles_supported": len(checks[1]),
                    "angles_total": sum(
                        a.GetDegree() * (a.GetDegree() - 1) // 2 for a in mol.GetAtoms()
                    ),
                }
            order = ("chematic", "rdkit") if index % 2 == 0 else ("rdkit", "chematic")
            result["order"] = order
            score_engines(
                result, order, mol, other, reference, checks,
                args.seconds, args.max_candidates, args.repair_iters,
            )
        except Exception as exc:
            result["error"] = repr(exc)
        output["results"].append(result)
        checkpoint = args.output.with_suffix(args.output.suffix + ".tmp")
        checkpoint.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        checkpoint.replace(args.output)
        print(row["name"], result.get("error") or result.get("engine_errors") or "ok", flush=True)


if __name__ == "__main__":
    main()
