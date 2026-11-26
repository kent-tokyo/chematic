#!/usr/bin/env python3
"""Run one-conformer ETKDG/MMFF attempts under a shared per-engine time budget."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import multiprocessing
import platform
import time
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from platinum_development_probe import graph_matches, rms_values, stereo_violations
from platinum_geometry import check_geometry, prepare_geometry_checks


def attempts(engine, mol, other, seconds, max_candidates, repair_iters=0, on_event=None,
             profile="corrected"):
    accepted = []
    statuses = []
    force_field_errors = []
    no_conformer = 0
    stereo_rejected = {"violated": 0, "unevaluable": 0, "invalid_coordinates": 0,
                       "unclassified_attempts": 0}
    attempted = 0
    start = time.perf_counter()
    deadline = start + seconds
    while attempted < max_candidates and time.perf_counter() < deadline:
        seed = 42 * (attempted + 1)
        attempted += 1
        if on_event is not None:
            on_event(("started", attempted))
        if engine == "chematic":
            if profile == "mixed-fallback":
                result, embedded, _, violated, unevaluable, invalid = (
                    other.etkdg_embed_multiple_confs_mixed_checked_diagnostics(
                        1, random_seed=seed, max_mmff_iters=200, repair_mmff_iters=repair_iters)
                )
                rejected = {"violated": violated or 0, "unevaluable": unevaluable or 0,
                            "invalid_coordinates": invalid or 0,
                            "unclassified_attempts": int(embedded > 0 and violated is None)}
                for key, value in rejected.items():
                    stereo_rejected[key] += value
                if on_event is not None:
                    on_event(("stereo_rejected", rejected))
            else:
                result = other.etkdg_embed_multiple_confs_mmff_corrected(
                    1, random_seed=seed, max_mmff_iters=200, repair_mmff_iters=repair_iters)
            coords = result[0][2] if result else None
            status = result[0][0] if result else 1
            force_field_error = result[0][3] if result and profile == "mixed-fallback" else None
        else:
            force_field_error = None
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
            force_field_errors.append(force_field_error)
            if on_event is not None:
                on_event(("accepted", attempted, coords, status, force_field_error))
        elif coords is None:
            no_conformer += 1
            if on_event is not None:
                on_event(("no_conformer", attempted))
    return accepted, {
        "attempted": attempted,
        "accepted": len(accepted),
        "elapsed": time.perf_counter() - start,
        "statuses": statuses,
        "force_field_errors": force_field_errors,
        "no_conformer": no_conformer,
        "stereo_rejected": stereo_rejected if engine == "chematic" and profile == "mixed-fallback" else None,
    }


def _hard_budget_worker(connection, engine, smiles, seconds, max_candidates, repair_iters, profile):
    """Prepare each engine before starting its timed work in a separate process."""
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        other = chematic.from_smiles(smiles).add_hydrogens()
        connection.send(("ready",))
        connection.recv()
        _, timing = attempts(engine, mol, other, seconds, max_candidates, repair_iters,
                             connection.send, profile)
        connection.send(("done", timing))
    except Exception as exc:
        connection.send(("error", repr(exc)))
    finally:
        connection.close()


def hard_attempts(engine, smiles, seconds, max_candidates, repair_iters=0, profile="corrected"):
    """Bound one engine by a wall clock deadline, retaining completed attempts."""
    context = multiprocessing.get_context("spawn")
    parent, child = context.Pipe()
    process = context.Process(target=_hard_budget_worker,
                              args=(child, engine, smiles, seconds, max_candidates, repair_iters, profile))
    process.start()
    child.close()
    accepted, statuses, force_field_errors = [], [], []
    attempted = 0
    no_conformer = 0
    stereo_rejected = {"violated": 0, "unevaluable": 0, "invalid_coordinates": 0,
                       "unclassified_attempts": 0}
    try:
        if not parent.poll(30):
            raise RuntimeError("engine worker did not initialize within 30 seconds")
        message = parent.recv()
        if message[0] != "ready":
            raise RuntimeError(f"engine worker initialization failed: {message}")
        start = time.perf_counter()
        deadline = start + seconds
        parent.send("start")
        finished = False
        while (remaining := deadline - time.perf_counter()) > 0 and parent.poll(remaining):
            message = parent.recv()
            if message[0] == "started":
                attempted = message[1]
            elif message[0] == "accepted":
                attempted = message[1]
                accepted.append(message[2])
                statuses.append(message[3])
                force_field_errors.append(message[4])
            elif message[0] == "no_conformer":
                no_conformer += 1
            elif message[0] == "stereo_rejected":
                for key, value in message[1].items():
                    stereo_rejected[key] += value
            elif message[0] == "done":
                attempted = message[1]["attempted"]
                finished = True
                break
            else:
                raise RuntimeError(f"engine worker failed: {message}")
        elapsed = time.perf_counter() - start
    finally:
        if process.is_alive():
            process.terminate()
        process.join()
        parent.close()
    return accepted, {"attempted": attempted, "accepted": len(accepted),
                      "elapsed": elapsed, "elapsed_with_cleanup": time.perf_counter() - start,
                      "statuses": statuses, "force_field_errors": force_field_errors,
                      "no_conformer": no_conformer,
                      "stereo_rejected": stereo_rejected if engine == "chematic" and profile == "mixed-fallback" else None,
                      "timed_out": not finished,
                      "deadline_reached": not finished or attempted < max_candidates,
                      "attempted_is_lower_bound": not finished}


def score_engines(result, order, mol, other, reference, checks, seconds, max_candidates, repair_iters,
                  hard_timeout_smiles=None, profile="corrected"):
    """Retain and score each engine independently if the other engine fails."""
    ensembles = {}
    engine_errors = {}
    for engine in order:
        try:
            if hard_timeout_smiles is None:
                conformers, timing = attempts(engine, mol, other, seconds, max_candidates, repair_iters,
                                              profile=profile) if profile != "corrected" else attempts(
                    engine, mol, other, seconds, max_candidates, repair_iters)
            else:
                conformers, timing = hard_attempts(
                    engine, hard_timeout_smiles, seconds, max_candidates, repair_iters, profile
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
    parser.add_argument("--hard-timeout", action="store_true",
                        help="run each engine in a process terminated at the deadline")
    parser.add_argument("--profile", choices=("corrected", "mixed-fallback"), default="corrected")
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
        "hard_timeout": args.hard_timeout,
        "profile": args.profile,
        "engine_error_policy": "independent generation and scoring per engine",
        "seed_schedule": "42 * (attempt + 1)",
        "timing": ("worker deadline after initialization; scoring and process startup excluded"
                   if args.hard_timeout else
                   "individual embedding and minimization attempts, accept only attempts completed before the deadline; scoring excluded"),
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
                hard_timeout_smiles=row["smiles"] if args.hard_timeout else None,
                profile=args.profile,
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
