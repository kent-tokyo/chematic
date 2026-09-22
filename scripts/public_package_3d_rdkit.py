#!/usr/bin/env python3
"""Run the published RDKit wheel through chematic's fixed 3D corpus."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RANDOM_SEED = 20260801
MAX_ATTEMPTS = 8
FORCE_FIELD_MAX_ITERATIONS = 200
BEST_OF_N = 10

ARMS = (
    "rdkit_etkdgv3_raw",
    "rdkit_etkdgv3_uff",
    "rdkit_etkdgv3_mmff94",
    "rdkit_etkdgv3_best_of_n",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_rows(tiers: str) -> list[dict[str, object]]:
    selected = set(tiers)
    rows: list[dict[str, object]] = []
    for tier, relative in (
        ("A", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
        ("B", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
    ):
        if tier not in selected:
            continue
        manifest = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        for molecule in manifest["molecules"]:
            rows.append({"tier": tier, **molecule})
    return rows


def params(AllChem: object) -> object:
    value = AllChem.ETKDGv3()  # type: ignore[attr-defined]
    value.useSmallRingTorsions = True
    value.useMacrocycleTorsions = True
    value.useMacrocycle14config = True
    value.enforceChirality = True
    return value


def embed_with_retry(
    AllChem: object, mol: object, parameters: object
) -> tuple[int, int]:
    for attempt in range(MAX_ATTEMPTS):
        parameters.randomSeed = RANDOM_SEED + attempt
        conformer_id = AllChem.EmbedMolecule(mol, parameters)  # type: ignore[attr-defined]
        if conformer_id >= 0:
            return conformer_id, attempt + 1
    return -1, MAX_ATTEMPTS


def coordinates(mol: object, conformer_id: int, heavy_atoms: int) -> list[list[float]]:
    conformer = mol.GetConformer(conformer_id)
    return [list(conformer.GetAtomPosition(index)) for index in range(heavy_atoms)]


def run_single(
    Chem: object, AllChem: object, mol_without_h: object, arm: str
) -> dict[str, object]:
    heavy_atoms = mol_without_h.GetNumAtoms()
    started = time.perf_counter_ns()
    mol = Chem.AddHs(mol_without_h)
    parameters = params(AllChem)
    conformer_id, attempts = embed_with_retry(AllChem, mol, parameters)
    if conformer_id < 0:
        return {
            "status": "typed_failure",
            "failure_cause": "EmbedMolecule_failed",
            "elapsed_ms": (time.perf_counter_ns() - started) / 1_000_000,
        }

    return_code = None
    if arm == "rdkit_etkdgv3_uff":
        return_code = AllChem.UFFOptimizeMolecule(
            mol, confId=conformer_id, maxIters=FORCE_FIELD_MAX_ITERATIONS
        )
    elif arm == "rdkit_etkdgv3_mmff94":
        properties = AllChem.MMFFGetMoleculeProperties(mol)
        if properties is None:
            return {
                "status": "typed_failure",
                "failure_cause": "MMFF_parameters_unavailable",
                "elapsed_ms": (time.perf_counter_ns() - started) / 1_000_000,
                "embed_attempts_used": attempts,
            }
        return_code = AllChem.MMFFOptimizeMolecule(
            mol, confId=conformer_id, maxIters=FORCE_FIELD_MAX_ITERATIONS
        )

    elapsed_ms = (time.perf_counter_ns() - started) / 1_000_000
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "embed_attempts_used": attempts,
        "coords": coordinates(mol, conformer_id, heavy_atoms),
        "force_field": {
            "rdkit_etkdgv3_raw": "none",
            "rdkit_etkdgv3_uff": "uff",
            "rdkit_etkdgv3_mmff94": "mmff94",
        }[arm],
        "force_field_return_code": return_code,
        "force_field_converged": return_code == 0 if return_code is not None else True,
    }


def run_best_of_n(
    Chem: object, AllChem: object, mol_without_h: object
) -> dict[str, object]:
    heavy_atoms = mol_without_h.GetNumAtoms()
    started = time.perf_counter_ns()
    mol = Chem.AddHs(mol_without_h)
    parameters = params(AllChem)
    parameters.randomSeed = RANDOM_SEED
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(mol, numConfs=BEST_OF_N, params=parameters)
    )
    if not conformer_ids:
        return {
            "status": "typed_failure",
            "failure_cause": "EmbedMultipleConfs_failed",
            "elapsed_ms": (time.perf_counter_ns() - started) / 1_000_000,
        }

    candidates: list[tuple[float, int]] = []
    for conformer_id in conformer_ids:
        try:
            AllChem.UFFOptimizeMolecule(
                mol,
                confId=conformer_id,
                maxIters=FORCE_FIELD_MAX_ITERATIONS,
            )
            force_field = AllChem.UFFGetMoleculeForceField(mol, confId=conformer_id)
            if force_field is not None:
                candidates.append((force_field.CalcEnergy(), conformer_id))
        except Exception:
            continue
    elapsed_ms = (time.perf_counter_ns() - started) / 1_000_000
    if not candidates:
        return {
            "status": "typed_failure",
            "failure_cause": "all_conformers_failed_uff_optimization",
            "elapsed_ms": elapsed_ms,
            "conformers_embedded": len(conformer_ids),
        }
    energy, conformer_id = min(candidates)
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "coords": coordinates(mol, conformer_id, heavy_atoms),
        "force_field": "uff",
        "conformers_embedded": len(conformer_ids),
        "conformers_optimized": len(candidates),
        "best_uff_energy": energy,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--tiers", choices=("A", "B", "AB"), default="AB")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int)
    parser.add_argument("--arms", nargs="+", choices=ARMS, default=list(ARMS))
    parser.add_argument("--wall-budget-seconds", type=float, default=3300.0)
    args = parser.parse_args()
    if args.start < 0 or (args.count is not None and args.count <= 0):
        parser.error("--start must be non-negative and --count must be positive")
    if args.wall_budget_seconds <= 0:
        parser.error("--wall-budget-seconds must be positive")
    if not args.wheel.is_file():
        parser.error(f"wheel not found: {args.wheel}")
    return args


def main() -> int:
    args = parse_args()
    import rdkit
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem

    all_rows = load_rows(args.tiers)
    selected = all_rows[args.start :]
    if args.count is not None:
        selected = selected[: args.count]
    session_started = time.monotonic()
    emitted = 0
    termination = "completed"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        for row_index, row in enumerate(selected, start=args.start):
            if time.monotonic() - session_started >= args.wall_budget_seconds:
                termination = "wall_budget_exhausted"
                break
            smiles = str(row["smiles"])
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                result_rows = [
                    {
                        "arm": arm,
                        "status": "parse_failure",
                        "failure_cause": "MolFromSmiles returned None",
                        "elapsed_ms": 0.0,
                    }
                    for arm in args.arms
                ]
            else:
                result_rows = []
                for arm in args.arms:
                    try:
                        result = (
                            run_best_of_n(Chem, AllChem, mol)
                            if arm == "rdkit_etkdgv3_best_of_n"
                            else run_single(Chem, AllChem, mol, arm)
                        )
                    except Exception as error:
                        result = {
                            "status": "internal_error",
                            "failure_cause": f"{type(error).__name__}: {error}",
                            "elapsed_ms": 0.0,
                        }
                    result_rows.append({"arm": arm, **result})
            for result in result_rows:
                output_row = {
                    "tier": row["tier"],
                    "row_index": row_index,
                    "name": row["name"],
                    "smiles": smiles,
                    "primary_category": row.get("primary_category", "unknown"),
                    "engine": "rdkit_public_python",
                    **result,
                }
                stream.write(json.dumps(output_row, separators=(",", ":")) + "\n")
                stream.flush()
                emitted += 1

    metadata = {
        "schema_version": 1,
        "benchmark": "public-package-3d-vs-rdkit",
        "engine": "rdkit_public_python",
        "package": {
            "name": "rdkit",
            "distribution_version": importlib.metadata.version("rdkit"),
            "runtime_version": rdBase.rdkitVersion,
            "module_path": rdkit.__file__,
            "wheel": {
                "path": str(args.wheel),
                "bytes": args.wheel.stat().st_size,
                "sha256": sha256(args.wheel),
            },
        },
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "python_executable": sys.executable,
        },
        "configuration": {
            "random_seed": RANDOM_SEED,
            "max_attempts": MAX_ATTEMPTS,
            "force_field_max_iterations": FORCE_FIELD_MAX_ITERATIONS,
            "best_of_n": BEST_OF_N,
            "add_hs_inside_timed_operation": True,
            "use_small_ring_torsions": True,
            "use_macrocycle_torsions": True,
            "use_macrocycle_14_bounds": True,
            "enforce_chirality": True,
            "tiers": args.tiers,
            "start": args.start,
            "requested_count": args.count,
            "arms": args.arms,
            "wall_budget_seconds": args.wall_budget_seconds,
        },
        "result": {
            "termination": termination,
            "selected_molecules": len(selected),
            "emitted_rows": emitted,
            "output": str(args.output),
            "output_sha256": sha256(args.output),
            "elapsed_seconds": time.monotonic() - session_started,
        },
        "collected_at_utc": datetime.now(timezone.utc).isoformat(),
        "boundary": "Published wheel only. Parsing is outside per-arm timing. AddHs is inside RDKit operation timing. Coordinates and failures are retained for the identical external scorer.",
    }
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata["result"], indent=2))
    return 0 if termination == "completed" else 3


if __name__ == "__main__":
    raise SystemExit(main())
