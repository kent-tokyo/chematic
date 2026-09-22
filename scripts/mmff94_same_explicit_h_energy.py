#!/usr/bin/env python3
"""Compare MMFF94 energies on one shared explicit-H coordinate set.

RDKit generates each conformer once, then RDKit and CheMatic evaluate the same
atom-ordered coordinates. This isolates force-field energy agreement from
conformer generation, minimization, stereo preservation, and runtime speed.
Every input row is retained with one terminal status.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

import chematic


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFESTS = (
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
)
TERMINAL_STATUSES = {
    "ok",
    "declared_unsupported",
    "parse_failure",
    "atom_order_mismatch",
    "embed_failure",
    "rdkit_unsupported",
    "rdkit_force_field_failure",
    "schematic_unsupported",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command: list[str]) -> str | None:
    try:
        return subprocess.run(
            command,
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def artifact_record(path: Path | None) -> dict | None:
    if path is None:
        return None
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValueError(f"artifact not found: {resolved}")
    return {
        "path": str(resolved),
        "bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def explicit_symbols(mol: Chem.Mol) -> list[str]:
    return [atom.GetSymbol() for atom in mol.GetAtoms()]


def terminal(row: dict, input_index: int, status: str, **details: object) -> dict:
    if status not in TERMINAL_STATUSES:
        raise ValueError(f"unknown terminal status: {status}")
    return {"input_index": input_index, **row, "status": status, **details}


def gradient_diagnostic(mol: object, coords: list[list[float]], delta: float) -> dict:
    analytic = mol.mmff94_bounded_analytic_gradient(coords)
    absolute_errors: list[float] = []
    scaled_errors: list[float] = []
    for atom_index, point in enumerate(coords):
        for axis in range(3):
            plus = [candidate.copy() for candidate in coords]
            minus = [candidate.copy() for candidate in coords]
            plus[atom_index][axis] += delta
            minus[atom_index][axis] -= delta
            finite_difference = (
                mol.mmff94_energy_breakdown(plus)["total"]
                - mol.mmff94_energy_breakdown(minus)["total"]
            ) / (2.0 * delta)
            error = abs(float(analytic[atom_index][axis]) - finite_difference)
            absolute_errors.append(error)
            scaled_errors.append(error / (1.0 + abs(finite_difference)))
    return {
        "components": len(absolute_errors),
        "central_difference_delta_angstrom": delta,
        "max_abs_error_kcal_mol_angstrom": max(absolute_errors, default=0.0),
        "max_scaled_error": max(scaled_errors, default=0.0),
        "rms_abs_error_kcal_mol_angstrom": (
            (sum(error * error for error in absolute_errors) / len(absolute_errors))
            ** 0.5
            if absolute_errors
            else 0.0
        ),
    }


def evaluate(
    row: dict,
    input_index: int,
    seed: int,
    gradient_indices: set[int],
    gradient_delta: float,
) -> dict:
    smiles = row["smiles"]
    if row.get("primary_category") == "force_field_unsupported":
        return terminal(row, input_index, "declared_unsupported")
    rdkit_base = Chem.MolFromSmiles(smiles)
    if rdkit_base is None:
        return terminal(row, input_index, "parse_failure", engine="rdkit")
    try:
        schematic_mol = chematic.from_smiles(smiles).add_hydrogens()
    except (RuntimeError, ValueError) as exc:
        return terminal(
            row,
            input_index,
            "parse_failure",
            engine="chematic",
            error={"type": type(exc).__name__, "message": str(exc)},
        )
    rdkit_mol = Chem.AddHs(rdkit_base)
    rdkit_symbols = explicit_symbols(rdkit_mol)
    schematic_symbols = [
        atom["element"] for atom in schematic_mol.depict_data()["atoms"]
    ]
    if rdkit_symbols != schematic_symbols:
        return terminal(
            row,
            input_index,
            "atom_order_mismatch",
            rdkit_atom_symbols=rdkit_symbols,
            schematic_atom_symbols=schematic_symbols,
        )
    if AllChem.EmbedMolecule(rdkit_mol, randomSeed=seed) != 0:
        return terminal(row, input_index, "embed_failure")
    properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
    if properties is None:
        return terminal(row, input_index, "rdkit_unsupported")
    force_field = AllChem.MMFFGetMoleculeForceField(rdkit_mol, properties, confId=0)
    if force_field is None:
        return terminal(row, input_index, "rdkit_force_field_failure")
    conformer = rdkit_mol.GetConformer()
    coords = [
        [position.x, position.y, position.z]
        for position in (
            conformer.GetAtomPosition(i) for i in range(rdkit_mol.GetNumAtoms())
        )
    ]
    try:
        breakdown = {
            key: float(value)
            for key, value in schematic_mol.mmff94_energy_breakdown(coords).items()
        }
    except (RuntimeError, ValueError) as exc:
        return terminal(
            row,
            input_index,
            "schematic_unsupported",
            error={"type": type(exc).__name__, "message": str(exc)},
        )
    schematic_energy = breakdown["total"]
    rdkit_energy = float(force_field.CalcEnergy())
    delta = schematic_energy - rdkit_energy
    coordinate_payload = json.dumps(coords, separators=(",", ":")).encode("utf-8")
    gradient = (
        gradient_diagnostic(schematic_mol, coords, gradient_delta)
        if input_index in gradient_indices
        else None
    )
    return terminal(
        row,
        input_index,
        "ok",
        atom_count=len(coords),
        coordinate_sha256=sha256_bytes(coordinate_payload),
        schematic_energy_breakdown_kcal_mol=breakdown,
        schematic_energy_kcal_mol=schematic_energy,
        rdkit_energy_kcal_mol=rdkit_energy,
        delta_kcal_mol=delta,
        abs_delta_kcal_mol=abs(delta),
        **({"gradient_diagnostic": gradient} if gradient is not None else {}),
    )


def percentile_nearest_rank(
    values: list[float], numerator: int, denominator: int
) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    rank = (len(ordered) * numerator + denominator - 1) // denominator
    return ordered[max(0, rank - 1)]


def summarize(results: list[dict]) -> dict:
    ok = [result for result in results if result["status"] == "ok"]
    deltas = [result["abs_delta_kcal_mol"] for result in ok]
    by_category: dict[str, list[float]] = defaultdict(list)
    for result in ok:
        by_category[result.get("primary_category", "unknown")].append(
            result["abs_delta_kcal_mol"]
        )
    gradient_rows = [result for result in ok if "gradient_diagnostic" in result]
    return {
        "row_accounting": {
            "input_count": len(results),
            "terminal_count": sum(
                result.get("status") in TERMINAL_STATUSES for result in results
            ),
            "status_counts": dict(Counter(result["status"] for result in results)),
        },
        "comparable_rows": len(ok),
        "median_abs_delta_kcal_mol": statistics.median(deltas) if deltas else None,
        "p90_abs_delta_kcal_mol": percentile_nearest_rank(deltas, 9, 10),
        "max_abs_delta_kcal_mol": max(deltas) if deltas else None,
        "within_1_kcal_mol": sum(delta <= 1.0 for delta in deltas),
        "within_5_kcal_mol": sum(delta <= 5.0 for delta in deltas),
        "gradient_diagnostic": {
            "rows": len(gradient_rows),
            "input_indices": [row["input_index"] for row in gradient_rows],
            "max_abs_error_kcal_mol_angstrom": max(
                (
                    row["gradient_diagnostic"]["max_abs_error_kcal_mol_angstrom"]
                    for row in gradient_rows
                ),
                default=None,
            ),
            "max_scaled_error": max(
                (
                    row["gradient_diagnostic"]["max_scaled_error"]
                    for row in gradient_rows
                ),
                default=None,
            ),
        },
        "by_primary_category": {
            category: {
                "rows": len(values),
                "median_abs_delta_kcal_mol": statistics.median(values),
                "max_abs_delta_kcal_mol": max(values),
            }
            for category, values in sorted(by_category.items())
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, action="append", default=[])
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--rdkit-artifact", type=Path)
    parser.add_argument("--schematic-artifact", type=Path)
    parser.add_argument("--schematic-build-command")
    parser.add_argument("--expected-rdkit")
    parser.add_argument("--expected-schematic")
    parser.add_argument("--gradient-input-index", type=int, action="append", default=[])
    parser.add_argument("--gradient-delta", type=float, default=1e-5)
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if any(index < 0 for index in args.gradient_input_index):
        parser.error("--gradient-input-index must be non-negative")
    if args.gradient_delta <= 0:
        parser.error("--gradient-delta must be positive")
    return args


def main() -> int:
    args = parse_args()
    manifests = args.manifest or list(DEFAULT_MANIFESTS)
    rows: list[dict] = []
    manifest_records: list[dict] = []
    for manifest in manifests:
        document = json.loads(manifest.read_text(encoding="utf-8"))
        molecules = document["molecules"]
        rows.extend(molecules)
        manifest_records.append(
            {
                "path": str(manifest.resolve()),
                "sha256": sha256_file(manifest),
                "rows": len(molecules),
            }
        )
    if args.limit is not None:
        rows = rows[: args.limit]
    gradient_indices = set(args.gradient_input_index)
    missing_gradient_indices = gradient_indices.difference(range(len(rows)))
    if missing_gradient_indices:
        raise ValueError(
            f"gradient input indices outside selected rows: {sorted(missing_gradient_indices)}"
        )

    results: list[dict] = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for index, row in enumerate(rows):
            result = evaluate(
                row,
                index,
                args.seed + index,
                gradient_indices,
                args.gradient_delta,
            )
            results.append(result)
            handle.write(json.dumps(result, sort_keys=True) + "\n")

    if args.summary:
        missing_dimensions = []
        if args.rdkit_artifact is None:
            missing_dimensions.append("rdkit_artifact")
        if args.schematic_artifact is None:
            missing_dimensions.append("schematic_artifact")
        if not args.schematic_build_command:
            missing_dimensions.append("schematic_build_command")
        git_status = command_output(
            ["git", "status", "--porcelain", "--untracked-files=no"]
        )
        measurements = summarize(results)
        gate = {
            "row_accounting_complete": measurements["row_accounting"]["terminal_count"]
            == len(results),
            "input_indices_contiguous": [row["input_index"] for row in results]
            == list(range(len(results))),
            "rdkit_version_matches": args.expected_rdkit is None
            or rdBase.rdkitVersion == args.expected_rdkit,
            "schematic_version_matches": args.expected_schematic is None
            or importlib.metadata.version("chematic") == args.expected_schematic,
            "provenance_complete": not missing_dimensions,
            "source_tree_clean": git_status == "",
        }
        summary = {
            "schema_version": 2,
            "profile": "mmff94_same_explicit_h_current_source_v2",
            "status": "current_source_diagnostic_not_release_gate",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "protocol": {
                "comparison": "RDKit and CheMatic MMFF94 total energy on the same RDKit-generated explicit-H coordinates",
                "coordinate_producer": "RDKit EmbedMolecule default parameters with per-row fixed random seed",
                "excludes": [
                    "conformer_quality",
                    "minimization_convergence",
                    "stereo_preservation",
                    "runtime_speed",
                    "RDKit_per_term_energy_parity",
                ],
                "seed": args.seed,
                "gradient_input_indices": sorted(gradient_indices),
                "gradient_central_difference_delta_angstrom": args.gradient_delta,
            },
            "versions": {
                "chematic": importlib.metadata.version("chematic"),
                "rdkit_distribution": importlib.metadata.version("rdkit"),
                "rdkit_runtime": rdBase.rdkitVersion,
                "python": platform.python_version(),
            },
            "source": {
                "git_revision": command_output(["git", "rev-parse", "HEAD"]),
                "tracked_tree_dirty": git_status != "",
                "rustc": command_output(["rustc", "-Vv"]),
                "maturin": command_output(["maturin", "--version"]),
                "build_command": args.schematic_build_command,
            },
            "host": {
                "platform": platform.platform(),
                "machine": platform.machine(),
                "processor": platform.processor(),
            },
            "artifacts": {
                "rdkit": artifact_record(args.rdkit_artifact),
                "schematic": artifact_record(args.schematic_artifact),
            },
            "manifests": manifest_records,
            "command": [sys.executable, *sys.argv],
            "rows_artifact": {
                "path": str(args.output.resolve()),
                "bytes": args.output.stat().st_size,
                "sha256": sha256_file(args.output),
            },
            "missing_dimensions": missing_dimensions,
            "measurements": measurements,
            "gate": gate,
            "gate_passed": all(gate.values()),
            "caveat": "RDKit generated the shared explicit-H coordinates. This measures total force-field energy agreement, not per-term RDKit parity, conformer quality, minimization convergence, stereo preservation, or speed.",
        }
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0 if summary["gate_passed"] else 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
