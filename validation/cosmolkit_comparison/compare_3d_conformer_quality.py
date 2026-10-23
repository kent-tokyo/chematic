#!/usr/bin/env python3
"""Compare converged MMFF94 minima from fixed ETKDGv3 ensembles.

The CheMatic source arm is produced by ``cargo run -p chematic-3d --example
cosmolkit_quality_probe``. This script runs the installed COSMolKit arm,
verifies that both arms cover the same manifest rows, and saves every row.
Lower MMFF94 energy is only a force-field search metric, not experimental
geometry accuracy. The quality profile may spend four times the 24-conformer
candidate budget, so only the compatibility profile is a same-budget lane.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_results(path: Path, names: set[str]) -> dict[tuple[str, str], dict]:
    result = {}
    for line in path.read_text().splitlines():
        row = json.loads(line)
        key = (row["name"], row["profile"])
        if key in result or row["name"] not in names or row["profile"] not in ("compat", "quality"):
            raise ValueError(f"unexpected or duplicate source row: {key}")
        if row["requested"] != 24 or not 0 <= row["converged"] <= row["retained"] <= 24:
            raise ValueError(f"invalid source count: {key}")
        energy = row["best_converged_energy"]
        if energy is not None and not math.isfinite(energy):
            raise ValueError(f"nonfinite source energy: {key}")
        result[key] = row
    expected = {(name, profile) for name in names for profile in ("compat", "quality")}
    if result.keys() != expected:
        raise ValueError(f"source coverage differs: {sorted(expected - result.keys())}")
    return result


def cosmolkit_minimum(smiles: str) -> dict:
    import cosmolkit as ck
    from rdkit import Chem

    reference = Chem.MolFromSmiles(smiles)
    candidate = ck.Molecule.from_smiles(smiles)
    candidate_readback = Chem.MolFromSmiles(candidate.to_smiles())
    if reference is None or candidate_readback is None or (
        Chem.MolToSmiles(reference) != Chem.MolToSmiles(candidate_readback)
    ):
        raise ValueError("COSMolKit graph differs on RDKit canonical read-back")

    mol = candidate.with_hydrogens()
    if mol.num_atoms() != Chem.AddHs(reference).GetNumAtoms():
        raise ValueError("explicit-hydrogen atom count differs from RDKit")
    params = ck.EmbedParams.etkdg_v3()
    params.random_seed = 42
    params.num_threads = 1
    embedded = mol.with_3d_conformers_result(24, params).molecule()
    optimized = embedded.with_mmff_optimized_conformers(
        num_threads=1, max_iterations=1000, non_bonded_threshold=100.0
    )
    outcomes = optimized.conformer_results()
    converged = [entry.energy() for entry in outcomes if entry.status_code() == 0]
    if any(not math.isfinite(value) for value in converged):
        raise ValueError("COSMolKit returned a nonfinite converged energy")
    return {
        "retained": len(outcomes),
        "converged": len(converged),
        "best_converged_energy": min(converged, default=None),
    }


def classify(cos_energy: float | None, source_energy: float | None) -> str:
    if cos_energy is None or source_energy is None:
        return "incomparable"
    difference = cos_energy - source_energy
    if difference > 0.01:
        return "chematic_lower"
    if difference < -0.01:
        return "cosmolkit_lower"
    return "within_0.01"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--source-results", type=Path, required=True)
    parser.add_argument("--source-bin", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=24)
    parser.add_argument("--indices", help="comma-separated manifest indices instead of start/limit")
    parser.add_argument("--cosmolkit-version", default="0.5.0")
    parser.add_argument("--rdkit-version", default="2026.09.1")
    args = parser.parse_args()

    import cosmolkit as ck
    from rdkit import rdBase

    if ck.__version__ != args.cosmolkit_version or rdBase.rdkitVersion != args.rdkit_version:
        raise ValueError("comparator version differs from the pinned benchmark")

    manifest = json.loads(args.manifest.read_text())
    if args.indices:
        indices = [int(value) for value in args.indices.split(",")]
        if len(indices) != len(set(indices)) or any(
            index < 0 or index >= len(manifest["molecules"]) for index in indices
        ):
            raise ValueError("duplicate or out-of-range manifest indices")
    else:
        indices = list(range(args.start, min(args.start + args.limit, len(manifest["molecules"]))))
    molecules = [manifest["molecules"][index] for index in indices]
    names = [row["name"] for row in molecules]
    if len(names) != len(set(names)):
        raise ValueError("duplicate names in manifest")
    source = source_results(args.source_results, set(names))
    rows = []
    for molecule in molecules:
        name = molecule["name"]
        cos = cosmolkit_minimum(molecule["smiles"])
        arms = {profile: source[(name, profile)] for profile in ("compat", "quality")}
        rows.append({
            "name": name,
            "smiles": molecule["smiles"],
            "cosmolkit": cos,
            "chematic": arms,
            "classification": {
                profile: classify(cos["best_converged_energy"], arm["best_converged_energy"])
                for profile, arm in arms.items()
            },
        })
    summary = {
        profile: {label: sum(row["classification"][profile] == label for row in rows)
                  for label in ("chematic_lower", "cosmolkit_lower", "within_0.01", "incomparable")}
        for profile in ("compat", "quality")
    }
    extension = next(Path(ck.__file__).parent.glob("cosmolkit*.so"))
    report = {
        "schema_version": 1,
        "metric": "minimum converged MMFF94 energy in kcal/mol; difference threshold 0.01",
        "conditions": {
            "embedding": "ETKDGv3, seed 42, one thread",
            "conformers_per_batch": 24,
            "max_mmff_iterations": 1000,
            "nonbonded_threshold": 100.0,
            "quality_budget": "2-4 batches of 24; conditional fourth batch",
        },
        "manifest": str(args.manifest),
        "manifest_indices": indices,
        "manifest_sha256": sha256(args.manifest),
        "source_results_sha256": sha256(args.source_results),
        "source_binary_sha256": sha256(args.source_bin),
        "cosmolkit_version": ck.__version__,
        "rdkit_version": rdBase.rdkitVersion,
        "cosmolkit_extension_sha256": sha256(extension),
        "python_version": platform.python_version(),
        "host": platform.platform(),
        "summary": summary,
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
