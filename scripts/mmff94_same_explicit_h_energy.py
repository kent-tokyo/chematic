#!/usr/bin/env python3
"""Compare MMFF94 on one shared explicit-H coordinate set.

RDKit generates the conformer once, then both RDKit and schematic evaluate
that same atom-ordered coordinate array.  This intentionally excludes
conformer-generation and minimisation quality from the numerical comparison.
Rows with different explicit-H atom ordering or unsupported MMFF typing are
retained with an explicit status instead of being silently dropped.
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

import chematic


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFESTS = (
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
)


def explicit_symbols(mol: Chem.Mol) -> list[str]:
    return [atom.GetSymbol() for atom in mol.GetAtoms()]


def evaluate(row: dict, seed: int) -> dict:
    smiles = row["smiles"]
    if row.get("primary_category") == "force_field_unsupported":
        return {**row, "status": "declared_unsupported"}
    rdkit_base = Chem.MolFromSmiles(smiles)
    if rdkit_base is None:
        return {**row, "status": "parse_failure"}
    rdkit_mol = Chem.AddHs(rdkit_base)
    schematic_mol = chematic.from_smiles(smiles).add_hydrogens()
    rdkit_symbols = explicit_symbols(rdkit_mol)
    schematic_symbols = [atom["element"] for atom in schematic_mol.depict_data()["atoms"]]
    if rdkit_symbols != schematic_symbols:
        return {
            **row,
            "status": "atom_order_mismatch",
            "rdkit_atom_symbols": rdkit_symbols,
            "schematic_atom_symbols": schematic_symbols,
        }
    if AllChem.EmbedMolecule(rdkit_mol, randomSeed=seed) != 0:
        return {**row, "status": "embed_failure"}
    properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
    if properties is None:
        return {**row, "status": "rdkit_unsupported"}
    force_field = AllChem.MMFFGetMoleculeForceField(
        rdkit_mol, properties, confId=0
    )
    if force_field is None:
        return {**row, "status": "rdkit_force_field_failure"}
    conformer = rdkit_mol.GetConformer()
    coords = [
        [position.x, position.y, position.z]
        for position in (conformer.GetAtomPosition(i) for i in range(rdkit_mol.GetNumAtoms()))
    ]
    try:
        schematic_energy = float(
            schematic_mol.mmff94_energy_breakdown(coords)["total"]
        )
    except (RuntimeError, ValueError):
        return {**row, "status": "schematic_unsupported"}
    rdkit_energy = float(force_field.CalcEnergy())
    return {
        **row,
        "status": "ok",
        "atom_count": len(coords),
        "schematic_energy_kcal_mol": schematic_energy,
        "rdkit_energy_kcal_mol": rdkit_energy,
        "delta_kcal_mol": schematic_energy - rdkit_energy,
        "abs_delta_kcal_mol": abs(schematic_energy - rdkit_energy),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, action="append", default=[])
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    manifests = args.manifest or list(DEFAULT_MANIFESTS)
    rows: list[dict] = []
    for manifest in manifests:
        rows.extend(json.loads(manifest.read_text(encoding="utf-8"))["molecules"])
    if args.limit is not None:
        rows = rows[: args.limit]
    results: list[dict] = []
    with args.output.open("w", encoding="utf-8") as handle:
        for index, row in enumerate(rows):
            result = evaluate(row, args.seed + index)
            results.append(result)
            handle.write(json.dumps(result, sort_keys=True) + "\n")
    if args.summary:
        ok = [result for result in results if result["status"] == "ok"]
        deltas = sorted(result["abs_delta_kcal_mol"] for result in ok)
        by_category: dict[str, list[float]] = defaultdict(list)
        for result in ok:
            by_category[result.get("primary_category", "unknown")].append(
                result["abs_delta_kcal_mol"]
            )
        args.summary.write_text(
            json.dumps(
                {
                    "status": "candidate_diagnostic_not_release_gate",
                    "protocol": "mmff94-same-explicit-h-coordinates",
                    "comparator": "RDKit MMFF94",
                    "seed": args.seed,
                    "rows": len(results),
                    "status_counts": dict(Counter(result["status"] for result in results)),
                    "comparable_rows": len(ok),
                    "median_abs_delta_kcal_mol": statistics.median(deltas) if deltas else None,
                    "p90_abs_delta_kcal_mol": deltas[max(0, (len(deltas) * 9 + 9) // 10 - 1)] if deltas else None,
                    "max_abs_delta_kcal_mol": max(deltas) if deltas else None,
                    "within_1_kcal_mol": sum(delta <= 1.0 for delta in deltas),
                    "within_5_kcal_mol": sum(delta <= 5.0 for delta in deltas),
                    "by_primary_category": {
                        category: {
                            "rows": len(values),
                            "median_abs_delta_kcal_mol": statistics.median(values),
                            "max_abs_delta_kcal_mol": max(values),
                        }
                        for category, values in sorted(by_category.items())
                    },
                    "caveat": "RDKit generated the shared explicit-H coordinates; this measures force-field energy parity, not conformer quality or minimization convergence.",
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
