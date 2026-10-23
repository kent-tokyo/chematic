#!/usr/bin/env python3
"""Score the next 21 eligible Tier B molecules after the original 12-row slice.

This is development evidence: the manifest was already used, and row 0022's
energy gap was observed before this durable record was created.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

import etkdg_quality_seed_budget_probe as prior


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def declared_stereo(smiles: str) -> tuple[dict[int, str], dict[int, str]]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    cip = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in mol.GetAtoms()
        if atom.HasProp("_CIPCode")
    }
    bonds = {
        bond.GetIdx(): str(bond.GetStereo())
        for bond in mol.GetBonds()
        if bond.GetStereo() != Chem.BondStereo.STEREONONE
    }
    return cip, bonds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    rows = prior.select_rows(start=12, count=21)
    assert len(rows) == 21
    binary = args.binary.resolve()
    manifest = prior.MANIFEST
    results = []
    for row in rows:
        run = subprocess.run(
            [str(binary)], input=json.dumps(row) + "\n", text=True,
            capture_output=True, check=True,
        )
        source = json.loads(run.stdout)
        assert source["name"] == row["name"]
        assert source["high_order_gradient_terms"] is False
        mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        prior.check_graph(source, mol)
        expected_cip, expected_bonds = declared_stereo(row["smiles"])
        candidate_errors = []
        cip_violations = 0
        bond_violations = 0
        for candidate in source["quality_candidates"]:
            energy, cip, bonds = prior.score_coords(mol, candidate["coords"])
            candidate_errors.append(abs(energy - candidate["energy"]))
            cip_violations += any(cip.get(i) != code for i, code in expected_cip.items())
            bond_violations += any(
                bonds.get(i) != stereo for i, stereo in expected_bonds.items()
            )
        assert len(candidate_errors) == source["quality"]["generated"]
        assert candidate_errors and max(candidate_errors) < 1e-6
        best_energy, _, _ = prior.score_coords(
            mol, source["quality"]["best_coords"]
        )
        oracle = prior.rdkit_ensemble(mol, 60)
        assert oracle["best_energy"] is not None
        delta = best_energy - oracle["best_energy"]
        result = {
            "name": row["name"],
            "smiles": row["smiles"],
            "quality_best_energy_kcal_mol": best_energy,
            "rdkit_60_best_energy_kcal_mol": oracle["best_energy"],
            "quality_minus_rdkit_kcal_mol": delta,
            "quality_generated": source["quality"]["generated"],
            "quality_converged": source["quality"]["converged"],
            "rdkit_generated": oracle["generated"],
            "rdkit_converged": oracle["converged"],
            "max_same_coordinate_energy_error_kcal_mol": max(candidate_errors),
            "quality_candidate_cip_violations": cip_violations,
            "quality_candidate_double_bond_violations": bond_violations,
        }
        print(f"{row['name']}: {delta:+.6f} kcal/mol", flush=True)
        results.append(result)
    tolerance = 1e-6
    deltas = [row["quality_minus_rdkit_kcal_mol"] for row in results]
    output = {
        "profile": "etkdg_quality_tier_b_followup_development_v1",
        "rdkit_version": rdBase.rdkitVersion,
        "selector": "same eligibility as third-seed batch; SHA-256 SMILES order; indices 12:33",
        "comparison": "CheMatic quality: up to 60 candidates, two optional 400-step retries, returns up to 20; RDKit: 60 candidates, 200 MMFF94 steps each",
        "hashes": {
            "manifest": sha256(manifest),
            "binary": sha256(binary),
            "followup_script": sha256(Path(__file__)),
            "prior_script": sha256(Path(prior.__file__)),
        },
        "summary": {
            "rows": len(results),
            "at_or_below_rdkit_within_1e_6": sum(x <= tolerance for x in deltas),
            "strictly_below_rdkit_by_more_than_1e_6": sum(x < -tolerance for x in deltas),
            "above_rdkit_by_more_than_1e_6": sum(x > tolerance for x in deltas),
            "largest_deficit_kcal_mol": max(deltas),
            "max_same_coordinate_energy_error_kcal_mol": max(
                row["max_same_coordinate_energy_error_kcal_mol"] for row in results
            ),
            "cip_violations": sum(row["quality_candidate_cip_violations"] for row in results),
            "double_bond_violations": sum(
                row["quality_candidate_double_bond_violations"] for row in results
            ),
        },
        "results": results,
    }
    assert all(math.isfinite(x) for x in deltas)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
