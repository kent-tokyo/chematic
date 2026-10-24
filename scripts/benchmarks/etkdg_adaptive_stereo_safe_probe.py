#!/usr/bin/env python3
"""Compare bounded adaptive MMFF94 retries with RDKit on 12 hash-selected ChEMBL rows.

Both engines initially generate 20 ETKDGv3 conformers and optimize for at
most 200 MMFF94 iterations. CheMatic's adaptive arm alone may spend up to
400 more iterations on its lowest-energy unfinished candidate. This measures
converged energy under disclosed unequal work, not experimental accuracy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

import etkdg_quality_seed_budget_probe as quality
import etkdg_torsion_search_probe as selection

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "scripts/fixtures/chembl_accuracy_corpus_4999.smi"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score_arm(mol: Chem.Mol, arm: dict, expected_cip: dict, expected_bonds: dict) -> dict:
    coords = arm["best_converged_coords"]
    energy = arm["best_converged_energy"]
    checked = None
    if coords is not None:
        checked, cip, bonds = quality.score_coords(mol, coords)
        if not math.isfinite(checked) or abs(checked - energy) > 1e-6:
            raise ValueError("independent RDKit MMFF94 re-score disagrees")
        if any(cip.get(index) != code for index, code in expected_cip.items()):
            raise ValueError("declared tetrahedral stereo changed")
        if any(bonds.get(index) != code for index, code in expected_bonds.items()):
            raise ValueError("declared double-bond stereo changed")
    elif energy is not None:
        raise ValueError("energy without converged coordinates")
    return {
        "generated": arm["generated"],
        "converged": arm["converged"],
        "best_converged_energy": checked,
        "same_coordinate_energy_error": None if checked is None else checked - energy,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chembl-start", type=int, default=24)
    args = parser.parse_args()
    if args.chembl_start < 0:
        parser.error("--chembl-start must be nonnegative")
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    selected = selection.select_chembl(args.chembl_start)
    if len(selected) != 12:
        raise ValueError("expected exactly 12 selected molecules")
    rows = []
    for row in selected:
        payload = {**row, "num_confs": 20, "max_iters": 200, "adaptive_retry_iters": 400}
        run = subprocess.run(
            [str(args.binary.resolve())], input=json.dumps(payload) + "\n",
            text=True, capture_output=True, check=True,
        )
        source = json.loads(run.stdout)
        if source["name"] != row["name"]:
            raise ValueError("source row order changed")
        mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        quality.check_graph(source, mol)
        declared = Chem.MolFromSmiles(row["smiles"])
        Chem.AssignStereochemistry(declared, cleanIt=True, force=True)
        expected_cip = {
            atom.GetIdx(): atom.GetProp("_CIPCode") for atom in declared.GetAtoms()
            if atom.HasProp("_CIPCode")
        }
        expected_bonds = {
            bond.GetIdx(): str(bond.GetStereo()) for bond in declared.GetBonds()
            if bond.GetStereo() != Chem.BondStereo.STEREONONE
        }
        arms = {
            arm: score_arm(mol, source[arm], expected_cip, expected_bonds)
            for arm in ("compatibility", "adaptive")
        }
        oracle = quality.rdkit_ensemble(mol, 20)
        for arm in arms.values():
            value, other = arm["best_converged_energy"], oracle["best_converged_energy"]
            arm["minus_rdkit"] = None if value is None or other is None else value - other
        before = arms["compatibility"]["best_converged_energy"]
        after = arms["adaptive"]["best_converged_energy"]
        if before is not None and after is not None and after > before + 1e-6:
            raise ValueError("adaptive retry regressed the converged minimum")
        rows.append({**row, "arms": arms, "rdkit": oracle})
        print(f"{row['name']}: {arms['adaptive']['minus_rdkit']}", flush=True)
    summary = {}
    for name in ("compatibility", "adaptive"):
        deltas = [row["arms"][name]["minus_rdkit"] for row in rows
                  if row["arms"][name]["minus_rdkit"] is not None]
        summary[name] = {
            "comparable": len(deltas),
            "chematic_lower": sum(delta < -1e-6 for delta in deltas),
            "rdkit_lower": sum(delta > 1e-6 for delta in deltas),
            "ties": sum(abs(delta) <= 1e-6 for delta in deltas),
        }
    report = {
        "schema_version": 1,
        "rdkit_version": rdBase.rdkitVersion,
        "selection": f"ChEMBL eligible SMILES SHA-256 positions {args.chembl_start}-{args.chembl_start + 11}",
        "conditions": (
            "20 ETKDGv3 seed-42 candidates and MMFF94 200 iterations per engine; "
            "CheMatic adaptive arm retries one unfinished candidate for up to "
            "400 additional iterations; lowest converged energy in kcal/mol"
        ),
        "hashes": {
            "binary": sha256(args.binary),
            "script": sha256(Path(__file__)),
            "corpus": sha256(CORPUS),
            "retry_source": sha256(ROOT / "crates/chematic-3d/src/rdkit_embed/molecule/quality_retry.rs"),
        },
        "summary": summary,
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
