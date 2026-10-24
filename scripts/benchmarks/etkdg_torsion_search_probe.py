#!/usr/bin/env python3
"""Compare bounded torsion moves with RDKit 2026.09.1 ETKDGv3/MMFF94.

Tier B's 21 rows are prior development examples. Tier A's nine eligible rows
are disjoint but also come from an established repository validation manifest.
The ChEMBL subset is selected by SMILES hash before either engine is run.
Minimum converged MMFF94 energy measures search within one model, not
experimental structure quality. CheMatic's torsion arm spends extra work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors

import etkdg_quality_seed_budget_probe as prior


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_tier_a() -> list[dict[str, str]]:
    manifest = ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"
    selected = []
    for row in json.loads(manifest.read_text())["molecules"]:
        mol = Chem.MolFromSmiles(row["smiles"])
        if (mol is None or mol.GetNumHeavyAtoms() > 30
                or row["primary_category"] == "macrocycle"
                or not 3 <= rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
                or len(Chem.GetMolFrags(mol)) != 1):
            continue
        with_h = Chem.AddHs(mol)
        if not AllChem.MMFFHasAllMoleculeParams(with_h):
            continue
        torsions = rdDistGeom.GetExperimentalTorsions(with_h)
        if any(float(torsion["V"][i]) != 0 for torsion in torsions for i in (4, 5)):
            continue
        selected.append({"name": row["name"], "smiles": row["smiles"]})
    return sorted(selected, key=lambda row: hashlib.sha256(row["smiles"].encode()).hexdigest())


def select_chembl() -> list[dict[str, str]]:
    known = {
        row["smiles"] for tier in ("a", "b")
        for row in json.loads((ROOT / f"validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{tier}.json").read_text())["molecules"]
    }
    corpus = ROOT / "scripts/fixtures/chembl_accuracy_corpus_4999.smi"
    selected = []
    for index, smiles in enumerate(corpus.read_text().splitlines()):
        mol = Chem.MolFromSmiles(smiles)
        if (not smiles or smiles in known or mol is None
                or mol.GetNumHeavyAtoms() > 25
                or not 3 <= rdMolDescriptors.CalcNumRotatableBonds(mol) <= 8
                or len(Chem.GetMolFrags(mol)) != 1):
            continue
        with_h = Chem.AddHs(mol)
        if not AllChem.MMFFHasAllMoleculeParams(with_h):
            continue
        torsions = rdDistGeom.GetExperimentalTorsions(with_h)
        if any(float(torsion["V"][i]) != 0 for torsion in torsions for i in (4, 5)):
            continue
        selected.append({"name": f"chembl_corpus_{index:04d}", "smiles": smiles})
    return sorted(selected, key=lambda row: hashlib.sha256(row["smiles"].encode()).hexdigest())[:12]


def score_arm(mol: Chem.Mol, arm: dict, expected_cip: dict, expected_bonds: dict) -> dict:
    coords = arm["best_converged_coords"]
    if coords is None:
        return {"returned": arm["returned"], "converged": arm["converged"],
                "best_converged_energy": None}
    energy, cip, bonds = prior.score_coords(mol, coords)
    if abs(energy - arm["best_converged_energy"]) > 1e-6:
        raise ValueError("independent RDKit MMFF energy disagrees")
    if any(cip.get(index) != code for index, code in expected_cip.items()):
        raise ValueError("declared tetrahedral stereo changed")
    if any(bonds.get(index) != code for index, code in expected_bonds.items()):
        raise ValueError("declared double-bond stereo changed")
    return {"returned": arm["returned"], "converged": arm["converged"],
            "best_converged_energy": energy,
            "same_coordinate_energy_error": energy - arm["best_converged_energy"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tier", choices=("a", "b", "chembl"), required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    selected = (select_tier_a() if args.tier == "a" else
                prior.select_rows(start=12, count=21) if args.tier == "b" else
                select_chembl())
    assert len(selected) == {"a": 9, "b": 21, "chembl": 12}[args.tier]
    rows = []
    for row in selected:
        run = subprocess.run([str(args.binary.resolve())], input=json.dumps(row) + "\n",
                             text=True, capture_output=True, check=True)
        source = json.loads(run.stdout)
        assert source["name"] == row["name"]
        mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        prior.check_graph(source, mol)
        declared = Chem.MolFromSmiles(row["smiles"])
        Chem.AssignStereochemistry(declared, cleanIt=True, force=True)
        expected_cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in declared.GetAtoms()
                        if a.HasProp("_CIPCode")}
        expected_bonds = {b.GetIdx(): str(b.GetStereo()) for b in declared.GetBonds()
                          if b.GetStereo() != Chem.BondStereo.STEREONONE}
        scored = {arm: score_arm(mol, source[arm], expected_cip, expected_bonds)
                  for arm in ("quality", "torsion")}
        rdkit = prior.rdkit_ensemble(mol, 60)
        baseline = scored["quality"]["best_converged_energy"]
        improved = scored["torsion"]["best_converged_energy"]
        if baseline is not None and improved is not None and improved > baseline + 1e-6:
            raise ValueError("torsion search regressed the converged minimum")
        delta = improved - rdkit["best_converged_energy"] if (
            improved is not None and rdkit["best_converged_energy"] is not None) else None
        print(f"{row['name']}: {delta:+.6f}" if delta is not None else
              f"{row['name']}: incomparable", flush=True)
        rows.append({**row, "quality": scored["quality"], "torsion": scored["torsion"],
                     "rdkit": rdkit, "torsion_minus_rdkit": delta})
    differences = [row["torsion_minus_rdkit"] for row in rows
                   if row["torsion_minus_rdkit"] is not None]
    report = {
        "schema_version": 1,
        "rdkit_version": rdBase.rdkitVersion,
        "tier": args.tier,
        "selection": "Tier A all nine eligible; Tier B SHA-256 SMILES positions 12:33; ChEMBL 12 lowest SHA-256 among 851 RDKit-eligible molecules excluding Tier A/B SMILES",
        "conditions": "ETKDGv3 seed 42; 20 CheMatic candidates per seed batch; RDKit 60 candidates; MMFF94 200 iterations; torsion arm adds at most 24 minimizations",
        "hashes": {
            "binary": sha256(args.binary),
            "script": sha256(Path(__file__)),
            "corpus": sha256(ROOT / ("scripts/fixtures/chembl_accuracy_corpus_4999.smi"
                if args.tier == "chembl" else
                f"validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{args.tier}.json")),
            "torsion_source": sha256(ROOT / "crates/chematic-3d/src/rdkit_embed/molecule/torsion_search.rs"),
        },
        "summary": {
            "rows": len(rows),
            "comparable": len(differences),
            "torsion_lower_than_rdkit_by_1e_6": sum(value < -1e-6 for value in differences),
            "rdkit_lower_than_torsion_by_1e_6": sum(value > 1e-6 for value in differences),
            "within_1e_6": sum(abs(value) <= 1e-6 for value in differences),
            "torsion_improved_over_quality": sum(
                row["quality"]["best_converged_energy"] is not None
                and row["torsion"]["best_converged_energy"] is not None
                and row["torsion"]["best_converged_energy"]
                < row["quality"]["best_converged_energy"] - 1e-6 for row in rows),
        },
        "rows": rows,
    }
    if not all(math.isfinite(value) for value in differences):
        raise ValueError("nonfinite comparison")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
