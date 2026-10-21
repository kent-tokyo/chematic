#!/usr/bin/env python3
"""Compare the ETKDG quality profile with RDKit on a fixed development slice.

The 12 rows are selected by structure and SMILES hash before any embedding or
energy result is inspected. This already-used Tier B corpus is development
evidence, not a sealed holdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"


def select_rows() -> list[dict[str, str]]:
    candidates = []
    for row in json.loads(MANIFEST.read_text())["molecules"]:
        if row["heavy_atom_count"] > 30 or row["primary_category"] == "macrocycle":
            continue
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None or not 3 <= rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10:
            continue
        if len(Chem.GetMolFrags(mol)) != 1:
            continue
        with_h = Chem.AddHs(mol)
        if not AllChem.MMFFHasAllMoleculeParams(with_h):
            continue
        torsions = rdDistGeom.GetExperimentalTorsions(with_h)
        if any(float(torsion["V"][i]) != 0 for torsion in torsions for i in (4, 5)):
            continue
        key = hashlib.sha256(row["smiles"].encode()).hexdigest()
        candidates.append((key, {"name": row["name"], "smiles": row["smiles"]}))
    return [row for _, row in sorted(candidates)[:12]]


def check_graph(source: dict, mol: Chem.Mol) -> None:
    assert source["atomic_numbers"] == [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    rust_bonds = sorted(
        (a, b, "SINGLE" if order in ("UP", "DOWN") else order)
        for a, b, order in source["bonds"]
    )
    rdkit_bonds = sorted(
        (
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            str(bond.GetBondType()).upper(),
        )
        for bond in mol.GetBonds()
    )
    assert rust_bonds == rdkit_bonds, source["name"]


def score_coords(
    mol: Chem.Mol, coords: list[list[float]]
) -> tuple[float, dict[int, str], dict[int, str]]:
    scored = Chem.Mol(mol)
    assert len(coords) == scored.GetNumAtoms()
    conformer = Chem.Conformer(scored.GetNumAtoms())
    for index, point in enumerate(coords):
        conformer.SetAtomPosition(index, point)
    scored.RemoveAllConformers()
    scored.AddConformer(conformer)
    props = AllChem.MMFFGetMoleculeProperties(scored, mmffVariant="MMFF94")
    assert props is not None
    field = AllChem.MMFFGetMoleculeForceField(scored, props)
    assert field is not None
    energy = field.CalcEnergy()
    Chem.AssignStereochemistryFrom3D(scored)
    Chem.AssignStereochemistry(scored, cleanIt=True, force=True)
    cip = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in scored.GetAtoms()
        if atom.HasProp("_CIPCode")
    }
    bond_stereo = {
        bond.GetIdx(): str(bond.GetStereo())
        for bond in scored.GetBonds()
        if bond.GetStereo() != Chem.BondStereo.STEREONONE
    }
    return energy, cip, bond_stereo


def rdkit_40(source: Chem.Mol) -> dict[str, float | int | None]:
    mol = Chem.Mol(source)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=40, params=params)
    optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=200)
    assert len(ids) == len(optimized)
    energies = [energy for _, energy in optimized if math.isfinite(energy)]
    return {
        "generated": len(ids),
        "converged": sum(status == 0 for status, _ in optimized),
        "best_energy": min(energies, default=None),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    selected = select_rows()
    assert len(selected) == 12
    payload = "".join(json.dumps(row) + "\n" for row in selected)
    run = subprocess.run(
        [str(args.binary.resolve())], input=payload, text=True,
        capture_output=True, check=True,
    )
    rust_rows = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(rust_rows) == len(selected)
    results = []
    for expected_row, rust in zip(selected, rust_rows):
        assert expected_row["name"] == rust["name"]
        assert rust["high_order_gradient_terms"] is False, rust["name"]
        mol = Chem.AddHs(Chem.MolFromSmiles(rust["smiles"]))
        check_graph(rust, mol)
        input_mol = Chem.MolFromSmiles(rust["smiles"])
        Chem.AssignStereochemistry(input_mol, cleanIt=True, force=True)
        expected_cip = {
            atom.GetIdx(): atom.GetProp("_CIPCode")
            for atom in input_mol.GetAtoms()
            if atom.HasProp("_CIPCode")
        }
        expected_bond_stereo = {
            bond.GetIdx(): str(bond.GetStereo())
            for bond in input_mol.GetBonds()
            if bond.GetStereo() != Chem.BondStereo.STEREONONE
        }
        scored = {}
        for arm in ("compatibility", "quality"):
            result = rust[arm]
            assert result["best_energy"] is not None, (rust["name"], arm)
            independent_energy, cip, bond_stereo = score_coords(mol, result["best_coords"])
            assert abs(independent_energy - result["best_energy"]) < 1e-6
            scored[arm] = {
                "generated": result["generated"],
                "converged": result["converged"],
                "best_energy": result["best_energy"],
                "same_coordinate_energy_error": independent_energy - result["best_energy"],
                "tetrahedral_cip_matches_input": all(
                    cip.get(i) == code for i, code in expected_cip.items()
                ),
                "double_bond_stereo_matches_input": all(
                    bond_stereo.get(i) == stereo
                    for i, stereo in expected_bond_stereo.items()
                ),
            }
        oracle = rdkit_40(mol)
        assert oracle["best_energy"] is not None
        assert scored["quality"]["best_energy"] <= scored["compatibility"]["best_energy"] + 1e-6
        results.append({
            **expected_row,
            "declared_stereo_centers": len(expected_cip),
            "declared_stereo_double_bonds": len(expected_bond_stereo),
            "compatibility": scored["compatibility"],
            "quality": scored["quality"],
            "rdkit_40": oracle,
            "quality_minus_rdkit": scored["quality"]["best_energy"] - oracle["best_energy"],
            "quality_minus_compatibility": scored["quality"]["best_energy"] - scored["compatibility"]["best_energy"],
        })
    deltas = [row["quality_minus_rdkit"] for row in results]
    report = {
        "profile": "etkdg_quality_seed_budget_development_v1",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": hashlib.sha256((ROOT / "crates/chematic-3d/src/rdkit_embed/molecule.rs").read_bytes()).hexdigest(),
        "rust_binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        "manifest": str(MANIFEST.relative_to(ROOT)),
        "manifest_sha256": hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
        "selection": "single-fragment non-macrocycles, <=30 heavy atoms, 3-10 RDKit rotatable bonds, MMFF parameters available, no ETKDG V5/V6 coefficients; sort by SHA-256(SMILES), take first 12",
        "rdkit_version": rdBase.rdkitVersion,
        "candidate_budget": "CheMatic quality: 20 + 20 distinct seeds; RDKit: 40; MMFF94 up to 200 iterations per candidate",
        "results": results,
        "summary": {
            "rows": len(results),
            "quality_improved_vs_prior_profile": sum(row["quality_minus_compatibility"] < -1e-6 for row in results),
            "quality_regressed_vs_prior_profile": sum(row["quality_minus_compatibility"] > 1e-6 for row in results),
            "quality_at_or_below_rdkit": sum(delta <= 1e-6 for delta in deltas),
            "median_quality_minus_rdkit": statistics.median(deltas),
            "max_quality_minus_rdkit": max(deltas),
            "tetrahedral_cip_violations": sum(
                not row["quality"]["tetrahedral_cip_matches_input"] for row in results
            ),
            "double_bond_stereo_violations": sum(
                not row["quality"]["double_bond_stereo_matches_input"] for row in results
            ),
            "max_same_coordinate_energy_error": max(
                abs(row[arm]["same_coordinate_energy_error"])
                for row in results for arm in ("compatibility", "quality")
            ),
        },
        "limits": [
            "The Tier B corpus was used in earlier development and is not an untouched holdout.",
            "Lowest MMFF94 energy is one quality metric, not experimental conformer accuracy or population.",
            "This run does not measure elapsed time, memory, conformer diversity, or public release artifacts.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
