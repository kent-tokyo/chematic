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
    converged_energies = [
        energy for status, energy in optimized
        if status == 0 and math.isfinite(energy)
    ]
    return {
        "generated": len(ids),
        "converged": sum(status == 0 for status, _ in optimized),
        "best_energy": min(energies, default=None),
        "best_converged_energy": min(converged_energies, default=None),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--prior-report", type=Path,
        default=ROOT / "validation/results/etkdg-quality-seed-budget-2026-10-22.json",
    )
    parser.add_argument("--profile", default="etkdg_quality_pooled_development_v2")
    args = parser.parse_args()
    prior_path = args.prior_report.resolve()
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
    prior = json.loads(prior_path.read_text())
    prior_rows = {row["name"]: row for row in prior["results"]}
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
                "best_converged_energy": result["best_converged_energy"],
                "tetrahedral_cip_matches_input": all(
                    cip.get(i) == code for i, code in expected_cip.items()
                ),
                "double_bond_stereo_matches_input": all(
                    bond_stereo.get(i) == stereo
                    for i, stereo in expected_bond_stereo.items()
                ),
            }
            if result["best_converged_coords"] is not None:
                converged_energy, _, _ = score_coords(mol, result["best_converged_coords"])
                assert abs(converged_energy - result["best_converged_energy"]) < 1e-6
        oracle = rdkit_40(mol)
        assert oracle["best_energy"] is not None
        assert scored["quality"]["best_energy"] <= scored["compatibility"]["best_energy"] + 1e-6
        candidate_checks = [
            (candidate, score_coords(mol, candidate["coords"]))
            for candidate in rust["quality_candidates"]
        ]
        assert len(candidate_checks) == scored["quality"]["generated"]
        all_cip_violations = 0
        all_double_bond_violations = 0
        all_energy_errors = []
        for candidate, (independent_energy, cip, bond_stereo) in candidate_checks:
            all_energy_errors.append(abs(independent_energy - candidate["energy"]))
            all_cip_violations += not all(
                cip.get(i) == code for i, code in expected_cip.items()
            )
            all_double_bond_violations += not all(
                bond_stereo.get(i) == stereo
                for i, stereo in expected_bond_stereo.items()
            )
        assert max(all_energy_errors, default=0.0) < 1e-6
        previous = prior_rows[rust["name"]]["quality"]
        assert scored["quality"]["best_energy"] <= previous["best_energy"] + 1e-6
        results.append({
            **expected_row,
            "declared_stereo_centers": len(expected_cip),
            "declared_stereo_double_bonds": len(expected_bond_stereo),
            "compatibility": scored["compatibility"],
            "quality": scored["quality"],
            "rdkit_40": oracle,
            "quality_minus_rdkit": scored["quality"]["best_energy"] - oracle["best_energy"],
            "quality_minus_compatibility": scored["quality"]["best_energy"] - scored["compatibility"]["best_energy"],
            "prior_quality_converged": previous["converged"],
            "prior_quality_best_energy": previous["best_energy"],
            "all_candidate_cip_violations": all_cip_violations,
            "all_candidate_double_bond_violations": all_double_bond_violations,
            "all_candidate_max_same_coordinate_energy_error": max(all_energy_errors, default=0.0),
        })
    deltas = [row["quality_minus_rdkit"] for row in results]
    report = {
        "profile": args.profile,
        "candidate_base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "candidate_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain", "--", "crates/chematic-3d", "scripts/benchmarks/etkdg_quality_seed_budget_probe.py"],
            cwd=ROOT, text=True,
        ).strip()),
        "source_sha256": hashlib.sha256((ROOT / "crates/chematic-3d/src/rdkit_embed/molecule.rs").read_bytes()).hexdigest(),
        "quality_retry_source_sha256": hashlib.sha256((ROOT / "crates/chematic-3d/src/rdkit_embed/molecule/quality_retry.rs").read_bytes()).hexdigest(),
        "probe_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "rust_binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        "manifest": str(MANIFEST.relative_to(ROOT)),
        "manifest_sha256": hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
        "prior_report": str(prior_path.relative_to(ROOT)),
        "prior_report_sha256": hashlib.sha256(prior_path.read_bytes()).hexdigest(),
        "selection": "single-fragment non-macrocycles, <=30 heavy atoms, 3-10 RDKit rotatable bonds, MMFF parameters available, no ETKDG V5/V6 coefficients; sort by SHA-256(SMILES), take first 12",
        "rdkit_version": rdBase.rdkitVersion,
        "candidate_budget": (
            "CheMatic quality: 20 + 20 distinct seeds, MMFF94 up to 200 iterations per "
            "candidate; up to 400 extra iterations for the lowest finite-energy "
            "unfinished returned candidate even if the overall minimum converged; "
            "RDKit: 40 seeds, MMFF94 up to 200 iterations per candidate"
            if args.profile.endswith("_v5") else
            "CheMatic quality: 20 + 20 distinct seeds, MMFF94 up to 200 iterations per "
            "candidate; if the lowest-energy candidate did not converge, up to 400 "
            "extra iterations for that candidate; RDKit: 40 seeds, MMFF94 up to "
            "200 iterations per candidate"
            if args.profile.endswith("_v4") else
            "CheMatic quality: 20 + 20 distinct seeds, MMFF94 up to 200 iterations per "
            "candidate; if none converges, up to 400 extra iterations for the lowest "
            "finite-energy candidate; RDKit: 40 seeds, MMFF94 up to 200 iterations "
            "per candidate"
            if args.profile.endswith("_v3") else
            "CheMatic quality: 20 + 20 distinct seeds; RDKit: 40; MMFF94 up to "
            "200 iterations per candidate"
        ),
        "results": results,
        "summary": {
            "rows": len(results),
            "quality_improved_vs_compatibility": sum(row["quality_minus_compatibility"] < -1e-6 for row in results),
            "quality_regressed_vs_compatibility": sum(row["quality_minus_compatibility"] > 1e-6 for row in results),
            "quality_improved_vs_prior_profile": sum(
                row["quality"]["best_energy"] < row["prior_quality_best_energy"] - 1e-6
                for row in results
            ),
            "quality_at_or_below_rdkit": sum(delta <= 1e-6 for delta in deltas),
            "median_quality_minus_rdkit": statistics.median(deltas),
            "max_quality_minus_rdkit": max(deltas),
            "quality_more_converged_than_prior": sum(
                row["quality"]["converged"] > row["prior_quality_converged"] for row in results
            ),
            "quality_fewer_converged_than_prior": sum(
                row["quality"]["converged"] < row["prior_quality_converged"] for row in results
            ),
            "quality_best_energy_regressed_vs_prior": sum(
                row["quality"]["best_energy"] > row["prior_quality_best_energy"] + 1e-6
                for row in results
            ),
            "quality_better_converged_energy_than_rdkit": sum(
                row["quality"]["best_converged_energy"] is not None
                and row["rdkit_40"]["best_converged_energy"] is not None
                and row["quality"]["best_converged_energy"]
                < row["rdkit_40"]["best_converged_energy"] - 1e-6
                for row in results
            ),
            "quality_converged_when_rdkit_did_not": sum(
                row["quality"]["best_converged_energy"] is not None
                and row["rdkit_40"]["best_converged_energy"] is None
                for row in results
            ),
            "tetrahedral_cip_violations": sum(
                not row["quality"]["tetrahedral_cip_matches_input"] for row in results
            ),
            "double_bond_stereo_violations": sum(
                not row["quality"]["double_bond_stereo_matches_input"] for row in results
            ),
            "all_candidate_cip_violations": sum(
                row["all_candidate_cip_violations"] for row in results
            ),
            "all_candidate_double_bond_violations": sum(
                row["all_candidate_double_bond_violations"] for row in results
            ),
            "max_all_candidate_same_coordinate_energy_error": max(
                row["all_candidate_max_same_coordinate_energy_error"] for row in results
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
            "The optional MMFF retry adds work, so its convergence result is not an equal-compute comparison.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
