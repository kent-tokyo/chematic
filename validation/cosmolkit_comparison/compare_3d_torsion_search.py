#!/usr/bin/env python3
"""Compare bounded CheMatic torsion search with pinned COSMolKit MMFF94.

The 12 inputs are chosen by SMILES hash before scoring either engine. This
measures the lowest converged force-field minimum, not physical accuracy or
performance at equal computational cost.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
from pathlib import Path

import cosmolkit as ck
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"
THRESHOLD = 0.01


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_rows(offset: int = 0, limit: int = 12) -> tuple[list[int], list[dict[str, str]], int]:
    rows = json.loads(MANIFEST.read_text())["molecules"]
    eligible = []
    for index, row in enumerate(rows[48:], 48):
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None or row.get("primary_category") == "macrocycle":
            continue
        if len(Chem.GetMolFrags(mol)) != 1:
            continue
        heavy = mol.GetNumHeavyAtoms()
        rotors = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if not (8 <= heavy <= 25 and 3 <= rotors <= 8):
            continue
        if not AllChem.MMFFHasAllMoleculeParams(Chem.AddHs(mol)):
            continue
        eligible.append((hashlib.sha256(row["smiles"].encode()).hexdigest(), index, row))
    eligible.sort()
    if len(eligible) != 35:
        raise ValueError(f"the pinned eligible count changed: {len(eligible)}")
    chosen = eligible[offset:offset + limit]
    return [index for _, index, _ in chosen], [row for _, _, row in chosen], len(eligible)


def graph_check(source: dict, cos_mol: ck.Molecule, rd_mol: Chem.Mol) -> None:
    rd_atoms = [
        (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope())
        for atom in rd_mol.GetAtoms()
    ]
    cos_atoms = [
        (atom.atomic_number(), atom.formal_charge(), atom.isotope() or 0)
        for atom in cos_mol.atoms()
    ]
    if cos_atoms != rd_atoms:
        raise ValueError(f"COSMolKit indexed atoms differ: {source['name']}")
    if source["atomic_numbers"] != [atom[0] for atom in rd_atoms]:
        raise ValueError(f"CheMatic indexed atoms differ: {source['name']}")
    rd_bonds = sorted(
        (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
         max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
         str(b.GetBondType()).upper())
        for b in rd_mol.GetBonds()
    )
    cos_bonds = sorted(
        (min(b.begin(), b.end()), max(b.begin(), b.end()), b.order().name)
        for b in cos_mol.bonds()
    )
    chem_bonds = sorted(
        (a, b, "SINGLE" if order in ("UP", "DOWN") else order)
        for a, b, order in source["bonds"]
    )
    if cos_bonds != rd_bonds or chem_bonds != rd_bonds:
        raise ValueError(f"indexed bonds differ: {source['name']}")


def score_coords(
    rd_mol: Chem.Mol, coords: list[list[float]]
) -> tuple[float, dict[int, str], dict[int, str]]:
    mol = Chem.Mol(rd_mol)
    if len(coords) != mol.GetNumAtoms() or not all(len(point) == 3 for point in coords) or not all(
        math.isfinite(value) for point in coords for value in point
    ):
        raise ValueError("invalid coordinates")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for index, point in enumerate(coords):
        conf.SetAtomPosition(index, point)
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    if props is None:
        raise ValueError("RDKit MMFF parameters unavailable")
    field = AllChem.MMFFGetMoleculeForceField(
        mol, props, nonBondedThresh=100.0
    )
    if field is None:
        raise ValueError("RDKit MMFF field unavailable")
    energy = field.CalcEnergy()
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    cip = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in mol.GetAtoms() if atom.HasProp("_CIPCode")
    }
    bonds = {
        bond.GetIdx(): str(bond.GetStereo())
        for bond in mol.GetBonds()
        if bond.GetStereo() != Chem.BondStereo.STEREONONE
    }
    return energy, cip, bonds


def check_stereo(
    expected_cip: dict[int, str], expected_bonds: dict[int, str],
    actual_cip: dict[int, str], actual_bonds: dict[int, str],
) -> None:
    if any(actual_cip.get(i) != value for i, value in expected_cip.items()):
        raise ValueError("declared tetrahedral stereo changed")
    if any(actual_bonds.get(i) != value for i, value in expected_bonds.items()):
        raise ValueError("declared double-bond stereo changed")


def score_source_arm(
    arm: dict, mol: Chem.Mol,
    expected_cip: dict[int, str], expected_bonds: dict[int, str],
) -> dict:
    energy = arm["best_converged_energy"]
    coords = arm["best_converged_coords"]
    if energy is None or coords is None:
        raise ValueError("no converged CheMatic result")
    independently_scored, cip, bonds = score_coords(mol, coords)
    check_stereo(expected_cip, expected_bonds, cip, bonds)
    error = independently_scored - energy
    if abs(error) > 1e-6:
        raise ValueError(f"CheMatic energy disagrees with RDKit: {error}")
    return {
        "returned": arm["returned"],
        "converged": arm["converged"],
        "best_converged_energy": independently_scored,
        "rdkit_rescore_error": error,
    }


def score_cosmolkit(
    cos_mol: ck.Molecule, rd_mol: Chem.Mol,
    expected_cip: dict[int, str], expected_bonds: dict[int, str],
) -> dict:
    params = ck.EmbedParams.etkdg_v3()
    params.random_seed = 42
    params.num_threads = 1
    embedded = cos_mol.with_3d_conformers_result(24, params).molecule()
    result = embedded.with_mmff_optimized_conformers(
        num_threads=1, max_iterations=1000, non_bonded_threshold=100.0
    )
    outcomes = result.conformer_results()
    conformers = result.molecule().conformers_3d()
    if len(outcomes) != len(conformers):
        raise ValueError("COSMolKit result/conformer count differs")
    converged = [
        i for i, outcome in enumerate(outcomes)
        if outcome.status_code() == 0 and math.isfinite(outcome.energy())
    ]
    if not converged:
        raise ValueError("no converged COSMolKit result")
    best = min(converged, key=lambda index: outcomes[index].energy())
    coords = conformers[best].coordinates().tolist()
    energy, cip, bonds = score_coords(rd_mol, coords)
    check_stereo(expected_cip, expected_bonds, cip, bonds)
    error = energy - outcomes[best].energy()
    if abs(error) > 1e-6:
        raise ValueError(f"COSMolKit energy disagrees with RDKit: {error}")
    return {
        "returned": len(outcomes),
        "converged": len(converged),
        "best_converged_energy": energy,
        "rdkit_rescore_error": error,
    }


def classify(source: float, cosmolkit: float) -> str:
    delta = source - cosmolkit
    if delta < -THRESHOLD:
        return "chematic_lower"
    if delta > THRESHOLD:
        return "cosmolkit_lower"
    return "within_0.01"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-output", required=True, type=Path)
    parser.add_argument("--source-results", type=Path)
    parser.add_argument("--selection-offset", type=int, default=0)
    parser.add_argument("--selection-limit", type=int, default=12)
    parser.add_argument("--candidate-confs", type=int,
                        help="CheMatic seed-batch candidate count; returns 24")
    args = parser.parse_args()
    if args.selection_offset < 0 or args.selection_limit <= 0:
        parser.error("selection offset must be nonnegative and limit positive")
    if args.candidate_confs is not None and args.candidate_confs < 24:
        parser.error("candidate conformers must be at least the 24 returned conformers")
    if ck.__version__ != "0.5.0" or rdBase.rdkitVersion != "2026.09.1":
        parser.error("expected COSMolKit 0.5.0 and RDKit 2026.09.1")
    RDLogger.DisableLog("rdApp.*")
    indices, selected, eligible_count = selected_rows(
        args.selection_offset, args.selection_limit
    )
    if len(selected) != args.selection_limit:
        parser.error("selection exceeds the eligible rows")
    payload = "".join(
        json.dumps({
            "name": row["name"], "smiles": row["smiles"],
            "num_confs": 24, "max_iters": 1000,
            **({"candidate_confs": args.candidate_confs}
               if args.candidate_confs is not None else {}),
        }) + "\n"
        for row in selected
    )
    if args.source_results is None:
        run = subprocess.run(
            [str(args.binary.resolve())], input=payload,
            text=True, capture_output=True, check=True,
        )
        source_text = run.stdout
    else:
        source_text = args.source_results.read_text()
    raw_rows = [json.loads(line) for line in source_text.splitlines()]
    if len(raw_rows) != len(selected):
        raise ValueError("CheMatic source row count differs")
    rows = []
    for expected, source in zip(selected, raw_rows, strict=True):
        if source["name"] != expected["name"] or source["smiles"] != expected["smiles"]:
            raise ValueError("CheMatic source order differs")
        if source["num_confs"] != 24 or source["max_iters"] != 1000:
            raise ValueError("CheMatic settings differ")
        if source.get("candidate_confs") != args.candidate_confs:
            raise ValueError("CheMatic candidate count differs")
        input_mol = Chem.MolFromSmiles(expected["smiles"])
        rd_mol = Chem.AddHs(input_mol)
        cos_mol = ck.Molecule.from_smiles(expected["smiles"]).with_hydrogens()
        graph_check(source, cos_mol, rd_mol)
        Chem.AssignStereochemistry(input_mol, cleanIt=True, force=True)
        expected_cip = {
            atom.GetIdx(): atom.GetProp("_CIPCode")
            for atom in input_mol.GetAtoms() if atom.HasProp("_CIPCode")
        }
        expected_bonds = {
            bond.GetIdx(): str(bond.GetStereo())
            for bond in input_mol.GetBonds()
            if bond.GetStereo() != Chem.BondStereo.STEREONONE
        }
        quality = score_source_arm(
            source["quality"], rd_mol, expected_cip, expected_bonds
        )
        torsion = score_source_arm(
            source["torsion"], rd_mol, expected_cip, expected_bonds
        )
        pool = (
            score_source_arm(source["pool"], rd_mol, expected_cip, expected_bonds)
            if args.candidate_confs is not None else None
        )
        if torsion["best_converged_energy"] > quality["best_converged_energy"] + 1e-6:
            raise ValueError("torsion search regressed the quality minimum")
        cosmolkit = score_cosmolkit(cos_mol, rd_mol, expected_cip, expected_bonds)
        rows.append({
            "name": expected["name"], "smiles": expected["smiles"],
            "quality": quality, "torsion": torsion, "cosmolkit": cosmolkit,
            **({"pool": pool,
                "pool_vs_cosmolkit": classify(
                    pool["best_converged_energy"],
                    cosmolkit["best_converged_energy"],
                )} if pool is not None else {}),
            "quality_vs_cosmolkit": classify(
                quality["best_converged_energy"],
                cosmolkit["best_converged_energy"],
            ),
            "torsion_vs_cosmolkit": classify(
                torsion["best_converged_energy"],
                cosmolkit["best_converged_energy"],
            ),
        })
        print(
            f"{expected['name']}: {rows[-1]['quality_vs_cosmolkit']} -> "
            f"{rows[-1]['torsion_vs_cosmolkit']}"
            + (f" -> {rows[-1]['pool_vs_cosmolkit']}" if pool is not None else ""),
            flush=True,
        )
    args.source_output.parent.mkdir(parents=True, exist_ok=True)
    args.source_output.write_text(source_text)
    summary = {
        arm: {
            label: sum(row[f"{arm}_vs_cosmolkit"] == label for row in rows)
            for label in ("chematic_lower", "cosmolkit_lower", "within_0.01")
        }
        for arm in (("quality", "torsion", "pool") if args.candidate_confs is not None
                    else ("quality", "torsion"))
    }
    extension = next(Path(ck.__file__).parent.glob("cosmolkit*.so"))
    report = {
        "schema_version": 1,
        "metric": "lowest converged MMFF94 kcal/mol; 0.01 kcal/mol threshold",
        "conditions": (
            "explicit H; ETKDGv3 seed 42; one COSMolKit thread; "
            "24 conformers per seed batch; 1000 MMFF iterations; "
            "100 A nonbonded threshold; CheMatic torsion arm up to "
            "24 additional MMFF minimizations"
        ),
        "candidate_confs": args.candidate_confs,
        "selection": {
            "source_manifest": str(MANIFEST.relative_to(ROOT)),
            "manifest_sha256": sha256(MANIFEST),
            "source_indices_at_least": 48,
            "single_fragment": True,
            "exclude_macrocycles": True,
            "heavy_atoms": [8, 25],
            "rotatable_bonds": [3, 8],
            "mmff94_parameters_required": True,
            "eligible": eligible_count,
            "order": "SHA-256 of source SMILES ascending",
            "offset": args.selection_offset,
            "limit": args.selection_limit,
            "selected_indices": indices,
        },
        "versions": {
            "cosmolkit": ck.__version__,
            "rdkit": rdBase.rdkitVersion,
            "python": platform.python_version(),
            "host": platform.platform(),
        },
        "hashes": {
            "binary": sha256(args.binary),
            "script": sha256(Path(__file__)),
            "source_stream": sha256(args.source_output),
            "torsion_source": sha256(
                ROOT / "crates/chematic-3d/src/rdkit_embed/molecule/torsion_search.rs"
            ),
            "cosmolkit_extension": sha256(extension),
        },
        "summary": summary,
        "torsion_improved_over_quality": sum(
            row["torsion"]["best_converged_energy"]
            < row["quality"]["best_converged_energy"] - 1e-6
            for row in rows
        ),
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
