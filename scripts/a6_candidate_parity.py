#!/usr/bin/env python3
"""Diagnose per-candidate ETKDG/MMFF differences against pinned RDKit.

The script uses one manifest row, the same explicit-hydrogen atom order, seed
schedule, conformer count and MMFF iteration budget in both engines.  It emits
initial-coordinate and optimized-energy differences without changing the
frozen corpus or selecting a preferred candidate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import chematic
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from public_package_3d_chematic import RANDOM_SEED


def aligned_rmsd(left: list[list[float]], right: list[list[float]]) -> float:
    left_array = np.asarray(left, dtype=float)
    right_array = np.asarray(right, dtype=float)
    left_array -= left_array.mean(axis=0)
    right_array -= right_array.mean(axis=0)
    covariance = left_array.T @ right_array
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = u_matrix @ vt_matrix
    if np.linalg.det(rotation) < 0:
        u_matrix[:, -1] *= -1
        rotation = u_matrix @ vt_matrix
    delta = left_array @ rotation - right_array
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def stage_profile(
    chematic_molecule: object,
    rdkit_source: Chem.Mol,
    conformers: int,
    *,
    use_exp_torsions: bool,
    use_basic_knowledge: bool,
) -> dict[str, object]:
    chematic_coords = chematic_molecule.rdkit_embed_multiple_confs(
        conformers,
        random_seed=RANDOM_SEED,
        use_exp_torsion_angle_prefs=use_exp_torsions,
        use_basic_knowledge=use_basic_knowledge,
    )
    rdkit_molecule = Chem.Mol(rdkit_source)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = RANDOM_SEED
    parameters.numThreads = 0
    parameters.useExpTorsionAnglePrefs = use_exp_torsions
    parameters.useBasicKnowledge = use_basic_knowledge
    identifiers = list(
        AllChem.EmbedMultipleConfs(
            rdkit_molecule, numConfs=conformers, params=parameters
        )
    )
    rdkit_coords = [
        [list(position) for position in rdkit_molecule.GetConformer(identifier).GetPositions()]
        for identifier in identifiers
    ]
    rmsds = [
        aligned_rmsd(left, right)
        for left, right in zip(chematic_coords, rdkit_coords)
    ]
    return {
        "use_exp_torsion_angle_prefs": use_exp_torsions,
        "use_basic_knowledge": use_basic_knowledge,
        "generated": {
            "chematic": len(chematic_coords),
            "rdkit": len(rdkit_coords),
        },
        "aligned_rmsd": rmsds,
        "max_aligned_rmsd": max(rmsds, default=None),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--row", type=int, required=True)
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--max-mmff-iters", type=int, default=200)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    source = manifest["molecules"][args.row]
    smiles = str(source["smiles"])

    chematic_molecule = chematic.from_smiles(smiles).add_hydrogens()
    chematic_initial = chematic_molecule.rdkit_embed_multiple_confs(
        args.conformers, random_seed=RANDOM_SEED
    )
    chematic_optimized = chematic_molecule.rdkit_embed_multiple_confs_mmff(
        args.conformers,
        random_seed=RANDOM_SEED,
        max_mmff_iters=args.max_mmff_iters,
    )

    rdkit_molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = RANDOM_SEED
    parameters.numThreads = 0
    identifiers = list(
        AllChem.EmbedMultipleConfs(
            rdkit_molecule, numConfs=args.conformers, params=parameters
        )
    )
    rdkit_initial = [
        [list(position) for position in rdkit_molecule.GetConformer(identifier).GetPositions()]
        for identifier in identifiers
    ]
    rdkit_optimized = AllChem.MMFFOptimizeMoleculeConfs(
        rdkit_molecule, numThreads=0, maxIters=args.max_mmff_iters
    )

    rows = []
    for index, (left, right, c_opt, r_opt) in enumerate(
        zip(chematic_initial, rdkit_initial, chematic_optimized, rdkit_optimized)
    ):
        seed = RANDOM_SEED * (index + 1)
        chematic_single = chematic_molecule.rdkit_embed(random_seed=seed)
        rdkit_single_molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        single_status = AllChem.EmbedMolecule(rdkit_single_molecule, randomSeed=seed)
        rdkit_single = (
            [list(position) for position in rdkit_single_molecule.GetConformer().GetPositions()]
            if single_status >= 0
            else None
        )
        rows.append(
            {
                "candidate": index,
                "seed": seed,
                "initial_max_abs_coordinate_delta": max(
                    abs(a - b)
                    for left_atom, right_atom in zip(left, right)
                    for a, b in zip(left_atom, right_atom)
                ),
                "initial_aligned_rmsd": aligned_rmsd(left, right),
                "chematic_multi_vs_single_aligned_rmsd": aligned_rmsd(
                    left, chematic_single
                ),
                "rdkit_multi_vs_single_aligned_rmsd": (
                    aligned_rmsd(right, rdkit_single)
                    if rdkit_single is not None
                    else None
                ),
                "single_engine_aligned_rmsd": (
                    aligned_rmsd(chematic_single, rdkit_single)
                    if rdkit_single is not None
                    else None
                ),
                "chematic_status": c_opt[0],
                "rdkit_status": r_opt[0],
                "chematic_energy": c_opt[1],
                "rdkit_energy": r_opt[1],
                "energy_delta_kcal_mol": c_opt[1] - r_opt[1],
            }
        )

    report = {
        "schema": "a6-candidate-parity/v1",
        "manifest": str(args.manifest),
        "corpus_sha256": manifest.get("corpus_sha256"),
        "row": args.row,
        "name": source.get("name"),
        "smiles": smiles,
        "chematic": getattr(chematic, "__version__", None),
        "rdkit": rdBase.rdkitVersion,
        "random_seed": RANDOM_SEED,
        "conformers": args.conformers,
        "max_mmff_iters": args.max_mmff_iters,
        "generated": {
            "chematic": len(chematic_initial),
            "rdkit": len(rdkit_initial),
        },
        "stage_profiles": [
            stage_profile(
                chematic_molecule,
                Chem.AddHs(Chem.MolFromSmiles(smiles)),
                args.conformers,
                use_exp_torsions=use_exp_torsions,
                use_basic_knowledge=use_basic_knowledge,
            )
            for use_exp_torsions, use_basic_knowledge in (
                (False, False),
                (False, True),
                (True, True),
            )
        ],
        "rows": rows,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
