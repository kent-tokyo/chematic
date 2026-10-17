"""Compare explicit-coordinate batch fits with RDKit 2026.09.1."""

import hashlib
import json
from pathlib import Path

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from chematic import align_conformers_checked


def molecule(coordinates):
    builder = Chem.RWMol()
    for _ in coordinates[0]:
        builder.AddAtom(Chem.Atom(6))
    mol = builder.GetMol()
    for frame in coordinates:
        conf = Chem.Conformer(len(frame))
        for i, xyz in enumerate(frame):
            conf.SetAtomPosition(i, xyz)
        mol.AddConformer(conf, assignId=True)
    return mol


def main():
    assert rdkit.__version__ == "2026.09.1"
    rng = np.random.default_rng(9601)
    cases = frames_compared = rms_compared = 0
    max_coordinate_difference = max_rms_difference = 0.0
    checksum = hashlib.sha256()
    for n_atoms in [4, 7, 16, 64]:
        for n_frames in [1, 2, 4, 7]:
            for sample in range(6):
                coordinates = rng.normal(size=(n_frames, n_atoms, 3))
                coordinates += rng.uniform(-5.0, 5.0, (n_frames, 1, 3))
                checksum.update(coordinates.tobytes())
                for atom_ids in [[], list(range(n_atoms - 1, -1, -1)), [0, 1, 2, 3]]:
                    selected_atoms = atom_ids or list(range(n_atoms))
                    for conf_ids in [
                        [],
                        list(range(n_frames - 1, -1, -1)),
                        [n_frames - 1],
                    ]:
                        for weights in [
                            [],
                            [1.0] * len(selected_atoms),
                            rng.uniform(0.1, 3.0, len(selected_atoms)).tolist(),
                        ]:
                            reference = molecule(coordinates)
                            rms_reference = [12.5]
                            rdMolAlign.AlignMolConformers(
                                reference,
                                atomIds=atom_ids,
                                confIds=conf_ids,
                                weights=weights,
                                RMSlist=rms_reference,
                            )
                            rms_candidate = [12.5]
                            aligned = align_conformers_checked(
                                coordinates,
                                atom_ids=atom_ids,
                                conformer_ids=conf_ids,
                                weights=weights,
                                rms_list=rms_candidate,
                            )
                            expected = np.array(
                                [
                                    conf.GetPositions()
                                    for conf in reference.GetConformers()
                                ]
                            )
                            # SVD and RDKit's iterative solver can choose slightly
                            # different rotations; coordinates are diagnostic.
                            np.testing.assert_allclose(
                                rms_candidate, rms_reference, rtol=1e-7, atol=1e-7
                            )
                            max_coordinate_difference = max(
                                max_coordinate_difference,
                                float(np.abs(aligned - expected).max()),
                            )
                            max_rms_difference = max(
                                max_rms_difference,
                                float(
                                    np.abs(
                                        np.array(rms_candidate) - rms_reference
                                    ).max()
                                ),
                            )
                            # Every frame receives a rigid transform, retaining all pair distances.
                            np.testing.assert_allclose(
                                np.linalg.norm(aligned[:, 1:] - aligned[:, :1], axis=2),
                                np.linalg.norm(
                                    coordinates[:, 1:] - coordinates[:, :1], axis=2
                                ),
                                rtol=0,
                                atol=1e-12,
                            )
                            frames_compared += n_frames
                            rms_compared += len(rms_candidate) - 1
                            cases += 1
    degenerate = []
    for n_atoms in [1, 2, 3]:
        source = rng.normal(size=(3, n_atoms, 3))
        mol = molecule(source)
        expected_rms = []
        rdMolAlign.AlignMolConformers(mol, RMSlist=expected_rms)
        observed_rms = []
        result = align_conformers_checked(source, rms_list=observed_rms)
        np.testing.assert_allclose(observed_rms, expected_rms, atol=1e-7)
        np.testing.assert_allclose(
            np.linalg.norm(result[:, 1:] - result[:, :1], axis=2),
            np.linalg.norm(source[:, 1:] - source[:, :1], axis=2),
            rtol=0,
            atol=1e-12,
        )
        degenerate.append(
            {"atoms": n_atoms, "rms_agrees": True, "coordinate_equality_claimed": False}
        )
    print(
        json.dumps(
            {
                "rdkit": rdkit.__version__,
                "seed": 9601,
                "cases": cases,
                "frames_compared": frames_compared,
                "rms_values_compared": rms_compared,
                "max_coordinate_absolute_difference": max_coordinate_difference,
                "max_rms_absolute_difference": max_rms_difference,
                "coordinate_equality_claimed": False,
                "rms_tolerance": {"rtol": 1e-7, "atol": 1e-7},
                "rigid_distance_checks": cases,
                "input_sha256": checksum.hexdigest(),
                "degenerate_diagnostics": degenerate,
                "scope": "paired explicit coordinate arrays; proper rotation, ordinary positive weights; no atom mapping, generation or Mol mutation",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
