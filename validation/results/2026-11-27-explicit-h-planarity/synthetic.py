"""Explicit H carbon with reconstructed and independent synthetic H directions."""
import collections
import hashlib
import json
import platform

import chematic
import numpy as np
from rdkit import Chem, rdBase

rdBase.DisableLog("rdApp.*")
smiles = "F[C@]([H])(Cl)Br"
mol = chematic.from_smiles(smiles)
parser = Chem.SmilesParserParams()
parser.removeHs = False
reference = Chem.MolFromSmiles(smiles, parser)
assert reference.GetNumAtoms() == 5 and mol.atom_count == 5
declared = reference.GetAtomWithIdx(1).GetChiralTag()
tetra = {Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
results = []
for label, seed, samples in [("reconstructed_h", 20261129, 30000),
                             ("independent_h", 20261130, 30000)]:
    rng = np.random.default_rng(seed)
    counts = collections.Counter()
    status_digest = hashlib.sha256()
    for _ in range(samples):
        directions = rng.normal(size=(3, 3))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        heavy = directions * rng.uniform(0.8, 2.0, size=(3, 1))
        if label == "reconstructed_h":
            h = -directions.sum(axis=0)
            length = np.linalg.norm(h)
            if length < 1e-10:
                counts[("unusable_phantom", "not_checked")] += 1
                status_digest.update(b"unusable_phantom:not_checked\n")
                continue
            h /= length
        else:
            h = rng.normal(size=3)
            h /= np.linalg.norm(h)
            h *= rng.uniform(0.8, 2.0)
        coords = [heavy[0].tolist(), [0., 0., 0.], h.tolist(),
                  heavy[1].tolist(), heavy[2].tolist()]
        raw = mol.stereo_3d_diagnostics(coords)["atom_stereo"][0]["status"]
        screen = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)["atom_stereo"][0]["status"]
        ref = Chem.Mol(reference)
        conformer = Chem.Conformer(5)
        for idx, point in enumerate(coords):
            conformer.SetAtomPosition(idx, point)
        ref.AddConformer(conformer)
        Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
        observed = ref.GetAtomWithIdx(1).GetChiralTag()
        if observed not in tetra:
            classification = "rdkit_unassigned"
        elif raw not in ("satisfied", "violated"):
            classification = "raw_unevaluable"
        elif (raw == "satisfied") == (observed == declared):
            classification = "raw_agrees"
        else:
            classification = "raw_disagrees"
        counts[(classification, screen)] += 1
        status_digest.update(f"{classification}:{screen}\n".encode())
    results.append({"label": label, "seed": seed, "samples": samples,
                    "status_sequence_sha256": status_digest.hexdigest(),
                    "counts": [{"oracle": k[0], "screen": k[1], "count": v}
                               for k, v in sorted(counts.items())]})

print(json.dumps({"smiles": smiles, "platform": platform.platform(),
                  "python": platform.python_version(), "numpy": np.__version__,
                  "rdkit": rdBase.rdkitVersion, "results": results}, indent=2))
