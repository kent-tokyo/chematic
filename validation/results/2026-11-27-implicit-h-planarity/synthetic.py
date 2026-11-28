"""Seeded synthetic implicit-H carbon directions compared with RDKit 3D tags."""
import collections
import json
import platform

import chematic
import numpy as np
from rdkit import Chem, rdBase

rdBase.DisableLog("rdApp.*")
smiles = "F[C@H](Cl)Br"
mol = chematic.from_smiles(smiles)
reference = Chem.MolFromSmiles(smiles)
declared = reference.GetAtomWithIdx(1).GetChiralTag()
tetra = {Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
results = []
for label, seed, samples in [("threshold_selection", 20261128, 30000),
                             ("independent_synthetic_seed", 920261127, 50000)]:
    rng = np.random.default_rng(seed)
    counts = collections.Counter()
    for _ in range(samples):
        v = rng.normal(size=(3, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        v *= rng.uniform(0.8, 2.0, size=(3, 1))
        coords = [v[0].tolist(), [0., 0., 0.], *v[1:].tolist()]
        raw = mol.stereo_3d_diagnostics(coords)["atom_stereo"][0]["status"]
        screen = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)["atom_stereo"][0]["status"]
        ref = Chem.Mol(reference)
        conformer = Chem.Conformer(4)
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
    results.append({"label": label, "seed": seed, "samples": samples,
                    "counts": [{"oracle": k[0], "screen": k[1], "count": v}
                               for k, v in sorted(counts.items())]})

print(json.dumps({"smiles": smiles, "platform": platform.platform(),
                  "python": platform.python_version(), "numpy": np.__version__,
                  "rdkit": rdBase.rdkitVersion, "results": results}, indent=2))
