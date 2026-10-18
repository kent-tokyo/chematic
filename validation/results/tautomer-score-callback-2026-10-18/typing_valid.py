from typing import SupportsIndex
import numpy as np
import chematic as C

def score(mol: C.Mol) -> SupportsIndex:
    return np.int64(-mol.rdkit_tautomer_score())

selected: C.Mol = C.rdkit_pick_tautomer_scored(
    (C.from_smiles(s) for s in ["CC(C)=O", "C=C(C)O"]), score,
    max_candidates=np.int64(2),
)
