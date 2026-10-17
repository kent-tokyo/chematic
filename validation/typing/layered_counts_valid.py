from typing import Any
import numpy as np
from chematic import rdkit_compat as Chem


def calculate(mol: Chem.Mol) -> Chem.ExplicitBitVect:
    counts: list[int] = [0] * mol.GetNumAtoms()
    mask = Chem.ExplicitBitVect(2048)
    mask.SetBit(0)
    fp = Chem.LayeredFingerprint(mol, atomCounts=counts, setOnlyBits=mask)
    output: np.ndarray[Any, np.dtype[np.uint32]] = np.zeros(
        mol.GetNumAtoms(), dtype=np.uint32
    )
    Chem.LayeredFingerprint(
        mol, atomCounts=output, fromAtoms=[0], fpSize=np.int64(2048)
    )
    return fp
