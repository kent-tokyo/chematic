from typing import Any
import numpy as np
from chematic import rdkit_compat as C


class Mask:
    def GetNumBits(self) -> int:
        return 2048

    def GetOnBits(self) -> tuple[int, ...]:
        return (0, 2)


def consume(mol: C.Mol) -> C.ExplicitBitVect:
    counts: C.PatternAtomCounts = [0] * mol.GetNumAtoms()
    C.PatternFingerprint(mol, atomCounts=counts, setOnlyBits=Mask())
    array: np.ndarray[Any, np.dtype[np.uint32]] = np.zeros(
        mol.GetNumAtoms(), dtype=np.uint32
    )
    return C.PatternFingerprint(mol, fpSize=128, atomCounts=array)
