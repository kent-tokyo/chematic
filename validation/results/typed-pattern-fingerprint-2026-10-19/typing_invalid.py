from typing import Any
import numpy as np
from chematic import rdkit_compat as C


def consume(mol: C.Mol, floats: np.ndarray[Any, np.dtype[np.float64]]) -> None:
    C.PatternFingerprint(mol, atomCounts=floats)
    C.PatternFingerprint(mol, fpSize=1.5)
    C.PatternFingerprint(mol, setOnlyBits=object())
    C.PatternFingerprint(mol, tautomerFingerprints="false")
