"""Public typed matrix usage, without depending on installed RDKit."""

from typing import Any
import numpy as np
from chematic import rdkit_compat as C

mol = C.MolFromSmiles("CCO")
assert mol is not None
adjacency = C.GetAdjacencyMatrix(mol, useBO=True, emptyVal=np.int64(0), force=True)
distance: np.ndarray[Any, np.dtype[np.float64]] = C.GetDistanceMatrix(
    mol, useAtomWts=True
)
