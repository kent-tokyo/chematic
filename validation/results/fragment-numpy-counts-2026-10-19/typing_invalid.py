from typing import Any
import numpy as np
from chematic import Mol, fragment_on_bonds_checked


def consume(mol: Mol, floats: np.ndarray[Any, np.dtype[np.float64]]) -> None:
    fragment_on_bonds_checked(mol, [0], cuts_per_atom=floats)
    fragment_on_bonds_checked(mol, [0], cuts_per_atom=(0, 0))
