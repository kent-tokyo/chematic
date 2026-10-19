from typing import Any
import numpy as np
from chematic import Mol, FragmentCutCounts, fragment_on_bonds_checked


def consume(mol: Mol) -> None:
    counts: FragmentCutCounts = [0] * mol.atom_count
    fragment_on_bonds_checked(mol, [0], cuts_per_atom=counts)
    array: np.ndarray[Any, np.dtype[np.uint32]] = np.zeros(
        mol.atom_count, dtype=np.uint32
    )
    fragment_on_bonds_checked(mol, [0], cuts_per_atom=array[::-1])
