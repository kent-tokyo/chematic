"""Bounded NumPy adapters for the existing native RDKit matrix profiles."""

from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any, SupportsIndex, Union

import numpy as np

if TYPE_CHECKING:
    from . import Mol as NativeMol
else:
    from .chematic import Mol as NativeMol

AdjacencyMatrix = Union[
    np.ndarray[Any, np.dtype[np.int32]], np.ndarray[Any, np.dtype[np.float64]]
]
DistanceMatrix = np.ndarray[Any, np.dtype[np.float64]]
_MAX_ATOMS = 1024


def _boolean(value: bool, name: str) -> None:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be bool")


def _options(mol: NativeMol, force: bool, prefix: str) -> int:
    _boolean(force, "force")
    if not isinstance(prefix, str):
        raise TypeError("prefix must be str")
    if prefix:
        raise NotImplementedError(
            "matrix cache prefixes are unsupported; use prefix=''"
        )
    count = mol.atom_count
    if count > _MAX_ATOMS:
        raise ValueError(f"matrix resource limit: {count} input atoms exceeds 1024")
    return count


def _check_rows(rows: list[list[float]], count: int) -> None:
    if len(rows) != count:
        raise ValueError(
            "native matrix preparation changed atom indices (for example by removing "
            "explicit H); use the native rdkit matrix API for its prepared indexing"
        )
    if any(len(row) != count for row in rows):
        raise ValueError("native matrix is not square")


def adjacency_matrix(
    mol: NativeMol,
    use_bo: bool,
    empty_value: SupportsIndex,
    force: bool,
    prefix: str,
) -> AdjacencyMatrix:
    """Return an independent int32/float64 adjacency; refuse nonzero fill."""
    _boolean(use_bo, "useBO")
    if isinstance(empty_value, (bool, np.bool_)):
        raise TypeError("emptyVal must be an integer, excluding bool")
    fill = operator.index(empty_value)
    if fill != 0:
        raise ValueError("only emptyVal=0 is supported")
    count = _options(mol, force, prefix)
    rows = mol.rdkit_adjacency_matrix(use_bo=use_bo, empty_value=0)
    _check_rows(rows, count)
    output: AdjacencyMatrix = np.empty(
        (count, count), dtype=np.float64 if use_bo else np.int32
    )
    if count:
        output[:, :] = rows
    return output


def distance_matrix(
    mol: NativeMol,
    use_bo: bool,
    use_atom_weights: bool,
    force: bool,
    prefix: str,
) -> DistanceMatrix:
    """Return an independent float64 topological distance matrix."""
    _boolean(use_bo, "useBO")
    _boolean(use_atom_weights, "useAtomWts")
    count = _options(mol, force, prefix)
    rows = mol.rdkit_distance_matrix(use_bo=use_bo, use_atom_wts=use_atom_weights)
    _check_rows(rows, count)
    output: DistanceMatrix = np.empty((count, count), dtype=np.float64)
    if count:
        output[:, :] = rows
    return output
