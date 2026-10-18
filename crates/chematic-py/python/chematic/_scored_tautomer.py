"""Bounded supplied-candidate selection with owned custom-scoring snapshots."""

from __future__ import annotations

import operator
from typing import (
    TYPE_CHECKING,
    Callable,
    Iterable,
    Protocol,
    Sequence,
    SupportsIndex,
    cast,
)
import numpy as np

if TYPE_CHECKING:
    from . import Mol
else:
    from .chematic import Mol

_MAX_CANDIDATES = 4096
_MAX_ATOMS = 1024
_MAX_BONDS = 4096
_MAX_TOTAL_ATOMS = 16_384
_MAX_PAYLOAD_BYTES = 1 << 20
_MAX_TOTAL_PAYLOAD_BYTES = 32 << 20


class _BondView(Protocol):
    @property
    def bond_table(self) -> Sequence[object]: ...


def _integer(value: SupportsIndex, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer, excluding bool")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc


def _snapshot(molecule: Mol, index: int) -> tuple[Mol, bytes, Callable[[bytes], Mol]]:
    if not isinstance(molecule, Mol):
        raise TypeError(f"candidate {index} must be a chematic Mol")
    if (
        molecule.atom_count > _MAX_ATOMS
        or len(cast(_BondView, molecule).bond_table) > _MAX_BONDS
    ):
        raise ValueError(f"candidate {index} exceeds atom/bond limits")
    # Use the native, versioned state codec directly. No pickle stream or
    # externally supplied payload is loaded, and subclasses cannot replace the
    # unbound native reducer. The codec applies its own allocation bounds too.
    restore_raw, arguments = Mol.__reduce__(molecule)
    payload = arguments[0]
    if len(payload) > _MAX_PAYLOAD_BYTES:
        raise ValueError(f"candidate {index} exceeds snapshot byte limit")
    restore = cast(Callable[[bytes], Mol], restore_raw)
    return restore(payload), payload, restore


def rdkit_pick_tautomer_scored(
    candidates: Iterable[Mol],
    score_func: Callable[[Mol], SupportsIndex],
    *,
    max_candidates: SupportsIndex = 1000,
) -> Mol:
    """Choose an owned supplied candidate using a typed integer score callback.

    Args:
        candidates: Native chematic molecules; generators are supported.
        score_func: Called once per candidate, in input order, on an isolated
            snapshot. Return a signed 32-bit integer (not bool or float).
        max_candidates: Integer in 1..4096, excluding bool.

    Returns:
        Independent copy of the highest-scoring source, with canonical-SMILES
        lexical ties and first occurrence for identical ties. No enumeration,
        tautomer-family check or final stereo cleanup is performed.

    Raises:
        TypeError: Wrong candidate, callback, limit or score type.
        ValueError: Empty input, graph/snapshot/collection limits, unsupported
            canonicalization or score outside signed 32-bit range.

    All inputs and canonical keys are prepared before callbacks start. Each
    yielded molecule is snapshotted before advancing the source iterator.
    Callback edits affect neither retained candidates nor the returned source
    copy. Source/callback exceptions propagate; external side effects cannot
    be rolled back. Limits bound data, not caller execution time. See the
    scored-tautomer guide for byte/atom limits and RDKit tie-policy differences.
    """
    limit = _integer(max_candidates, "max_candidates")
    if not 1 <= limit <= _MAX_CANDIDATES:
        raise ValueError("max_candidates must be in 1..4096")
    if not callable(score_func):
        raise TypeError("score_func must be callable")
    prepared: list[tuple[Mol, bytes, Callable[[bytes], Mol], str]] = []
    atoms = payload_bytes = 0
    for index, molecule in enumerate(candidates):
        if index >= limit:
            raise ValueError("candidates exceed max_candidates")
        snapshot, payload, restore = _snapshot(molecule, index)
        atoms += snapshot.atom_count
        payload_bytes += len(payload)
        if atoms > _MAX_TOTAL_ATOMS or payload_bytes > _MAX_TOTAL_PAYLOAD_BYTES:
            raise ValueError("candidates exceed combined snapshot/atom limits")
        try:
            key = snapshot.rdkit_smiles
        except ValueError as exc:
            raise ValueError(
                f"candidate {index} cannot be canonicalized: {exc}"
            ) from exc
        prepared.append((snapshot, payload, restore, key))
    if not prepared:
        raise ValueError("at least one candidate is required")
    winner = 0
    best: tuple[int, str] | None = None
    for index, (_, payload, restore, key) in enumerate(prepared):
        score = _integer(score_func(restore(payload)), f"score for candidate {index}")
        if not -(1 << 31) <= score < (1 << 31):
            raise ValueError(
                f"score for candidate {index} must fit signed 32-bit range"
            )
        rank = (-score, key)
        if best is None or rank < best:
            best, winner = rank, index
    return prepared[winner][0]
