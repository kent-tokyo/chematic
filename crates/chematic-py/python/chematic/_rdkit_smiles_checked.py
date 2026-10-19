"""Bounded stability checks for the legacy RDKit-compatible writer."""

from __future__ import annotations

from collections.abc import Iterator
import itertools
from typing import TYPE_CHECKING, Literal, Optional

if TYPE_CHECKING:
    from . import Mol, from_smiles
else:
    from .chematic import Mol, from_smiles


class RdkitSmilesStabilityError(ValueError):
    """A candidate failed a fixed-point or atom-order stability check.

    The candidate strings are diagnostic evidence, not accepted identifiers.
    ``atom_order`` is populated only for permutation-related failures.
    """

    def __init__(
        self,
        code: Literal[
            "unreadable_output",
            "rewrite_failed",
            "non_idempotent",
            "permutation_write_failed",
            "permutation_dependent",
        ],
        first_smiles: str,
        roundtrip_smiles: Optional[str],
        atom_order: Optional[tuple[int, ...]] = None,
    ) -> None:
        self.code = code
        self.first_smiles = first_smiles
        self.roundtrip_smiles = roundtrip_smiles
        self.atom_order = atom_order
        super().__init__(
            f"RDKit-compatible canonical SMILES failed stability check: {code}"
        )


def _bounded_atom_orders(atom_count: int, limit: int) -> Iterator[tuple[int, ...]]:
    """Yield broad, deterministic new-to-old orders without factorial growth."""
    identity = tuple(range(atom_count))
    yielded = {identity}
    yield identity
    if atom_count < 2 or limit == 1:
        return

    def emit(order: tuple[int, ...]) -> Iterator[tuple[int, ...]]:
        if order not in yielded and len(yielded) < limit:
            yielded.add(order)
            yield order

    yield from emit(tuple(reversed(identity)))
    for shift in range(1, atom_count):
        yield from emit(identity[shift:] + identity[:shift])
        if len(yielded) >= limit:
            return
    yield from emit(identity[::2] + identity[1::2])
    yield from emit(identity[1::2] + identity[::2])

    if atom_count <= 8:
        for order in itertools.permutations(identity):
            yield from emit(order)
            if len(yielded) >= limit:
                return
        return

    # Local xorshift32 keeps this sequence stable across Python versions and
    # does not mutate the process-wide random state.
    state = (0x9E3779B9 ^ atom_count) & 0xFFFFFFFF
    attempts = 0
    while len(yielded) < limit and attempts < limit * 32:
        attempts += 1
        order = list(identity)
        for index in range(atom_count - 1, 0, -1):
            state ^= (state << 13) & 0xFFFFFFFF
            state ^= state >> 17
            state ^= (state << 5) & 0xFFFFFFFF
            state &= 0xFFFFFFFF
            swap = state % (index + 1)
            order[index], order[swap] = order[swap], order[index]
        yield from emit(tuple(order))


def _write_permuted(mol: Mol, order: tuple[int, ...]) -> str:
    """Write one atom-renumbered copy; separated for fault-injection tests."""
    return mol.renumber_atoms(order).rdkit_smiles


def rdkit_smiles_checked(mol: Mol, /, *, max_permutations: int = 32) -> str:
    """Return RDKit-compatible SMILES only after bounded stability checks.

    Write ``mol.rdkit_smiles``, parse that string, then write it again. Reject
    unequal strings with :class:`RdkitSmilesStabilityError`; never select a
    member of a cycle, remove stereo tags or retry until convergence. The
    Then write deterministic atom-renumbered copies and require the identical
    string. The input molecule is not modified. Initial writer failures
    propagate. ``max_permutations`` counts the identity order and is bounded
    to 2--256 for graphs with at least two atoms.

    This checks textual fixed-point and bounded atom-order invariance. It does
    not prove preservation of chemical identity, complete stereochemistry or
    uniqueness over every possible graph representation. RDKit issue #9391 is
    one known legacy writer failure caught by the fixed-point stage.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic.Mol")
    if isinstance(max_permutations, bool) or not isinstance(max_permutations, int):
        raise TypeError("max_permutations must be an integer")
    if not 2 <= max_permutations <= 256:
        raise ValueError("max_permutations must be between 2 and 256")
    first = mol.rdkit_smiles
    try:
        readback = from_smiles(first)
    except (ValueError, RuntimeError) as error:
        raise RdkitSmilesStabilityError("unreadable_output", first, None) from error
    try:
        rewritten = readback.rdkit_smiles
    except (ValueError, RuntimeError) as error:
        raise RdkitSmilesStabilityError("rewrite_failed", first, None) from error
    if first != rewritten:
        raise RdkitSmilesStabilityError("non_idempotent", first, rewritten)
    for order in _bounded_atom_orders(mol.atom_count, max_permutations):
        if order == tuple(range(mol.atom_count)):
            continue
        try:
            candidate = _write_permuted(mol, order)
        except (ValueError, RuntimeError) as error:
            raise RdkitSmilesStabilityError(
                "permutation_write_failed", first, None, order
            ) from error
        if candidate != first:
            raise RdkitSmilesStabilityError(
                "permutation_dependent", first, candidate, order
            )
    return first
