"""Typed, cancellable execution of the Rust RDKit-profile tautomer enumerator."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Literal, Optional, Tuple, cast

if TYPE_CHECKING:
    from . import Mol
else:
    from .chematic import Mol

TautomerStatus = Literal[
    "Completed", "MaxTautomersReached", "MaxTransformsReached", "Canceled"
]


@dataclass(frozen=True)
class TautomerProgress:
    """Checkpoint before a transform match; keys are current enumeration keys."""

    transforms: int
    keys: Tuple[str, ...]


@dataclass(frozen=True)
class TautomerEnumeration:
    """Final or normalized partial result, in enumeration-key order."""

    smiles: Tuple[str, ...]
    keys: Tuple[str, ...]
    status: TautomerStatus


def _limit(value: int, name: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")
    if not 1 <= value <= 10000:
        raise ValueError(f"{name} must be between 1 and 10000")
    return value


def enumerate_tautomers_checked(
    mol: Mol,
    *,
    max_tautomers: int = 1000,
    max_transforms: int = 1000,
    callback: Optional[Callable[[TautomerProgress], bool]] = None,
) -> TautomerEnumeration:
    """Enumerate in Rust; return False at a checkpoint to cancel immediately.

    Args:
        mol: Concrete chematic molecule, at most 1024 graph atoms.
        max_tautomers: Positive limit, at most 10000; default 1000.
        max_transforms: Positive limit, at most 10000; default 1000.
        callback: Callable receiving an immutable snapshot before each matched
            transform. Must return an actual bool. False cancels further
            transforms; its exception propagates unchanged. Called synchronously
            on this thread. No call is made when no transform matches.

    Returns:
        Immutable SMILES, enumeration keys, and completion/limit/cancel status.
        Cancellation still runs final stereo normalization on the partial set.
        Limit checks precede callback invocation. Current keys may differ from
        final canonical SMILES after stereo normalization.

    Raises:
        TypeError: Wrong molecule, limit, callback, or callback return type.
        ValueError: Limit/atom cap exceeded or unsupported native chemistry.

    This is a chematic API, not RDKit's SetCallback signature. A transform
    checkpoint is not a wall-clock timeout or preemption of a SMARTS match.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic Mol")
    max_tautomers = _limit(max_tautomers, "max_tautomers")
    max_transforms = _limit(max_transforms, "max_transforms")
    if callback is not None and not callable(callback):
        raise TypeError("callback must be callable or None")
    if mol.atom_count > 1024:
        raise ValueError("mol must contain at most 1024 graph atoms")

    def checkpoint(transforms: int, keys: Tuple[str, ...]) -> bool:
        if callback is None:
            return True
        return callback(TautomerProgress(transforms, keys))

    smiles, keys, status = mol._rdkit_tautomers_controlled(
        max_tautomers, max_transforms, checkpoint
    )
    return TautomerEnumeration(tuple(smiles), tuple(keys), cast(TautomerStatus, status))
