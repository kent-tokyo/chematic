"""Typed, checked single-bond fragmentation with caller-owned cut counters."""

from __future__ import annotations

from collections.abc import Sequence as RuntimeSequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, List, Optional, Sequence, Tuple

if TYPE_CHECKING:
    from . import Mol
else:
    from .chematic import Mol


@dataclass(frozen=True)
class FragmentAttachment:
    """Original bond/anchor and appended dummy atom index and isotope label."""

    source_bond: int
    anchor: int
    dummy: int
    isotope: int


@dataclass(frozen=True)
class BondFragmentation:
    """New independent Mol, immutable counts and attachment provenance."""

    molecule: Mol
    cuts_per_atom: Tuple[int, ...]
    attachments: Tuple[FragmentAttachment, ...]


def _sequence(value: object, name: str) -> None:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(
        value, RuntimeSequence
    ):
        raise TypeError(f"{name} must be a sequence")


def _integer(value: int, maximum: int, name: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer, not bool")
    if not 0 <= value <= maximum:
        raise ValueError(f"{name} must be between 0 and {maximum}")
    return value


def fragment_on_bonds_checked(
    mol: Mol,
    bond_indices: Sequence[int],
    *,
    dummy_labels: Optional[Sequence[Tuple[int, int]]] = None,
    cuts_per_atom: Optional[List[int]] = None,
) -> BondFragmentation:
    """Cut bridge single bonds and append two isotope-labelled dummy atoms each.

    Args:
        mol: Native concrete molecule, at most 4096 atoms and 8192 bonds.
        bond_indices: Distinct original bond indices, at most 1024. Each must
            be a plain single bond outside a ring.
        dummy_labels: One (begin_label, end_label) pair per cut, each 0..65535.
            The begin label's dummy attaches to the original end atom, and the
            end label's dummy to the begin atom. Defaults to their atom indices.
            Zero means an unlabelled dummy.
        cuts_per_atom: Optional plain list[int], exactly one entry per original
            atom. Entries must be integers (not bool); initial values are ignored.
            Replaced with counts only after all validation and fragmentation
            succeed. Native/type errors leave it unchanged.

    Returns:
        Frozen result with a new Mol, counts and attachment records. Original
        atom indices, atom attributes and hydrogen counts remain. Dummies are
        appended in input cut order. Python molecule properties are copied into
        an independent store. Empty cuts clone the graph and write zero counts.

    Raises:
        TypeError: Wrong argument or counter element type.
        ValueError: Invalid/duplicate indices or labels, size cap, ring/non-single
            cut, or unsupported graph. Nonempty edits reject stereo annotations,
            existing dummy/query atoms and nonstandard bond types anywhere in
            the input. A no-op with no cuts retains existing graph annotations.

    This is a bounded chematic API, not the entire RDKit FragmentOnBonds API.
    It always adds single-bond dummies, never infers new bonds or generates
    coordinates. Sequence arguments accept list/tuple/range, not lazy iterators.
    Use SMILES for fragment interchange. Existing V2000 and MolJSON writers do
    not reliably preserve isotope-labelled dummy atoms; see the checked bond
    fragmentation guide for the validated V3000 fields and export limitations.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic Mol")
    if mol.atom_count > 4096:
        raise ValueError("bond fragmentation requires at most 4096 atoms")
    _sequence(bond_indices, "bond_indices")
    if len(bond_indices) > 1024:
        raise ValueError("at most 1024 cuts are supported")
    indices = [
        _integer(bond_indices[i], 2**32 - 1, "bond index")
        for i in range(len(bond_indices))
    ]
    labels = None
    if dummy_labels is not None:
        _sequence(dummy_labels, "dummy_labels")
        if len(dummy_labels) != len(indices):
            raise ValueError("dummy_labels must have one pair per cut")
        labels = []
        for i in range(len(indices)):
            pair = dummy_labels[i]
            _sequence(pair, "dummy label pair")
            if len(pair) != 2:
                raise ValueError("each dummy label pair must contain two integers")
            labels.append(
                (
                    _integer(pair[0], 65535, "dummy label"),
                    _integer(pair[1], 65535, "dummy label"),
                )
            )
    if cuts_per_atom is not None:
        if type(cuts_per_atom) is not list:
            raise TypeError("cuts_per_atom must be a plain list[int]")
        if len(cuts_per_atom) != mol.atom_count:
            raise ValueError("cuts_per_atom must have one entry per original atom")
        if any(type(value) is not int for value in cuts_per_atom):
            raise TypeError("cuts_per_atom entries must be integers, not bool")
    molecule, counts, raw_attachments = mol._fragment_single_bonds_checked(
        indices, labels
    )
    result = BondFragmentation(
        molecule, tuple(counts), tuple(FragmentAttachment(*a) for a in raw_attachments)
    )
    if cuts_per_atom is not None:
        cuts_per_atom[:] = counts
    return result
