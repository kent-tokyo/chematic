"""Checked concrete core removal with explicit, typed atom correspondence."""
from __future__ import annotations

from collections.abc import Sequence as RuntimeSequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional, Sequence, Tuple

if TYPE_CHECKING:
    from . import Mol
else:
    from .chematic import Mol


@dataclass(frozen=True)
class CoreAttachment:
    """Core/query index, original endpoints/bond and new dummy/anchor indices."""

    core_atom: int
    source_core_atom: int
    source_anchor: int
    source_bond: int
    anchor: int
    dummy: int


@dataclass(frozen=True)
class CoreReplacement:
    """New independent Mol, original-atom remapping and attachment provenance."""

    molecule: Mol
    original_to_output: Tuple[Optional[int], ...]
    attachments: Tuple[CoreAttachment, ...]


def replace_core_checked(mol: Mol, core: Mol, match: Sequence[int]) -> CoreReplacement:
    """Remove an explicitly mapped concrete core and label its side chains.

    Args:
        mol: Concrete source graph, at most 4096 atoms and 8192 bonds.
        core: Concrete nonempty graph, at most 1024 atoms and 8192 bonds.
        match: Eager sequence of distinct source atom indices. Entry i maps
            core atom i to the source. Integers exclude bool. No search is run.

    Returns:
        Frozen result with a new Mol, original-to-output indices (None for
        removed atoms), and one attachment per cut in source bond order.
        Retained source atoms keep order, attributes, H counts and caller tags.
        Dummies follow them, labelled by the core atom index (zero unlabelled).
        Python properties are copied independently. Unrelated components stay;
        removing the entire source returns an empty molecule.

    Raises:
        TypeError: Incorrect molecule, sequence or index type.
        ValueError: Bounds, duplicate/out-of-range indices, mismatched atoms or
            core bonds, or an unsupported graph/boundary. Atoms must agree in
            element, aromaticity, charge and isotope; H counts and maps are
            ignored in matching. Only bridge single bonds may cross the core
            boundary, at most 1024 cuts. Both graphs must be stereo-free with
            no pre-existing dummy/query atoms or nonstandard bond orders.

    A bounded native counterpart of RDKit ReplaceCore with labelByIndex=True.
    Automatic matching, query cores, dummy toggles and stereo transfer are
    outside this profile. No aromaticity perception is run; Kekule source and
    aromatic core representations can be refused. Neither input nor the
    supplied match is modified.
    Use SMILES for fragment interchange: existing V2000/MolJSON writers have
    known isotope-dummy preservation limits (see checked bond fragmentation).
    """
    if not isinstance(mol, Mol) or not isinstance(core, Mol):
        raise TypeError("mol and core must be chematic Mol objects")
    if mol.atom_count > 4096 or not 1 <= core.atom_count <= 1024:
        raise ValueError("requires at most 4096 source atoms and 1..1024 core atoms")
    if isinstance(match, (str, bytes, bytearray)) or not isinstance(match, RuntimeSequence):
        raise TypeError("match must be an eager sequence of integer source indices")
    if len(match) != core.atom_count:
        raise ValueError("match must contain one source atom per core atom")
    indices = []
    for i in range(len(match)):
        index = match[i]
        if type(index) is not int:
            raise TypeError("match indices must be integers, not bool")
        if not 0 <= index <= 2**32 - 1:
            raise ValueError("match index must fit an unsigned 32-bit integer")
        indices.append(index)
    output, remap, attachments = mol._replace_core_checked(core, indices)
    return CoreReplacement(output, tuple(remap), tuple(CoreAttachment(*a) for a in attachments))
