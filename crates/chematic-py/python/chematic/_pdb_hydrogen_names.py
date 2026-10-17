"""Explicit, source-preserving PDB hydrogen names (Open Babel #2945)."""

from __future__ import annotations

from dataclasses import dataclass
import operator
from typing import Dict, List, Optional, Sequence, Set, SupportsIndex, Tuple

import numpy as np

from .chematic import select_pdb_atoms

_MAX_BYTES = 16 << 20
_MAX_ATOMS = 99_999
_MAX_SUFFIX = 999
_Group = Tuple[int, str, int, str, str]


@dataclass(frozen=True)
class PdbHydrogenName:
    """One selected H row, with original index/serial and raw four-column names."""

    source_index: int
    serial: int
    old_name: str
    new_name: str


@dataclass(frozen=True)
class PdbHydrogenNaming:
    """Complete in-memory result; no coordinates, atom types or files are changed."""

    source_text: str
    pdb: str
    atom_count: int
    assignments: Tuple[PdbHydrogenName, ...]
    renamed_count: int
    updated_anisou_lines: Tuple[int, ...]


@dataclass(frozen=True)
class _Atom:
    index: int
    serial: int
    group: _Group
    name: str
    element: str
    alt_loc: str


def _atoms(lines: Sequence[str]) -> List[_Atom]:
    result: List[_Atom] = []
    segment = 0
    for number, line in enumerate(lines, 1):
        record = line[:6].strip()
        if record == "TER":
            segment += 1
        elif record in ("ATOM", "HETATM"):
            element = line[76:78].strip()
            if not element:
                raise ValueError(f"PDB naming line {number}: explicit element required")
            group = (segment, line[21], int(line[22:26]), line[26], line[17:20].strip())
            result.append(_Atom(
                len(result), int(line[6:11]), group, line[12:16], element, line[16]
            ))
    return result


def _selection(
    atoms: Sequence[_Atom], values: Optional[Sequence[SupportsIndex]]
) -> Set[int]:
    if values is None:
        return {atom.index for atom in atoms if atom.element == "H"}
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise TypeError("hydrogen_indices must be an indexed integer sequence or None")
    if len(values) > len(atoms):
        raise ValueError("hydrogen_indices exceeds the atom count")
    result: Set[int] = set()
    for offset in range(len(values)):
        value = values[offset]
        if isinstance(value, (bool, np.bool_)):
            raise TypeError("hydrogen index must be an integer, excluding bool")
        index = operator.index(value)
        if not 0 <= index < len(atoms) or index in result:
            raise ValueError("hydrogen indices must be unique and in range")
        if atoms[index].element != "H":
            raise ValueError("selected atom must have explicit element H")
        result.add(index)
    return result


def _assign(
    atoms: Sequence[_Atom], selected: Set[int]
) -> Tuple[PdbHydrogenName, ...]:
    used: Dict[_Group, Set[str]] = {}
    next_suffix: Dict[_Group, int] = {}
    for atom in atoms:
        names = used.setdefault(atom.group, set())
        if atom.index not in selected:
            names.add(atom.name.strip())
    assignments = []
    for atom in atoms:
        if atom.index not in selected:
            continue
        if atom.alt_loc != " ":
            raise ValueError("selected hydrogen alternate locations are unsupported")
        names = used[atom.group]
        suffix = next_suffix.get(atom.group, 1)
        while suffix <= _MAX_SUFFIX and f"H{suffix}" in names:
            suffix += 1
        if suffix > _MAX_SUFFIX:
            raise ValueError("PDB hydrogen name space H1..H999 exhausted in residue")
        name = f"H{suffix}"
        names.add(name)
        next_suffix[atom.group] = suffix + 1
        field = (" " + name).ljust(4) if len(name) < 4 else name
        assignments.append(PdbHydrogenName(atom.index, atom.serial, atom.name, field))
    return tuple(assignments)


def name_pdb_hydrogens_checked(
    text: str, hydrogen_indices: Optional[Sequence[SupportsIndex]] = None
) -> PdbHydrogenNaming:
    """Assign H1..H999 within residues, preserving all other source bytes.

    Args:
        text: A checked single-model PDB coordinate block, at most 16 MiB and
            99999 atoms. Every atom requires a nonempty explicit element column.
        hydrogen_indices: Zero-based ATOM/HETATM row indices. None selects every
            explicit H; an empty sequence changes nothing. No element inference.

    Returns:
        Frozen report with original/output text and selected H name assignments
        in source order, including unchanged names. ANISOU names follow their H.

    Raises:
        TypeError: Wrong text/selection type or noninteger/bool index.
        ValueError: Invalid coordinate records/indices, non-H selection, selected
            H altloc, missing element, exhausted names, unsupported metadata or
            resource limits. No partial result or file writes occur.

    Residues are keyed by TER segment, chain, residue number, insertion code and
    residue name. Assigned names avoid all unselected names in the same residue.
    Only columns 13..16 of selected ATOM/HETATM and associated ANISOU change.
    Line endings, serials, coordinates, connectivity and ATOM/HETATM kinds stay
    exact. This does not add H, infer parents or impose residue-template names.
    D/T element records are unselected. Chemical naming is not inferred.
    Unknown metadata (including LINK/SSBOND/SIGATM) refuses; remove it only via
    a separate explicitly loss-consenting selection. Multiple models refuse.
    """
    if not isinstance(text, str):
        raise TypeError("text must be str")
    if len(text) > _MAX_BYTES:
        raise ValueError("PDB naming input exceeds 16 MiB")
    # Validate complete source records and links before parsing or editing names.
    select_pdb_atoms(text, max_input_bytes=_MAX_BYTES, max_atoms=_MAX_ATOMS)
    lines = text.splitlines(keepends=True)
    atoms = _atoms(lines)
    if not atoms:
        raise ValueError("PDB naming requires at least one ATOM/HETATM record")
    assignments = _assign(atoms, _selection(atoms, hydrogen_indices))
    names = {item.serial: item.new_name for item in assignments}
    output = []
    updated = []
    for number, line in enumerate(lines, 1):
        record = line[:6].strip()
        if record in ("ATOM", "HETATM", "ANISOU"):
            name = names.get(int(line[6:11]))
            if name is not None and line[12:16] != name:
                line = line[:12] + name + line[16:]
                if record == "ANISOU":
                    updated.append(number)
        output.append(line)
    pdb = "".join(output)
    select_pdb_atoms(pdb, max_input_bytes=_MAX_BYTES, max_atoms=_MAX_ATOMS)
    return PdbHydrogenNaming(
        text, pdb, len(atoms), assignments,
        sum(item.old_name != item.new_name for item in assignments), tuple(updated)
    )
