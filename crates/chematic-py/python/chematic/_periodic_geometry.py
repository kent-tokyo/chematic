"""Checked 3D TURBOMOLE/POSCAR conversion using native geometry parsers."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    from . import PeriodicStructure, parse_turbomole, write_turbomole
else:
    from .chematic import PeriodicStructure, parse_turbomole, write_turbomole

_MAX_BYTES = 16 << 20
_MAX_ATOMS = 50_000
_MAX_LINE = 1024


@dataclass(frozen=True)
class PeriodicGeometryConversion:
    """Output geometry and complete source, with explicit atom correspondence.

    source_indices[i] is the zero-based source coordinate row for output atom i.
    POSCAR groups species in first-appearance order, stably within each species.
    Input text, comments and original units remain available in source_text;
    output uses canonical units/formatting. No files are written.
    """

    source_text: str
    text: str
    source_format: Literal["turbomole", "poscar"]
    target_format: Literal["turbomole", "poscar"]
    source_indices: tuple[int, ...]


def _text(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("text must be str")
    if len(value) > _MAX_BYTES:
        raise ValueError("periodic geometry input exceeds 16 MiB")
    try:
        if len(value.encode("utf-8")) > _MAX_BYTES:
            raise ValueError("periodic geometry input exceeds 16 MiB")
    except UnicodeError as exc:
        raise ValueError("text must be valid UTF-8") from exc
    if any(len(line.encode("utf-8")) > _MAX_LINE for line in value.splitlines()):
        raise ValueError("periodic geometry line exceeds 1024 bytes")
    return value


def _groups(elements: list[str]) -> tuple[list[str], list[list[int]]]:
    groups: dict[str, list[int]] = {}
    for index, element in enumerate(elements):
        groups.setdefault(element, []).append(index)
    return list(groups), list(groups.values())


def turbomole_to_poscar_checked(text: str) -> PeriodicGeometryConversion:
    """Convert a standalone periodic-3 TURBOMOLE geometry to VASP 5 POSCAR.

    Args:
        text: At most 16 MiB, 50000 atoms and 1024 bytes per line. Native
            TURBOMOLE unit, lattice/cell and coordinate checks apply.

    Returns:
        A frozen report containing original text, Direct-coordinate POSCAR and
        output-to-source atom indices. Lattice rows retain their orientation;
        positions are not wrapped. Frozen atoms become F F F, others T T T.

    Raises:
        TypeError: Input is not str.
        ValueError: Invalid/unsupported input, missing 3D lattice, singular cell,
            unknown groups, unsupported periodicity or resource limit.

    No cell, topology, charge or chemical bonding is inferred. Native parsers
    perform the input interpretation and Cartesian-to-fractional conversion.
    """
    source = _text(text)
    geometry: dict[str, Any] = parse_turbomole(
        source,
        max_input_bytes=_MAX_BYTES,
        max_atoms=_MAX_ATOMS,
        max_line_len=_MAX_LINE,
    )
    if geometry["periodicity"] != 3 or geometry["lattice_angstrom"] is None:
        raise ValueError("POSCAR conversion requires an explicit 3D lattice")
    species, groups = _groups(geometry["elements"])
    order = tuple(index for group in groups for index in group)
    rows = ["chematic TURBOMOLE conversion", "1.0"]
    rows += [" ".join(map(str, row)) for row in geometry["lattice_angstrom"]]
    rows += [
        " ".join(species),
        " ".join(str(len(group)) for group in groups),
        "Selective dynamics",
        "Cartesian",
    ]
    for index in order:
        flags = "F F F" if geometry["frozen"][index] else "T T T"
        rows.append(
            " ".join(map(str, geometry["positions_angstrom"][index])) + " " + flags
        )
    intermediate = "\n".join(rows) + "\n"
    output = PeriodicStructure.from_poscar(intermediate).to_poscar()
    restored = PeriodicStructure.from_poscar(output)
    if any(
        not math.isfinite(float(v))
        for row in restored.cartesian_positions()
        for v in row
    ):
        raise ValueError("POSCAR coordinate conversion overflow")
    return PeriodicGeometryConversion(source, output, "turbomole", "poscar", order)


def poscar_to_turbomole_checked(text: str) -> PeriodicGeometryConversion:
    """Convert VASP 5 geometry to standalone periodic-3 TURBOMOLE in bohr.

    Args:
        text: At most 16 MiB, 50000 atoms and 1024 bytes per line. Explicit
            species names are required. Native POSCAR scale/coordinate checks
            apply; source atom order is reported through source_indices.

    Returns:
        Original text and TURBOMOLE with lattice rows, Cartesian positions and
        fully fixed atoms retained. Canonical species grouping is explicit.

    Raises:
        TypeError: Input is not str.
        ValueError: Invalid/unsupported geometry, resource limit, velocities,
            predictor-corrector data or partially fixed atom directions.

    Partial directional constraints and restart data cannot be represented by
    this TURBOMOLE profile and are refused instead of silently discarded.
    A comment/title is kept in source_text, not interpreted as a data group.
    """
    source = _text(text)
    lines = source.splitlines()
    if len(lines) < 7:
        raise ValueError("VASP 5 input requires species and atom-count header lines")
    # Bound native allocations before from_poscar expands the count header.
    fields = lines[6].split()
    if not fields or len(fields) != len(lines[5].split()):
        raise ValueError("POSCAR species/count lengths must match")
    counts = []
    for value in fields:
        if not value.isascii() or not value.isdecimal():
            raise ValueError("POSCAR counts must be nonnegative decimal integers")
        count = (
            int(value.lstrip("0") or "0")
            if len(value.lstrip("0")) <= 5
            else _MAX_ATOMS + 1
        )
        if count > _MAX_ATOMS:
            raise ValueError("periodic conversion supports at most 50000 atoms")
        counts.append(count)
    if not 1 <= sum(counts) <= _MAX_ATOMS:
        raise ValueError("periodic conversion requires 1..50000 atoms")
    structure = PeriodicStructure.from_poscar(source)
    elements = [site.species[0][0] for site in structure.sites]
    _, groups = _groups(elements)
    order = tuple(index for group in groups for index in group)
    canonical = structure.to_poscar().splitlines()
    selective = canonical[7].startswith("Selective")
    start = 9 if selective else 8
    if any(line.strip() for line in canonical[start + len(order) :]):
        raise ValueError("TURBOMOLE conversion cannot retain POSCAR restart data")
    frozen = []
    for row in canonical[start : start + len(order)]:
        flags = row.split()[3:] if selective else ["T", "T", "T"]
        if flags not in (["T", "T", "T"], ["F", "F", "F"]):
            raise ValueError("TURBOMOLE requires fully free or fully fixed atoms")
        frozen.append(flags == ["F", "F", "F"])
    positions = structure.cartesian_positions()
    output = write_turbomole(
        [elements[i] for i in order],
        [[float(value) for value in positions[i]] for i in order],
        lattice_angstrom=structure.lattice.matrix.tolist(),
        frozen=frozen,
    )
    return PeriodicGeometryConversion(source, output, "poscar", "turbomole", order)
