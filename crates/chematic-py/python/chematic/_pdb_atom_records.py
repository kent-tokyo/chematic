"""Typed single-model PDB attributes without graph or conformer inference."""

from __future__ import annotations

from dataclasses import dataclass
import operator
from typing import Optional, SupportsIndex, TypedDict, cast

from .chematic import select_pdb_atoms

_MAX_BYTES = 8 * 1024 * 1024
_MAX_ATOMS = 50_000
_MAX_RECORDS = 250_000
_MAX_LINE = 1024


@dataclass(frozen=True)
class PdbAtomRecord:
    """One original ATOM/HETATM row, including optional occupancy/B factor.

    atom_name and character identifiers retain fixed-column spaces. element_field
    and charge_field are source slices, not inferred element/charge assignments.
    record_text excludes the line ending; source_text on the read result is exact.
    """

    source_index: int
    line_number: int
    serial: int
    record_type: str
    atom_name: str
    alt_loc: str
    residue_name: str
    chain_id: str
    residue_number: int
    insertion_code: str
    position_angstrom: tuple[float, float, float]
    occupancy: Optional[float]
    b_factor: Optional[float]
    element_field: str
    charge_field: str
    record_text: str


@dataclass(frozen=True)
class PdbAtomRead:
    """Immutable atom records, model number and complete unchanged source.

    unhandled_records contains (one-based line, record name) entries reported
    as unhandled metadata by the native coordinate selection validator. All
    source records, including ANISOU/CONECT and other metadata, stay in source_text.
    """

    source_text: str
    model_num: Optional[int]
    atoms: tuple[PdbAtomRecord, ...]
    unhandled_records: tuple[tuple[int, str], ...]


class _Unhandled(TypedDict):
    line: int
    record_name: str


class _Validation(TypedDict):
    source_atom_count: int
    model_num: Optional[int]
    losses: list[_Unhandled]


def _limit(value: SupportsIndex, name: str, maximum: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} requires an integer, excluding bool")
    number = operator.index(value)
    if not 1 <= number <= maximum:
        raise ValueError(f"{name} must be in 1..{maximum}")
    return number


def _optional_number(field: str) -> Optional[float]:
    # Native validation already checked complete, finite nonblank fields.
    return float(field) if field.strip() else None


def read_pdb_atom_records(
    text: str,
    *,
    max_input_bytes: SupportsIndex = _MAX_BYTES,
    max_atoms: SupportsIndex = _MAX_ATOMS,
    max_records: SupportsIndex = _MAX_RECORDS,
    max_line_len: SupportsIndex = _MAX_LINE,
) -> PdbAtomRead:
    """Read typed atom attributes in the native single-model PDB profile.

    Args:
        text: ASCII PDB contents, not a path. Original line endings are retained.
        max_input_bytes: Input byte limit, 1..8388608.
        max_atoms: ATOM/HETATM count limit, 1..50000.
        max_records: Physical source line limit, 1..250000.
        max_line_len: Per-line byte limit, 1..1024.

    Returns:
        Immutable source-ordered records and original text. Blank/absent
        occupancy and B-factor fields are None; explicit zero stays zero.
        Coordinates are Cartesian angstroms. Character fields keep spaces.

    Raises:
        TypeError: Non-string contents or noninteger limits (including bool).
        ValueError: Non-ASCII input, exceeded limits or malformed/native-
            unsupported framing, atom fields or connectivity references.

    Uses the native checked coordinate-selection reader for validation with an
    empty selection. No selection output is exported and no atom is discarded
    from this read view. Additional records remain in source_text. No bonds,
    residue templates, altloc selection, chemical types or conformers are inferred.
    Multiple MODEL blocks and extended/hybrid-36 serials are outside this profile.
    """
    if not isinstance(text, str):
        raise TypeError("text must be PDB str contents")
    byte_limit = _limit(max_input_bytes, "max_input_bytes", _MAX_BYTES)
    atom_limit = _limit(max_atoms, "max_atoms", _MAX_ATOMS)
    record_limit = _limit(max_records, "max_records", _MAX_RECORDS)
    line_limit = _limit(max_line_len, "max_line_len", _MAX_LINE)
    if len(text) > byte_limit:
        raise ValueError("PDB atom-record input byte limit exceeded")
    if not text.isascii():
        raise ValueError("PDB atom records require ASCII input")
    validation = cast(
        _Validation,
        select_pdb_atoms(
            text,
            [],
            allow_metadata_loss=True,
            max_input_bytes=byte_limit,
            max_atoms=atom_limit,
            max_records=record_limit,
            max_line_len=line_limit,
        ),
    )
    atoms: list[PdbAtomRecord] = []
    for number, line in enumerate(text.splitlines(), 1):
        record = line[:6].rstrip(" ")
        if record not in ("ATOM", "HETATM"):
            continue
        atoms.append(
            PdbAtomRecord(
                len(atoms),
                number,
                int(line[6:11]),
                record,
                line[12:16],
                line[16],
                line[17:20].strip(),
                line[21],
                int(line[22:26]),
                line[26],
                (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                _optional_number(line[54:60]),
                _optional_number(line[60:66]),
                line[76:78],
                line[78:80],
                line,
            )
        )
    if len(atoms) != validation["source_atom_count"]:
        raise ValueError("PDB atom row correspondence disagrees with native validation")
    return PdbAtomRead(
        text,
        validation["model_num"],
        tuple(atoms),
        tuple((row["line"], row["record_name"]) for row in validation["losses"]),
    )
