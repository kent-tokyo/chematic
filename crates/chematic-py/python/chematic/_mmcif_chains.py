"""Immutable model/label-chain views over the native mmCIF residue reader."""

from __future__ import annotations

from dataclasses import dataclass
import operator
from typing import Optional, SupportsIndex, TypedDict, cast

from .chematic import parse_mmcif_residues

_MAX_BYTES = 8 * 1024 * 1024
_MAX_ATOMS = 50_000
_MAX_CHAINS = 4_096
_MAX_LINE = 8_192


@dataclass(frozen=True)
class MmcifChainResidue:
    """Full author/label residue identity; placeholders remain None."""

    model_num: int
    label_asym_id: str
    auth_asym_id: Optional[str]
    entity_id: Optional[str]
    label_seq_id: Optional[int]
    auth_seq_id: Optional[str]
    insertion_code: Optional[str]
    label_comp_id: Optional[str]
    auth_comp_id: Optional[str]


@dataclass(frozen=True)
class MmcifChainAtom:
    """One source atom_site row, with zero-based source index and angstroms."""

    source_index: int
    serial: int
    group_pdb: str
    element: str
    label_atom_id: Optional[str]
    auth_atom_id: Optional[str]
    alt_loc: Optional[str]
    position_angstrom: tuple[float, float, float]
    occupancy: float
    b_iso: float
    formal_charge: Optional[int]
    residue: MmcifChainResidue


@dataclass(frozen=True)
class MmcifChain:
    """A model/label_asym_id partition; no connectivity is inferred."""

    model_num: int
    label_asym_id: str
    atoms: tuple[MmcifChainAtom, ...]


@dataclass(frozen=True)
class MmcifChainRead:
    """Ordered immutable chains and the complete unmodified input text."""

    source_text: str
    chains: tuple[MmcifChain, ...]
    unhandled_columns: tuple[str, ...]


class _Residue(TypedDict):
    model_num: int
    label_asym_id: Optional[str]
    auth_asym_id: Optional[str]
    entity_id: Optional[str]
    label_seq_id: Optional[int]
    auth_seq_id: Optional[str]
    insertion_code: Optional[str]
    label_comp_id: Optional[str]
    auth_comp_id: Optional[str]


class _Atom(TypedDict):
    serial: int
    group_pdb: str
    element: str
    label_atom_id: Optional[str]
    auth_atom_id: Optional[str]
    alt_loc: Optional[str]
    x: float
    y: float
    z: float
    occupancy: float
    b_iso: float
    formal_charge: Optional[int]


class _Document(TypedDict):
    atoms: list[_Atom]
    residues: list[_Residue]
    atom_to_residue: list[int]
    unhandled_columns: list[str]


def _limit(value: SupportsIndex, name: str, maximum: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} requires an integer, excluding bool")
    number = operator.index(value)
    if not 1 <= number <= maximum:
        raise ValueError(f"{name} must be in 1..{maximum}")
    return number


def read_mmcif_chains(
    text: str,
    *,
    max_input_bytes: SupportsIndex = _MAX_BYTES,
    max_atoms: SupportsIndex = _MAX_ATOMS,
    max_chains: SupportsIndex = _MAX_CHAINS,
    max_line_len: SupportsIndex = _MAX_LINE,
) -> MmcifChainRead:
    """Read model/label-chain records from the first native atom_site loop.

    Args:
        text: mmCIF contents, not a path. Complete source text is retained.
        max_input_bytes: UTF-8 byte limit, 1..8388608.
        max_atoms: Parsed atom limit, 1..50000, across all models.
        max_chains: Number of model/label-chain pairs, 1..4096.
        max_line_len: Native parser line-byte limit, 1..8192.

    Returns:
        Immutable chains ordered by first occurrence, with atoms in source
        row order. Different models never merge. Both author and label
        identities, alternate locations and native numeric fields are kept.
        Unhandled atom_site column names and complete source remain available.

    Raises:
        TypeError: Non-string input or noninteger limits (including bool).
        ValueError: Exceeded limits, malformed/native-unsupported input,
            missing label_asym_id or ambiguous/duplicate residue atoms.

    Requires a label chain ID on every atom; author aliases are never used to
    merge chains. Uses the existing bounded native residue parser, including
    its first-loop scope and numeric defaults. This is an eager atom-record
    view, not a molecule supplier, bond perception, CIF writer or biological
    assembly generator. Other loops/categories are retained only in source_text.
    """
    if not isinstance(text, str):
        raise TypeError("text must be mmCIF str contents")
    byte_limit = _limit(max_input_bytes, "max_input_bytes", _MAX_BYTES)
    atom_limit = _limit(max_atoms, "max_atoms", _MAX_ATOMS)
    chain_limit = _limit(max_chains, "max_chains", _MAX_CHAINS)
    line_limit = _limit(max_line_len, "max_line_len", _MAX_LINE)
    if len(text) > byte_limit or len(text.encode("utf-8")) > byte_limit:
        raise ValueError("mmCIF chain input byte limit exceeded")
    doc = cast(
        _Document,
        parse_mmcif_residues(
            text,
            max_input_bytes=byte_limit,
            max_atoms=atom_limit,
            max_line_len=line_limit,
        ),
    )
    identities = []
    for row in doc["residues"]:
        label = row["label_asym_id"]
        if label is None or not label.strip():
            raise ValueError("mmCIF chain reading requires label_asym_id on every atom")
        identities.append(
            MmcifChainResidue(
                row["model_num"],
                label,
                row["auth_asym_id"],
                row["entity_id"],
                row["label_seq_id"],
                row["auth_seq_id"],
                row["insertion_code"],
                row["label_comp_id"],
                row["auth_comp_id"],
            )
        )
    groups: dict[tuple[int, str], list[MmcifChainAtom]] = {}
    for index, atom in enumerate(doc["atoms"]):
        residue = identities[doc["atom_to_residue"][index]]
        key = (residue.model_num, residue.label_asym_id)
        if key not in groups:
            if len(groups) >= chain_limit:
                raise ValueError("mmCIF chain count limit exceeded")
            groups[key] = []
        groups[key].append(
            MmcifChainAtom(
                index,
                atom["serial"],
                atom["group_pdb"],
                atom["element"],
                atom["label_atom_id"],
                atom["auth_atom_id"],
                atom["alt_loc"],
                (atom["x"], atom["y"], atom["z"]),
                atom["occupancy"],
                atom["b_iso"],
                atom["formal_charge"],
                residue,
            )
        )
    return MmcifChainRead(
        text,
        tuple(
            MmcifChain(model, label, tuple(atoms))
            for (model, label), atoms in groups.items()
        ),
        tuple(doc["unhandled_columns"]),
    )
