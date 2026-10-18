"""Coordinate-only PDB assembly using strict native selection validation."""

from __future__ import annotations

from collections.abc import Mapping as RuntimeMapping
from dataclasses import dataclass, field
import string
from typing import TYPE_CHECKING, Iterable, Mapping, Optional

if TYPE_CHECKING:
    from . import select_pdb_atoms as _select
else:
    from .chematic import select_pdb_atoms as _select

_MAX_SOURCES = 32
_MAX_SOURCE_BYTES = 8 * 1024 * 1024
_MAX_BYTES = 16 * 1024 * 1024
_MAX_ATOMS = 50_000
_MAX_RECORDS = 250_000
_CHAINS = frozenset(string.ascii_letters + string.digits + " ")


@dataclass(frozen=True)
class PdbMergeAtom:
    """Atom correspondence; source and atom indices are zero-based input order."""

    source_index: int
    source_atom_index: int
    old_serial: int
    new_serial: int
    source_chain: str
    output_chain: str
    record_type: str


@dataclass(frozen=True)
class PdbMergeLoss:
    """Excluded full-entry metadata, with one-based original physical line."""

    source_index: int
    line: int
    record_name: str


@dataclass(frozen=True)
class PdbMergeResult:
    """Merged coordinate block and immutable original-source/mapping sidecars.

    Only atom/ANISOU/TER serial and chain fields change. CONECT is remapped
    within each source and moved after all coordinates, retaining multiplicity.
    Source boundaries lacking a final TER receive one, including HET-only
    blocks; this is a coordinate separator, not polymer/residue inference.
    """

    pdb: str = field(repr=False)
    source_texts: tuple[str, ...] = field(repr=False)
    atoms: tuple[PdbMergeAtom, ...]
    losses: tuple[PdbMergeLoss, ...]
    inserted_ter_after_sources: tuple[int, ...]


def _chain(value: str) -> str:
    if type(value) is not str:
        raise TypeError("chain identifiers must be strings")
    if value not in _CHAINS:
        raise ValueError("chain identifiers must be one ASCII letter, digit or space")
    return value


def _maps(
    values: Optional[Iterable[Mapping[str, str]]], count: int
) -> list[dict[str, str]]:
    if values is None:
        return [{} for _ in range(count)]
    result: list[dict[str, str]] = []
    for mapping in values:
        if len(result) >= _MAX_SOURCES:
            raise ValueError("chain_maps exceeds 32 entries")
        if not isinstance(mapping, RuntimeMapping):
            raise TypeError("chain_maps entries must be mappings")
        copied: dict[str, str] = {}
        for key in mapping:
            if len(copied) >= len(_CHAINS):
                raise ValueError("chain map exceeds the supported chain namespace")
            key = _chain(key)
            if key in copied:
                raise ValueError("chain mapping iteration repeats a key")
            copied[key] = _chain(mapping[key])
        result.append(copied)
    if len(result) != count:
        raise ValueError("chain_maps requires exactly one mapping per source")
    return result


def _sources(values: Iterable[str]) -> list[str]:
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError("sources requires an iterable of PDB text blocks")
    result: list[str] = []
    size = records = 0
    for text in values:
        if len(result) >= _MAX_SOURCES:
            raise ValueError("PDB merge exceeds 32 sources")
        if type(text) is not str:
            raise TypeError("each PDB source must be text")
        if len(text) > _MAX_SOURCE_BYTES:
            raise ValueError("PDB source exceeds 8 MiB")
        if not text.isascii():
            raise ValueError("PDB merge requires ASCII source text")
        size += len(text)
        records += text.count("\n") + int(bool(text) and not text.endswith("\n"))
        if size > _MAX_BYTES or records > _MAX_RECORDS:
            raise ValueError("PDB merge aggregate input byte/record limit exceeded")
        result.append(text)
    if not result:
        raise ValueError("PDB merge requires at least one source")
    return result


def _serial(value: int) -> int:
    if value > 99_999:
        raise ValueError("merged atom/TER serial exceeds the five-column PDB field")
    return value


def _patch(line: str, serial: int, chain: str) -> str:
    return line[:6] + f"{serial:5}" + line[11:21] + chain + line[22:]


def merge_pdb_coordinate_blocks(
    sources: Iterable[str],
    *,
    chain_maps: Optional[Iterable[Mapping[str, str]]] = None,
    allow_metadata_loss: bool = False,
) -> PdbMergeResult:
    """Merge 1..32 implicit single-model coordinate blocks without reclassification.

    Args:
        sources: ASCII PDB text blocks, each with at least one atom. Original
            input order and complete source texts are retained.
        chain_maps: Optional one mapping per input. Unmapped chains retain their
            ID; different source chain namespaces must remain disjoint. Blank
            chain is one space. Unknown keys/noninjective maps refuse.
        allow_metadata_loss: Explicitly consent to excluding unsupported entry
            records, with their original source/line/name in losses.

    Returns:
        An immutable result with renumbered atoms/TER, remapped ANISOU/CONECT,
        source-boundary TER records and one final END. Atom record kind, name,
        residue fields, element/charge text and numeric spelling are retained.

    Native selection validates each source and the complete output. MODEL
    wrappers, malformed supported records, hybrid-36 and colliding chains
    refuse, even with metadata-loss consent. No graph, bond, element, residue
    or hydrogen perception, alignment or cross-source connectivity is added.
    Caps: 8 MiB/source, 16 MiB aggregate input/output, 50000 atoms, 250000
    physical input/output records, 1024 bytes/line and five-column serials.
    Caller iterator/mapping exceptions propagate; no files are written.
    """
    if type(allow_metadata_loss) is not bool:
        raise TypeError("allow_metadata_loss must be bool")
    texts = _sources(sources)
    mappings = _maps(chain_maps, len(texts))
    coordinates: list[str] = []
    connections: list[str] = []
    atoms: list[PdbMergeAtom] = []
    losses: list[PdbMergeLoss] = []
    inserted: list[int] = []
    used_chains: set[str] = set()
    next_serial = 1
    output_size = 0

    def append(target: list[str], line: str) -> None:
        nonlocal output_size
        output_size += len(line) + 1
        if output_size > _MAX_BYTES:
            raise ValueError("merged PDB output exceeds 16 MiB")
        if len(coordinates) + len(connections) >= _MAX_RECORDS:
            raise ValueError("merged PDB output exceeds 250000 records")
        target.append(line)

    for source_index, (text, mapping) in enumerate(zip(texts, mappings)):
        report = _select(
            text,
            allow_metadata_loss=allow_metadata_loss,
            max_input_bytes=_MAX_SOURCE_BYTES,
            max_atoms=_MAX_ATOMS - len(atoms),
            max_records=_MAX_RECORDS,
            max_line_len=1024,
        )
        if report["model_num"] is not None:
            raise ValueError(
                "PDB merge requires implicit models; MODEL wrappers refuse"
            )
        if report["source_atom_count"] == 0:
            raise ValueError("each PDB merge source must have at least one atom")
        source_chains = {_chain(atom["chain_id"]) for atom in report["atoms"]}
        if not mapping.keys() <= source_chains:
            raise ValueError("chain map contains a chain absent from its source")
        renamed = {chain: mapping.get(chain, chain) for chain in source_chains}
        output_chains = set(renamed.values())
        if len(output_chains) != len(source_chains) or used_chains & output_chains:
            raise ValueError(
                "source chain namespaces collide; supply disjoint chain_maps"
            )
        used_chains.update(output_chains)
        losses.extend(
            PdbMergeLoss(source_index, loss["line"], loss["record_name"])
            for loss in report["losses"]
        )
        lines = report["pdb"].splitlines()
        serials: dict[int, int] = {}
        ters: dict[int, int] = {}
        atom_index = 0
        last_atom: Optional[str] = None
        for line_index, line in enumerate(lines):
            record = line[:6].strip()
            if record in ("ATOM", "HETATM"):
                old = int(line[6:11])
                serials[old] = _serial(next_serial)
                atoms.append(
                    PdbMergeAtom(
                        source_index,
                        atom_index,
                        old,
                        next_serial,
                        line[21],
                        renamed[line[21]],
                        record,
                    )
                )
                next_serial += 1
                atom_index += 1
                last_atom = line
            elif record == "TER":
                ters[line_index] = _serial(next_serial)
                next_serial += 1
                last_atom = None
        for line_index, line in enumerate(lines):
            record = line[:6].strip()
            if record in ("ATOM", "HETATM", "ANISOU"):
                append(
                    coordinates,
                    _patch(line, serials[int(line[6:11])], renamed[line[21]]),
                )
            elif record == "TER":
                append(coordinates, _patch(line, ters[line_index], renamed[line[21]]))
            elif record == "CONECT":
                refs = [
                    int(line[start : start + 5])
                    for start in range(6, 31, 5)
                    if line[start : start + 5].strip()
                ]
                append(
                    connections, "CONECT" + "".join(f"{serials[ref]:5}" for ref in refs)
                )
            elif record != "END":
                raise ValueError("unexpected record from native coordinate selection")
        if last_atom is not None:
            ter = f"TER   {_serial(next_serial):5}      {last_atom[17:27]}"
            append(coordinates, ter[:21] + renamed[last_atom[21]] + ter[22:])
            next_serial += 1
            inserted.append(source_index)
    append(connections, "END")
    pdb = "\n".join(coordinates + connections) + "\n"
    validated = _select(
        pdb,
        max_input_bytes=_MAX_BYTES,
        max_atoms=_MAX_ATOMS,
        max_records=_MAX_RECORDS,
        max_line_len=1024,
    )
    if validated["pdb"] != pdb:
        raise ValueError("native validation unexpectedly changed merged output")
    return PdbMergeResult(
        pdb, tuple(texts), tuple(atoms), tuple(losses), tuple(inserted)
    )
