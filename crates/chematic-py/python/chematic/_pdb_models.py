"""Checked coordinate-only PDB model splitting (Open Babel #2476)."""

from __future__ import annotations

from dataclasses import dataclass
import operator
import re
from typing import List, Optional, SupportsIndex, Tuple, TypedDict, cast

from .chematic import select_pdb_atoms

_MAX_BYTES = 64 * 1024 * 1024
_MAX_RECORDS = 1_000_000
_MAX_MODELS = 9999
_COORDINATES = {"ATOM", "HETATM", "ANISOU", "TER"}


@dataclass(frozen=True)
class PdbSplitLoss:
    """Metadata omitted from all outputs, identified in the original input."""

    line: int
    record_name: str


@dataclass(frozen=True)
class PdbModelFile:
    """One validated PDB coordinate document and its relative file name."""

    filename: str
    model_num: Optional[int]
    pdb: str
    atom_count: int


@dataclass(frozen=True)
class PdbModelSplit:
    """A complete, immutable in-memory split; no files have been written."""

    source_text: str
    files: Tuple[PdbModelFile, ...]
    losses: Tuple[PdbSplitLoss, ...]


class _Selection(TypedDict):
    source_atom_count: int


def _limit(value: SupportsIndex, name: str, maximum: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer, excluding bool")
    number = operator.index(value)
    if not 1 <= number <= maximum:
        raise ValueError(f"{name} must be in 1..{maximum}")
    return number


def _fail(line: int, message: str) -> ValueError:
    return ValueError(f"PDB split line {line}: {message}")


def _model_number(line: str, number: int) -> int:
    field = line[10:14].strip()
    if (
        len(line) < 14
        or line[6:10].strip()
        or not field.isdecimal()
        or int(field) == 0
        or line[14:].strip()
    ):
        raise _fail(number, "invalid positive four-column MODEL number")
    return int(field)


@dataclass
class _Parts:
    models: List[Tuple[Optional[int], List[str]]]
    conect: List[str]
    losses: List[PdbSplitLoss]


def _partition(text: str, allow_loss: bool, max_models: int) -> _Parts:
    parts = _Parts([], [], [])
    current: Optional[List[str]] = None
    implicit: List[str] = []
    model_num: Optional[int] = None
    seen = set()
    connected = ended = False
    # Bound record count before allocating the line list, including blank lines.
    count = text.count("\n") + int(bool(text) and not text.endswith("\n"))
    if count > _MAX_RECORDS:
        raise ValueError("PDB split input record limit exceeded")
    for number, raw in enumerate(text.split("\n"), 1):
        line = raw[:-1] if raw.endswith("\r") and number <= count else raw
        if len(line) > 1024 or any(not " " <= c <= "~" for c in line):
            raise _fail(number, "records require printable ASCII, at most 1024 bytes")
        if not line.strip():
            continue
        if ended:
            raise _fail(number, "records after END are unsupported")
        name = line[:6].rstrip()
        if name == "MODEL":
            if current is not None or implicit or connected:
                raise _fail(
                    number, "nested MODEL or mixed implicit/explicit coordinates"
                )
            if len(parts.models) >= max_models:
                raise _fail(number, "model count limit exceeded")
            model_num = _model_number(line, number)
            if model_num in seen:
                raise _fail(number, "duplicate MODEL number")
            seen.add(model_num)
            current = [line]
        elif name == "ENDMDL":
            if current is None or line[6:].strip():
                raise _fail(number, "unmatched or malformed ENDMDL")
            current.append(line)
            parts.models.append((model_num, current))
            current = None
        elif name in _COORDINATES:
            if connected or (parts.models and current is None):
                raise _fail(number, "coordinate record outside model")
            (current if current is not None else implicit).append(line)
        elif name == "CONECT":
            if current is not None:
                raise _fail(number, "CONECT must follow all ENDMDL records")
            connected = True
            parts.conect.append(line)
        elif name == "END":
            if current is not None or line[3:].strip():
                raise _fail(number, "unclosed MODEL or malformed END")
            ended = True
        else:
            if not allow_loss:
                raise _fail(number, f"record {name!r} requires metadata-loss consent")
            parts.losses.append(PdbSplitLoss(number, name))
    if current is not None:
        raise ValueError("PDB split MODEL lacks ENDMDL")
    if implicit:
        parts.models.append((None, implicit))
    if not parts.models:
        raise ValueError("PDB split requires at least one nonempty coordinate model")
    return parts


def split_pdb_models(
    text: str,
    *,
    prefix: str = "split",
    allow_metadata_loss: bool = False,
    max_models: SupportsIndex = _MAX_MODELS,
    max_output_bytes: SupportsIndex = _MAX_BYTES,
) -> PdbModelSplit:
    """Split checked coordinate models into padded names, without file I/O.

    Width is the digit count of the total, numbering outputs from one in source
    order. Original MODEL numbers and coordinate/CONECT lines are retained.
    Other records need explicit loss consent; their original line numbers and
    the entire input remain in the result. Output uses LF and one final END.
    All models and aggregate output limits are validated before returning.
    """
    if not isinstance(text, str):
        raise TypeError("text must be str")
    if len(text) > _MAX_BYTES:
        raise ValueError("PDB split input byte limit exceeded")
    if (
        not isinstance(prefix, str)
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", prefix) is None
    ):
        raise ValueError(
            "prefix must be a 1..64 character ASCII file stem without paths"
        )
    if not isinstance(allow_metadata_loss, bool):
        raise TypeError("allow_metadata_loss must be bool")
    model_limit = _limit(max_models, "max_models", _MAX_MODELS)
    output_limit = _limit(max_output_bytes, "max_output_bytes", _MAX_BYTES)
    parts = _partition(text, allow_metadata_loss, model_limit)
    common_size = sum(len(line) + 1 for line in parts.conect) + 4  # END\n
    output_size = len(parts.models) * common_size + sum(
        len(line) + 1 for _, rows in parts.models for line in rows
    )
    if output_size > output_limit:
        raise ValueError("PDB split aggregate output byte limit exceeded")
    tail = "\n".join(parts.conect) + "\nEND\n" if parts.conect else "END\n"
    width = len(str(len(parts.models)))
    files = []
    for index, (model_num, rows) in enumerate(parts.models, 1):
        pdb = "\n".join(rows) + "\n" + tail
        try:
            checked = cast(_Selection, select_pdb_atoms(pdb))
        except ValueError as error:
            raise ValueError(
                f"PDB split output {index} (MODEL {model_num}): {error}"
            ) from error
        atom_count = checked["source_atom_count"]
        if atom_count == 0:
            raise ValueError(f"PDB split output {index} has no ATOM/HETATM records")
        files.append(
            PdbModelFile(f"{prefix}{index:0{width}d}.pdb", model_num, pdb, atom_count)
        )
    return PdbModelSplit(text, tuple(files), tuple(parts.losses))
