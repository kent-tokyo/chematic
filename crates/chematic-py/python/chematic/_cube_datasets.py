"""Bounded CUBE dataset demultiplexing; independently based on h5cube's grammar."""

from __future__ import annotations

from dataclasses import dataclass, field
import io
import math
import operator
import re
from typing import TYPE_CHECKING, Iterator, Optional, SupportsIndex

import numpy as np

if TYPE_CHECKING:
    from . import VolumetricGrid
else:
    from .chematic import VolumetricGrid

_MAX_BYTES = 8 * 1024 * 1024
_MAX_ATOMS = 50_000
_MAX_POINTS = 262_144
_MAX_VALUES = 262_144
_MAX_DATASETS = 64
_MAX_LINE = 16_384
_MAX_OUTPUT = 16 * 1024 * 1024
_INT = re.compile(r"[+-]?[0-9]+\Z")
_FLOAT = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")


def _integer(token: str) -> int:
    if not _INT.fullmatch(token) or len(token) > 20:
        raise ValueError("CUBE requires a bounded decimal integer")
    value = int(token)
    if not -(2**63) <= value < 2**63:
        raise ValueError("CUBE integer is outside signed 64-bit range")
    return value


def _limit(value: int, maximum: int, name: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be a Python int, excluding bool")
    if not 1 <= value <= maximum:
        raise ValueError(f"{name} must be in 1..{maximum}")
    return value


def _lines(text: str) -> Iterator[str]:
    for line in io.StringIO(text, newline=None):
        line = line.rstrip("\n")
        if len(line) > _MAX_LINE:
            raise ValueError("CUBE line exceeds 16384 characters")
        yield line


def _required(lines: Iterator[str]) -> str:
    try:
        return next(lines)
    except StopIteration:
        raise ValueError("CUBE header/geometry is truncated") from None


def _tokens(lines: Iterator[str]) -> Iterator[str]:
    for line in lines:
        yield from line.split()


def _token(tokens: Iterator[str]) -> str:
    try:
        return next(tokens)
    except StopIteration:
        raise ValueError("CUBE dataset identifiers/data are truncated") from None


def _render(header: tuple[str, ...], values: tuple[float, ...]) -> str:
    lines = list(header)
    size = sum(len(line.encode("utf-8")) + 1 for line in lines)
    for offset in range(0, len(values), 6):
        line = " ".join(format(value, ".17g") for value in values[offset : offset + 6])
        size += len(line) + 1
        if size > _MAX_OUTPUT:
            raise ValueError("selected CUBE output exceeds 16 MiB")
        lines.append(line)
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class CubeDatasetSelection:
    """Explicit single-dataset export, retaining complete source and omission IDs.

    text has nonnegative NATOMS, no NVAL or dataset-ID record, and unchanged
    comments/geometry lines. Float64 scalar values are rendered with 17 digits.
    Dataset IDs are opaque identifiers, not assigned physical orbital meanings.
    """

    source_text: str = field(repr=False)
    dataset_index: int
    dataset_id: Optional[int]
    omitted_indices: tuple[int, ...]
    text: str = field(repr=False)

    def to_grid(self) -> VolumetricGrid:
        """Parse the selected scalar field into a fresh native grid, in bohr."""
        return VolumetricGrid.from_cube(
            self.text,
            max_input_bytes=_MAX_OUTPUT,
            max_atoms=_MAX_ATOMS,
            max_grid_points=_MAX_POINTS,
        )


@dataclass(frozen=True)
class CubeDatasetDocument:
    """Immutable parsed datasets, one tuple per dataset, third axis fastest.

    Positive NATOMS/NVAL input has dataset_ids=None: positional indices are
    not invented orbital identifiers. Original source text is retained exactly.
    Parser results use nested tuples and cannot be mutated through NumPy views.
    """

    source_text: str = field(repr=False)
    shape: tuple[int, int, int]
    dataset_ids: Optional[tuple[int, ...]]
    datasets: tuple[tuple[float, ...], ...] = field(repr=False)
    _header: tuple[str, ...] = field(repr=False)

    @property
    def dataset_count(self) -> int:
        return len(self.datasets)

    def select(self, index: SupportsIndex) -> CubeDatasetSelection:
        """Export one zero-based dataset; no implicit first-dataset truncation.

        index must be in 0..dataset_count-1; bool/negative indices are refused.
        The selection records omitted positional indices and the original text.
        """
        if isinstance(index, (bool, np.bool_)):
            raise TypeError("dataset index must be an integer, excluding bool")
        selected = operator.index(index)
        if not 0 <= selected < self.dataset_count:
            raise ValueError("dataset index is outside this CUBE document")
        result = CubeDatasetSelection(
            self.source_text,
            selected,
            None if self.dataset_ids is None else self.dataset_ids[selected],
            tuple(i for i in range(self.dataset_count) if i != selected),
            _render(self._header, self.datasets[selected]),
        )
        # Keep native atom/header/finite geometry validation, including on
        # explicitly constructed document instances; no cached mutable grid.
        result.to_grid()
        return result


def read_cube_datasets(
    text: str,
    *,
    max_input_bytes: int = _MAX_BYTES,
    max_atoms: int = _MAX_ATOMS,
    max_grid_points: int = _MAX_POINTS,
    max_total_values: int = _MAX_VALUES,
) -> CubeDatasetDocument:
    """Read positive NATOMS/NVAL or negative NATOMS/dataset-ID CUBE input.

    Dataset index is the fastest-varying input dimension. IDs may wrap across
    lines, must be unique signed 64-bit integers, and are retained in order.
    Negative NATOMS permits absent NVAL or NVAL=1. At most 64 datasets, 8 MiB
    UTF-8 input, 50000 atoms and 262144 points/total values; limits may only be
    lowered. Finite decimal/e-notation only, with positive axis counts (bohr).
    Native geometry validation applies; ghost atoms and ambiguous negative-axis
    units are refused. No resampling, scalar unit conversion or orbital inference.
    """
    byte_limit = _limit(max_input_bytes, _MAX_BYTES, "max_input_bytes")
    atom_limit = _limit(max_atoms, _MAX_ATOMS, "max_atoms")
    point_limit = _limit(max_grid_points, _MAX_POINTS, "max_grid_points")
    value_limit = _limit(max_total_values, _MAX_VALUES, "max_total_values")
    if type(text) is not str:
        raise TypeError("CUBE input must be UTF-8 text, not a path or bytes")
    if len(text) > byte_limit:
        raise ValueError("CUBE input exceeds byte limit")
    try:
        if len(text.encode("utf-8")) > byte_limit:
            raise ValueError("CUBE input exceeds byte limit")
    except UnicodeError as exc:
        raise ValueError("CUBE text must be valid UTF-8") from exc
    if any(ord(c) < 32 and c not in "\t\r\n" for c in text):
        raise ValueError("CUBE input contains unsupported control characters")
    lines = _lines(text)
    header = [_required(lines), _required(lines)]
    count = _required(lines).split()
    if len(count) not in (4, 5):
        raise ValueError("CUBE count line requires NATOMS origin [NVAL]")
    natoms = _integer(count[0])
    if count[0].startswith("-") and natoms == 0:
        raise ValueError("negative-zero NATOMS has ambiguous dataset semantics")
    if abs(natoms) > atom_limit:
        raise ValueError("CUBE atom count exceeds limit")
    nval = _integer(count[4]) if len(count) == 5 else 1
    if not 1 <= nval <= _MAX_DATASETS or (natoms < 0 and nval != 1):
        raise ValueError("invalid CUBE NVAL for atom-count sign or dataset limit")
    header.append(" ".join([str(abs(natoms)), *count[1:4]]))
    shape = []
    for _ in range(3):
        line = _required(lines)
        fields = line.split()
        if len(fields) != 4:
            raise ValueError("CUBE axis line requires count and three components")
        axis = _integer(fields[0])
        if axis <= 0:
            raise ValueError("CUBE dataset profile requires positive bohr axis counts")
        shape.append(axis)
        header.append(line)
    points = math.prod(shape)
    if points > point_limit:
        raise ValueError("CUBE grid point count exceeds limit")
    header.extend(_required(lines) for _ in range(abs(natoms)))
    tokens = _tokens(lines)
    ids: Optional[tuple[int, ...]] = None
    if natoms < 0:
        nval = _integer(_token(tokens))
        if not 1 <= nval <= _MAX_DATASETS:
            raise ValueError("CUBE dataset count must be in 1..64")
        ids = tuple(_integer(_token(tokens)) for _ in range(nval))
        if len(set(ids)) != nval:
            raise ValueError("CUBE dataset identifiers must be unique")
    if points * nval > value_limit:
        raise ValueError("CUBE total scalar count exceeds limit")
    datasets: list[list[float]] = [[] for _ in range(nval)]
    for scalar in range(points * nval):
        token = _token(tokens)
        if not _FLOAT.fullmatch(token):
            raise ValueError("CUBE scalar requires finite decimal/e notation")
        value = float(token)
        if not math.isfinite(value):
            raise ValueError("CUBE scalar must be finite")
        datasets[scalar % nval].append(value)
    if next(tokens, None) is not None:
        raise ValueError("CUBE contains trailing data after expected scalars")
    result = CubeDatasetDocument(
        text,
        (shape[0], shape[1], shape[2]),
        ids,
        tuple(tuple(values) for values in datasets),
        tuple(header),
    )
    result.select(0)  # Validate the common header/geometry with the native reader.
    return result
