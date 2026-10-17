"""Bounded binary fingerprint projection into a validated NumPy destination."""

from __future__ import annotations

import operator
from typing import Iterable, Protocol, SupportsIndex

import numpy as np

_MAX_BITS = 1_048_576
_MAX_INPUT_IDS = 1_048_576


class _BitVectorSource(Protocol):
    def GetNumBits(self) -> SupportsIndex: ...
    def GetOnBits(self) -> Iterable[SupportsIndex]: ...


def _index(value: SupportsIndex) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("fingerprint widths and bit IDs must be integers, not booleans")
    return operator.index(value)


def _destination(dest: np.ndarray, width: int) -> None:
    # Subclass setters/array hooks must not run halfway through committing bits.
    if type(dest) is not np.ndarray:
        raise TypeError("dest must be a plain numpy.ndarray")
    if dest.ndim != 1 or dest.shape[0] != width:
        raise ValueError(
            "dest must be one-dimensional with exactly GetNumBits() elements"
        )
    if dest.dtype.kind not in "biuf":
        raise TypeError("dest dtype must be boolean, integer or floating-point")
    if not dest.flags.writeable:
        raise ValueError("dest must be writable")
    if width > 1 and abs(dest.strides[0]) < dest.dtype.itemsize:
        raise ValueError("dest elements must not overlap")


def _convert(fp: _BitVectorSource, dest: np.ndarray | None) -> np.ndarray:
    try:
        size = fp.GetNumBits
        on_bits = fp.GetOnBits
    except AttributeError as error:
        raise TypeError("fp must provide GetNumBits() and GetOnBits()") from error
    if not callable(size) or not callable(on_bits):
        raise TypeError("fp must provide callable GetNumBits() and GetOnBits()")
    width = _index(size())
    if not 0 <= width <= _MAX_BITS:
        raise ValueError(f"fingerprint width must be between 0 and {_MAX_BITS}")
    if dest is not None:
        _destination(dest, width)
    # Extract/validate every source ID before touching a caller's destination.
    # Keep fp alive while consuming a method-produced lazy view.
    bits = np.zeros(width, dtype=np.int8)
    for count, value in enumerate(on_bits()):
        if count == _MAX_INPUT_IDS:
            raise ValueError(f"fingerprint exceeds {_MAX_INPUT_IDS} input bit IDs")
        bit = _index(value)
        if not 0 <= bit < width:
            raise ValueError("on-bit ID is outside GetNumBits()")
        bits[bit] = 1
    if dest is None:
        return bits
    # Every value is 0/1 and is exactly representable by all accepted dtypes.
    dest[...] = bits
    return dest
