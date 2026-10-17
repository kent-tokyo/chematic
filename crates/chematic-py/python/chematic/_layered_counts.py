"""Typed, staged caller-owned counters for the Layered fingerprint adapter."""

from __future__ import annotations

from dataclasses import dataclass
import operator
from typing import (
    Any,
    Iterable,
    List,
    Optional,
    Protocol,
    Sequence,
    SupportsIndex,
    TypedDict,
    Union,
)

import numpy as np

LayeredAtomCounts = Union[List[int], np.ndarray[Any, np.dtype[np.uint32]]]


class LayeredBitMask(Protocol):
    """Read-only explicit or sparse binary bit vector used as an output mask."""

    def GetNumBits(self) -> SupportsIndex: ...
    def GetOnBits(self) -> Iterable[SupportsIndex]: ...


class LayeredOptions(TypedDict):
    layer_flags: int
    min_path: int
    max_path: int
    nbits: int
    branched_paths: bool
    from_atoms: Optional[List[int]]
    set_only_bits: Optional[bytes]
    atom_counts: Optional[List[int]]
    max_visits: int
    max_paths: int


class LayeredResult(TypedDict):
    bits: bytes
    atom_counts: List[int]


def integer(value: SupportsIndex, name: str, lower: int, upper: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer, excluding bool")
    result = operator.index(value)
    if not lower <= result <= upper:
        raise ValueError(f"{name} must be in {lower}..{upper}")
    return result


def _destination(value: LayeredAtomCounts, minimum: int) -> int:
    if type(value) is list:
        length = len(value)
    elif type(value) is np.ndarray:
        if value.ndim != 1 or value.dtype.kind != "u" or value.dtype.itemsize != 4:
            raise TypeError("atomCounts requires a one-dimensional uint32 ndarray")
        if not value.flags.writeable:
            raise ValueError("atomCounts array must be writable")
        length = len(value)
        if length > 1 and abs(value.strides[0]) < value.dtype.itemsize:
            raise ValueError("atomCounts array elements must not overlap")
    else:
        raise TypeError("atomCounts must be a plain list or plain uint32 ndarray")
    if not minimum <= length <= 65536:
        raise ValueError(
            "atomCounts must cover all atoms and have at most 65536 entries"
        )
    return length


@dataclass
class Counters:
    destination: Optional[LayeredAtomCounts]
    seed: Optional[List[int]]

    @classmethod
    def prepare(cls, destination: Optional[LayeredAtomCounts], atoms: int) -> Counters:
        if destination is None:
            return cls(None, None)
        _destination(destination, atoms)
        return cls(
            destination,
            [integer(v, "atomCounts entry", 0, 2**32 - 1) for v in destination],
        )

    def commit(self, values: Sequence[int]) -> None:
        if self.destination is None:
            return
        expected = len(self.seed) if self.seed is not None else 0
        if len(values) != expected or _destination(self.destination, 0) != expected:
            raise ValueError("atomCounts length changed during calculation")
        staged = [integer(v, "computed atomCounts entry", 0, 2**32 - 1) for v in values]
        if type(self.destination) is list:
            self.destination[:] = staged
        else:
            # Allocate/convert the complete result before the first write. Assignment
            # to a plain, non-overlapping uint32 view has no user-defined callbacks.
            array = np.asarray(staged, dtype=np.uint32)
            self.destination[:] = array


def _roots(
    values: Optional[Sequence[SupportsIndex]], atoms: int
) -> Optional[List[int]]:
    if values is None:
        return None
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise TypeError("fromAtoms must be an indexed integer sequence")
    size = len(values)
    if size > 4096:
        raise ValueError("fromAtoms has more than 4096 entries")
    return [integer(values[i], "fromAtoms entry", 0, atoms - 1) for i in range(size)]


def _mask(value: Optional[LayeredBitMask], width: int) -> Optional[bytes]:
    if value is None:
        return None
    actual = integer(value.GetNumBits(), "mask width", 1, 65536)
    if actual != width:
        raise ValueError("setOnlyBits width must equal fpSize")
    packed = bytearray((width + 7) // 8)
    seen = set()
    for raw in value.GetOnBits():
        bit = integer(raw, "mask bit", 0, width - 1)
        if bit in seen:
            raise ValueError("setOnlyBits contains a duplicate bit")
        seen.add(bit)
        packed[bit // 8] |= 1 << (bit % 8)
    return bytes(packed)


def prepare(
    atoms: int,
    layer_flags: SupportsIndex,
    min_path: SupportsIndex,
    max_path: SupportsIndex,
    fp_size: SupportsIndex,
    destination: Optional[LayeredAtomCounts],
    mask: Optional[LayeredBitMask],
    branched: bool,
    roots: Optional[Sequence[SupportsIndex]],
    max_visits: SupportsIndex,
    max_paths: SupportsIndex,
) -> tuple[LayeredOptions, Counters]:
    if not isinstance(branched, bool):
        raise TypeError("branchedPaths must be bool")
    width = integer(fp_size, "fpSize", 1, 65536)
    low = integer(min_path, "minPath", 1, 12)
    high = integer(max_path, "maxPath", low, 12)
    options: LayeredOptions = {
        "layer_flags": integer(layer_flags, "layerFlags", 0, 2**32 - 1),
        "min_path": low,
        "max_path": high,
        "nbits": width,
        "branched_paths": branched,
        "from_atoms": _roots(roots, atoms),
        "set_only_bits": _mask(mask, width),
        "atom_counts": None,
        "max_visits": integer(max_visits, "maxVisits", 1, 10_000_000),
        "max_paths": integer(max_paths, "maxPaths", 1, 1_000_000),
    }
    counters = Counters.prepare(destination, atoms)
    options["atom_counts"] = counters.seed
    return options, counters
