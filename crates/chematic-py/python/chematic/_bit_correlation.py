"""Bounded binary bit-pair co-occurrence counts, without normalization."""

from __future__ import annotations

import operator
from typing import Iterable, SupportsIndex

import numpy as np

from ._numpy_fingerprint import _BitVectorSource as BinaryBitVector

_MAX_BITS = 1_024
_MAX_IDS = 65_536
_MAX_WIDTH = 2**32 - 1
_MAX_RECORDS = 1_000_000
_MAX_PAIR_UPDATES = 50_000_000


def _index(value: SupportsIndex) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("bit IDs and widths must be integers, excluding bool")
    return operator.index(value)


class BitCorrMatGenerator:
    """Checked InfoTheory co-occurrence counts for explicit/sparse binary sources.

    Select 2..1024 unique bit IDs with SetBitList before collecting votes.
    Each valid record adds one to every selected pair whose two bits are on.
    Output is an owned float64 strict lower triangle, with pair positions
    following the caller's selected-bit order. This is neither a normalized
    Pearson correlation matrix nor a chemical fingerprint generator.
    """

    def __init__(self) -> None:
        self._positions: dict[int, int] = {}
        self._counts: np.ndarray[tuple[int], np.dtype[np.uint64]] = np.zeros(
            0, dtype=np.uint64
        )
        self._width: int | None = None
        self._records = 0
        self._pair_updates = 0

    def SetBitList(self, bitList: Iterable[SupportsIndex]) -> None:
        """Select ordered unique IDs; success resets all votes and source width.

        IDs are in 0..4294967294, excluding bool. Iteration consumes at most
        1025 IDs to detect the cap. A failed setter leaves the previous list,
        counts and resource counters unchanged. The input is not retained.
        """
        if isinstance(bitList, (str, bytes, bytearray)):
            raise TypeError("bitList must be an iterable of integer IDs")
        positions: dict[int, int] = {}
        for value in bitList:
            if len(positions) >= _MAX_BITS:
                raise ValueError("bitList exceeds 1024 selected bits")
            bit = _index(value)
            if not 0 <= bit < _MAX_WIDTH or bit in positions:
                raise ValueError("selected IDs must be unique and in 0..4294967294")
            positions[bit] = len(positions)
        count = len(positions)
        if count < 2:
            raise ValueError("bitList requires at least two selected bits")
        counts = np.zeros(count * (count - 1) // 2, dtype=np.uint64)
        self._positions = positions
        self._counts = counts
        self._width = None
        self._records = 0
        self._pair_updates = 0

    def CollectVotes(self, bitVect: BinaryBitVector) -> None:
        """Count selected pairs from one fully validated binary vector.

        Sources provide GetNumBits()/GetOnBits(). Widths must be consistent,
        in 1..4294967295, and contain every selected bit. Raw IDs must be
        unique and within that width; at most 65536 are consumed, plus one to
        detect overflow. All IDs, including unselected IDs, are checked.
        Bool widths/IDs refuse. Source callbacks/iteration errors propagate;
        AttributeError during method lookup is treated as a missing method.
        Validation/limit errors do not change votes or counters.

        Per configuration, at most one million records and 50 million pair
        increments. The source remains alive during extraction, but is not
        retained after return. Caller methods have no timeout; concurrent
        mutations of a generator instance are outside this contract.
        """
        if not self._positions:
            raise ValueError("SetBitList must configure selected bits first")
        if self._records >= _MAX_RECORDS:
            raise ValueError("co-occurrence record limit exceeded")
        get_width = getattr(bitVect, "GetNumBits", None)
        get_bits = getattr(bitVect, "GetOnBits", None)
        if not callable(get_width) or not callable(get_bits):
            raise TypeError("bitVect requires GetNumBits() and GetOnBits()")
        width = _index(get_width())
        if not 1 <= width <= _MAX_WIDTH or max(self._positions) >= width:
            raise ValueError("source width must contain selected IDs and fit u32")
        if self._width is not None and width != self._width:
            raise ValueError("source bit widths must remain consistent")
        values = get_bits()
        if isinstance(values, (str, bytes, bytearray)):
            raise TypeError("GetOnBits() must return integer positions, not text")
        seen: set[int] = set()
        for value in values:
            if len(seen) >= _MAX_IDS:
                raise ValueError("GetOnBits exceeds 65536 raw IDs")
            bit = _index(value)
            if not 0 <= bit < width or bit in seen:
                raise ValueError("on-bit IDs must be unique and within source width")
            seen.add(bit)
        active = sorted(self._positions[bit] for bit in seen if bit in self._positions)
        updates = len(active) * (len(active) - 1) // 2
        pair_updates = self._pair_updates + updates
        if pair_updates > _MAX_PAIR_UPDATES:
            raise ValueError("co-occurrence exceeds 50000000 pair updates")
        records = self._records + 1
        indices = np.fromiter(
            (
                i * (i - 1) // 2 + j
                for pos, i in enumerate(active)
                for j in active[:pos]
            ),
            dtype=np.intp,
            count=updates,
        )
        # Allocate new values before changing state. Pair indices are unique;
        # uint64 counts <= one million remain exact in the float64 output.
        changed = self._counts[indices] + np.uint64(1)
        self._counts[indices] = changed
        self._width = width
        self._records = records
        self._pair_updates = pair_updates

    def GetCorrMatrix(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Return independent pair counts (1,0), (2,0), (2,1), ... .

        Requires a selected-bit configuration; before votes the array is zero.
        At most 523776 entries, without diagonal or implicit unit normalization.
        No source object or writable view of internal counts is returned.
        """
        if not self._positions:
            raise ValueError("SetBitList must configure selected bits first")
        return self._counts.astype(np.float64)
