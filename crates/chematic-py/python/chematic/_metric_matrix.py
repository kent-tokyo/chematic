"""Typed metric-matrix facade over existing bounded native calculators."""

from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Iterable, Literal, Protocol, SupportsIndex

import numpy as np

if TYPE_CHECKING:
    from . import euclidean_distance_matrix_checked as _euclidean
    from . import sparse_tanimoto_matrix_checked as _tanimoto
else:
    from .chematic import euclidean_distance_matrix_checked as _euclidean
    from .chematic import sparse_tanimoto_matrix_checked as _tanimoto

_MAX_ROWS = 1_024
_MAX_IDS_PER_ROW = 65_536
_MAX_TOTAL_IDS = 1_000_000
_MAX_WIDTH = 2**32 - 1


class MetricBitVector(Protocol):
    """Explicit/sparse binary source, with width and integer on-bit positions."""

    def GetNumBits(self) -> SupportsIndex: ...

    def GetOnBits(self) -> Iterable[SupportsIndex]: ...


def _index(value: SupportsIndex) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("bit widths and IDs must be integers, excluding bool")
    return operator.index(value)


def _rows(sources: Iterable[MetricBitVector]) -> list[list[int]]:
    if isinstance(sources, (str, bytes, bytearray)):
        raise TypeError("bitVectList requires binary vector sources, not text")
    result: list[list[int]] = []
    common_width = None
    total = 0
    for source in sources:
        if len(result) >= _MAX_ROWS:
            raise ValueError("bitVectList exceeds 1024 rows")
        width_method = getattr(source, "GetNumBits", None)
        bits_method = getattr(source, "GetOnBits", None)
        if not callable(width_method) or not callable(bits_method):
            raise TypeError("binary vectors require GetNumBits() and GetOnBits()")
        width = _index(width_method())
        if not 0 <= width <= _MAX_WIDTH:
            raise ValueError("bit width must be in 0..4294967295")
        if common_width is not None and width != common_width:
            raise ValueError("binary vectors must have equal declared bit widths")
        common_width = width
        ids: list[int] = []
        positions = bits_method()
        if isinstance(positions, (str, bytes, bytearray)):
            raise TypeError("GetOnBits() must return integer positions, not text")
        # source stays alive until its lazy method-produced iterator is consumed.
        for value in positions:
            if len(ids) >= _MAX_IDS_PER_ROW:
                raise ValueError("binary vector exceeds 65536 raw bit IDs")
            if total >= _MAX_TOTAL_IDS:
                raise ValueError("bitVectList exceeds 1000000 raw bit IDs")
            bit = _index(value)
            if bit < 0 or bit > _MAX_WIDTH:
                raise OverflowError("bit ID is outside unsigned 32-bit range")
            if bit >= width:
                raise ValueError("on-bit ID is outside the declared bit width")
            ids.append(bit)
            total += 1
        result.append(ids)
    return result


def _binary(
    bit_vectors: Iterable[MetricBitVector], metric: Literal["distance", "similarity"]
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    rows = _rows(bit_vectors)
    # Staged owned positions reach native code only after all source calls pass.
    values = _tanimoto(rows, metric=metric)
    # The native sparse-set API defines empty/empty similarity as 1. This named
    # RDKit facade deliberately uses RDKit's 0 similarity / 1 distance instead.
    empty_value = 1.0 if metric == "distance" else 0.0
    for i, row in enumerate(rows):
        if not row:
            for j in range(i):
                if not rows[j]:
                    values[i * (i - 1) // 2 + j] = empty_value
    return np.array(values, dtype=np.float64)


def GetEuclideanDistMat(
    descripMat: Iterable[Iterable[float]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Return condensed Euclidean distances as an independent float64 array.

    descripMat accepts rectangular numeric rows, NumPy arrays/views and nested
    iterables. Existing native checks/numerics apply: 1024 rows, 256 features,
    50 million pair-feature operations, finite inputs and representable norms.
    Empty/singleton input returns an empty array; zero features give zeros.
    Stable scaled norms retain representable tiny/large distances. No scaling,
    descriptor generation or molecule distance calculation is performed.
    """
    return np.array(_euclidean(descripMat), dtype=np.float64)


def GetTanimotoDistMat(
    bitVectList: Iterable[MetricBitVector],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Return condensed binary Tanimoto distances, empty/empty distance 1.

    Explicit/sparse sources provide GetNumBits()/GetOnBits(); raw index rows
    and count vectors are not accepted. Widths must match in 0..4294967295.
    Maxima: 1024 rows, 65536 raw IDs per row, one million total IDs, 50 million
    native pair-bit work units. Duplicate IDs are idempotent but count toward
    raw limits. No dense vector expansion occurs. Bool widths/IDs are refused.
    Wrong types raise TypeError, invalid widths/limits ValueError, out-of-u32
    IDs OverflowError. Source exceptions propagate; AttributeError during
    method lookup is treated as a missing method. Iterator calls have no timeout.
    Empty/singleton input returns an empty independent float64 array.
    """
    return _binary(bitVectList, "distance")


def GetTanimotoSimMat(
    bitVectList: Iterable[MetricBitVector],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Return condensed binary similarities, with the same checked input profile.

    Uses intersection/union, except empty/empty similarity is 0 as in RDKit.
    Inputs are copied before calculation, never retained/mutated. The fresh
    float64 array orders pairs (1,0), (2,0), (2,1), ... without the diagonal.
    See GetTanimotoDistMat for types, bounds and errors.
    """
    return _binary(bitVectList, "similarity")


class rdMetricMatrixCalc:
    """Bounded descriptor/binary metric namespace; no RDKit runtime dependency."""

    GetEuclideanDistMat = staticmethod(GetEuclideanDistMat)
    GetTanimotoDistMat = staticmethod(GetTanimotoDistMat)
    GetTanimotoSimMat = staticmethod(GetTanimotoSimMat)
