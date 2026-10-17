"""Bounded entropy and Pearson chi-square ranking of binary fingerprints."""

from __future__ import annotations

from enum import IntEnum
import operator
from typing import Any, Iterable, Protocol, SupportsIndex

import numpy as np

MAX_RANKER_BITS = 65_536
MAX_RANKER_CLASSES = 256
MAX_RANKER_ENTRIES = 1_048_576
MAX_RANKER_RECORDS = 1_048_576


class BinaryBitVector(Protocol):
    """Structural interface shared by explicit and sparse binary bit vectors."""

    def GetNumBits(self) -> SupportsIndex: ...

    def GetOnBits(self) -> Iterable[SupportsIndex]: ...


class InfoType(IntEnum):
    ENTROPY = 1
    BIASENTROPY = 2
    CHISQUARE = 3
    BIASCHISQUARE = 4


def _integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer, excluding bool")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc


def _on_bits(value: BinaryBitVector, width: int) -> list[int]:
    get_width = getattr(value, "GetNumBits", None)
    get_on = getattr(value, "GetOnBits", None)
    if not callable(get_width) or not callable(get_on):
        raise TypeError("bitVect requires callable GetNumBits and GetOnBits")
    if _integer(get_width(), "bit-vector width") != width:
        raise ValueError("bit-vector width must match nBits")
    try:
        values = iter(get_on())
    except TypeError as exc:
        raise TypeError("GetOnBits must return an iterable of integer bit IDs") from exc
    seen: set[int] = set()
    result: list[int] = []
    for value in values:
        if len(result) >= width:
            raise ValueError("GetOnBits exceeds the declared width")
        bit = _integer(value, "on-bit ID")
        if not 0 <= bit < width or bit in seen:
            raise ValueError("on-bit IDs must be unique and within the declared width")
        seen.add(bit)
        result.append(bit)
    return result


def _entropy(
    counts: np.ndarray[Any, np.dtype[np.float64]]
) -> np.ndarray[Any, np.dtype[np.float64]]:
    totals = counts.sum(axis=1, keepdims=True)
    probabilities = np.divide(
        counts, totals, out=np.zeros_like(counts), where=totals > 0
    )
    logarithms = np.log2(
        probabilities, out=np.zeros_like(probabilities), where=probabilities > 0
    )
    result: np.ndarray[Any, np.dtype[np.float64]] = np.negative(
        np.sum(probabilities * logarithms, axis=1)
    )
    return result


def _chi_square(
    on: np.ndarray[Any, np.dtype[np.float64]],
    totals: np.ndarray[Any, np.dtype[np.float64]],
    records: int,
) -> np.ndarray[Any, np.dtype[np.float64]]:
    """Pearson statistic; empty rows/classes contribute zero, without Yates."""
    scores: np.ndarray[Any, np.dtype[np.float64]] = np.zeros(
        on.shape[0], dtype=np.float64
    )
    for observed in (on, totals[None, :] - on):
        expected = observed.sum(axis=1, keepdims=True) * (totals / records)
        differences = observed - expected
        terms = np.divide(
            differences * differences,
            expected,
            out=np.zeros_like(expected),
            where=expected > 0,
        )
        scores += terms.sum(axis=1)
    return scores


class InfoBitRanker:
    """Checked ENTROPY/CHISQUARE binary fingerprint ranker.

    Accepts structural explicit/sparse binary vectors. Counts are committed
    only after all bit IDs and the label have passed validation. No chemical
    fingerprint algorithm or dataset ordering is supplied by this class.
    """

    def __init__(
        self,
        nBits: SupportsIndex,
        nClasses: SupportsIndex,
        infoType: SupportsIndex = InfoType.ENTROPY,
    ) -> None:
        width = _integer(nBits, "nBits")
        classes = _integer(nClasses, "nClasses")
        metric = _integer(infoType, "infoType")
        if metric not in (InfoType.ENTROPY, InfoType.CHISQUARE):
            raise ValueError("only InfoType.ENTROPY and CHISQUARE are supported")
        if not 1 <= width <= MAX_RANKER_BITS or not 1 <= classes <= MAX_RANKER_CLASSES:
            raise ValueError("nBits must be 1..65536 and nClasses must be 1..256")
        if width * classes > MAX_RANKER_ENTRIES:
            raise ValueError("nBits*nClasses must not exceed 1048576 count entries")
        self._width = width
        self._classes = classes
        self._metric = metric
        self._records = 0
        self._class_counts: np.ndarray[Any, np.dtype[np.uint64]] = np.zeros(
            classes, dtype=np.uint64
        )
        self._on_counts: np.ndarray[Any, np.dtype[np.uint64]] = np.zeros(
            (width, classes), dtype=np.uint64
        )

    def AccumulateVotes(self, bitVect: BinaryBitVector, label: SupportsIndex) -> None:
        """Count one labelled fingerprint without retaining the vector.

        Args:
            bitVect: Explicit/sparse binary vector with the configured width.
            label: Integer class ID in [0, nClasses).

        Raises:
            TypeError: Unsupported vector interface or noninteger values.
            ValueError: Invalid label, width, duplicate/out-of-range bit ID,
                or record limit. Invalid input leaves all counts unchanged.
        """
        category = _integer(label, "label")
        if not 0 <= category < self._classes:
            raise ValueError("label must be in [0, nClasses)")
        if self._records >= MAX_RANKER_RECORDS:
            raise ValueError("ranker record limit exceeded")
        bits = _on_bits(bitVect, self._width)
        self._on_counts[bits, category] += 1
        self._class_counts[category] += 1
        self._records += 1

    def GetTopN(self, num: SupportsIndex) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return rows [bit ID, score, on-count for each class].

        ENTROPY scores are information gain in bits; CHISQUARE scores are the
        uncorrected Pearson statistic of the 2-by-nClasses on/off table, not
        p-values. Empty classes/constant features contribute zero. Scores rank
        descending; exact ties use ascending bit ID, unlike RDKit's heap.
        At least one record is required. The returned float64 array is a copy.
        """
        count = _integer(num, "num")
        if not 0 <= count <= self._width:
            raise ValueError("num must be in [0, nBits]")
        if self._records == 0:
            raise ValueError("information ranking requires at least one record")
        output: np.ndarray[Any, np.dtype[np.float64]] = np.empty(
            (count, self._classes + 2), dtype=np.float64
        )
        if count == 0:
            return output
        totals = self._class_counts.astype(np.float64)
        on = self._on_counts.astype(np.float64)
        if self._metric == InfoType.CHISQUARE:
            gains = _chi_square(on, totals, self._records)
        else:
            off = totals[None, :] - on
            parent = float(_entropy(totals[None, :])[0])
            probability_on = on.sum(axis=1) / self._records
            gains = np.asarray(
                parent
                - probability_on * _entropy(on)
                - (1.0 - probability_on) * _entropy(off)
            )
            gains = np.clip(gains, 0.0, parent)
        order: np.ndarray[Any, np.dtype[np.intp]] = np.lexsort(
            (np.arange(self._width), -gains)
        )
        order = order[:count]
        output[:, 0] = order
        output[:, 1] = gains[order]
        output[:, 2:] = on[order]
        return output

    def SetBiasList(self, classList: Iterable[int]) -> None:
        """Biasing is outside this ENTROPY/CHISQUARE profile."""
        raise NotImplementedError("biased metrics are unsupported")

    def SetMaskBits(self, maskBits: Iterable[int]) -> None:
        """Feature masking is outside this profile."""
        raise NotImplementedError("bit masks are unsupported")


class rdInfoTheory:
    """Typed namespace exposing the bounded entropy/chi-square ranker."""

    InfoType = InfoType
    InfoBitRanker = InfoBitRanker
