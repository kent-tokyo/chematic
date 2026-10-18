"""Bounded entropy and Pearson ranking, with optional class bias."""

from __future__ import annotations

from enum import IntEnum
import operator
from typing import Any, Iterable, Protocol, SupportsIndex

import numpy as np

from ._bit_correlation import BitCorrMatGenerator as BitCorrMatGenerator

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
    """Checked binary fingerprint ranker with optional class-biased scores.

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
        if metric not in tuple(InfoType):
            raise ValueError("infoType must be a supported InfoType value")
        if not 1 <= width <= MAX_RANKER_BITS or not 1 <= classes <= MAX_RANKER_CLASSES:
            raise ValueError("nBits must be 1..65536 and nClasses must be 1..256")
        if width * classes > MAX_RANKER_ENTRIES:
            raise ValueError("nBits*nClasses must not exceed 1048576 count entries")
        self._width = width
        self._classes = classes
        self._metric = metric
        self._bias: tuple[int, ...] = ()
        self._records = 0
        self._mask_bits: frozenset[int] | None = None
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
        if self._mask_bits is not None:
            bits = [bit for bit in bits if bit in self._mask_bits]
        self._on_counts[bits, category] += 1
        self._class_counts[category] += 1
        self._records += 1

    def GetTopN(self, num: SupportsIndex) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return rows [bit ID, score, on-count for each class].

        ENTROPY scores are information gain in bits; CHISQUARE scores are the
        uncorrected Pearson statistic of the 2-by-nClasses on/off table, not
        p-values. BIAS modes zero scores failing SetBiasList preferences.
        Empty classes/constant features contribute zero. Scores rank
        descending; exact ties use ascending bit ID, unlike RDKit's heap.
        Only currently masked-in bits are eligible; num cannot exceed their
        unique count. At least one record is required. The returned float64
        array is a copy.
        """
        count = _integer(num, "num")
        available = self._width if self._mask_bits is None else len(self._mask_bits)
        if not 0 <= count <= available:
            raise ValueError("num must be in [0, number of eligible bits]")
        if self._records == 0:
            raise ValueError("information ranking requires at least one record")
        output: np.ndarray[Any, np.dtype[np.float64]] = np.empty(
            (count, self._classes + 2), dtype=np.float64
        )
        if count == 0:
            return output
        totals = self._class_counts.astype(np.float64)
        on = self._on_counts.astype(np.float64)
        if self._metric in (InfoType.CHISQUARE, InfoType.BIASCHISQUARE):
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
        if (
            self._metric in (InfoType.BIASENTROPY, InfoType.BIASCHISQUARE)
            and 0 < len(self._bias) < self._classes
        ):
            fractions = np.divide(on, totals, out=np.zeros_like(on), where=totals > 0)
            others = [c for c in range(self._classes) if c not in self._bias]
            accepted = fractions[:, self._bias].max(axis=1) >= fractions[:, others].max(
                axis=1
            )
            gains = np.where(accepted, gains, 0.0)
        candidates = (
            np.arange(self._width)
            if self._mask_bits is None
            else np.asarray(sorted(self._mask_bits), dtype=np.intp)
        )
        order = candidates[np.lexsort((candidates, -gains[candidates]))][:count]
        output[:, 0] = order
        output[:, 1] = gains[order]
        output[:, 2:] = on[order]
        return output

    def SetBiasList(self, classList: Iterable[SupportsIndex]) -> None:
        """Replace the preferred classes after validating the complete iterable.

        Args:
            classList: Unique integer IDs in [0, nClasses). An empty
                list accepts every bit. Fewer than nClasses IDs are allowed.
                The list is copied, not retained.

        Raises:
            TypeError: Noniterable input or noninteger IDs (including bool).
            ValueError: Duplicate/out-of-range IDs or nClasses or more IDs.

        In biased modes, a bit retains its base score if any preferred class
        has an on-fraction at least as high as every nonpreferred class;
        otherwise its score becomes zero. Empty classes have fraction zero.
        Counts remain intact; zero-scored bits still participate in ranking.
        Unbiased modes ignore the preference. Replacing it after votes
        reranks the existing data.
        Invalid input leaves the previous preference and counts unchanged.
        """
        selected: set[int] = set()
        for value in classList:
            if len(selected) >= self._classes - 1:
                raise ValueError("classList must contain fewer than nClasses IDs")
            category = _integer(value, "bias class ID")
            if not 0 <= category < self._classes or category in selected:
                raise ValueError("bias class IDs must be unique and in [0, nClasses)")
            selected.add(category)
        self._bias = tuple(sorted(selected))

    def SetMaskBits(self, maskBits: Iterable[SupportsIndex]) -> None:
        """Replace the eligible feature set for future votes and ranking.

        Args:
            maskBits: At most nBits integer IDs in [0, nBits). Duplicates
                collapse to a set. Empty selects no features; range(nBits)
                restores all features. Text, bytes and bool IDs are rejected.

        Raises:
            TypeError: The source is not an integer iterable.
            ValueError: An ID is out of range or the raw item limit is exceeded.

        The complete iterable is validated before replacing the previous mask.
        Existing counts and class totals remain unchanged. Later votes update
        on-counts only for selected bits, but always update class totals. A
        newly selected bit does not regain on-votes skipped by an earlier mask.
        Exceptions from user iteration propagate; caller side effects are not
        rolled back. This method does not retain the source iterable.
        """
        if isinstance(maskBits, (str, bytes, bytearray)):
            raise TypeError("maskBits must be an iterable of integer bit IDs")
        selected: set[int] = set()
        for position, value in enumerate(maskBits):
            if position >= self._width:
                raise ValueError("maskBits raw item count must not exceed nBits")
            bit = _integer(value, "mask bit ID")
            if not 0 <= bit < self._width:
                raise ValueError("mask bit ID must be in [0, nBits)")
            selected.add(bit)
        self._mask_bits = frozenset(selected)


class rdInfoTheory:
    """Typed bounded bit ranking and pair co-occurrence namespace."""

    InfoType = InfoType
    InfoBitRanker = InfoBitRanker
    BitCorrMatGenerator = BitCorrMatGenerator
