"""Typed ranked-record adapter for the existing native enrichment calculator."""

from __future__ import annotations

import math
import numbers
import operator
from typing import TYPE_CHECKING, Iterable, Protocol, SupportsIndex

import numpy as np

if TYPE_CHECKING:
    from . import enrichment_factors as _enrichment_factors
else:
    from .chematic import enrichment_factors as _enrichment_factors

_MAX_RECORDS = 100_000
_MAX_FRACTIONS = 4_096
_MAX_COLUMNS = 1_024


class RankedScoreRow(Protocol):
    """Heterogeneous indexed score record; only the active column is read."""

    def __len__(self) -> int: ...

    def __getitem__(self, index: int, /) -> object: ...


def _column(value: SupportsIndex) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("col must be an integer, excluding bool")
    column = operator.index(value)
    if not -_MAX_COLUMNS <= column < _MAX_COLUMNS:
        raise ValueError("col must be in -1024..1023")
    return column


def _fractions(values: Iterable[float]) -> list[float]:
    result: list[float] = []
    for value in values:
        if len(result) >= _MAX_FRACTIONS:
            raise ValueError("fractions exceed the 4096-entry limit")
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
            raise TypeError("fractions require real numeric values, excluding bool")
        fraction = float(value)
        if not math.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError("fractions must be finite and in (0, 1]")
        result.append(fraction)
    if not result:
        raise ValueError("fraction list is empty")
    return result


def _active(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if not isinstance(value, numbers.Real):
        raise TypeError("active labels must be bool or real numeric 0/1")
    if value != 0 and value != 1:
        raise ValueError("active labels must be binary 0/1, not weights or NaN")
    return bool(value == 1)


def CalcEnrichment(
    scores: Iterable[RankedScoreRow],
    col: SupportsIndex,
    fractions: Iterable[float],
) -> list[float]:
    """Calculate EF for caller-ranked records via the bounded native API.

    Args:
        scores: Best-first iterable of indexed records, at most 100000 rows
            and 1..1024 columns per row. Only the selected column is inspected.
        col: Activity column, including Python-style negative indexing.
            Labels must be bool or real numeric 0/1; truthy objects are refused.
        fractions: Iterable of 1..4096 finite real fractions in (0, 1].
            Order and duplicates are preserved. Bool fractions are refused.

    Returns:
        One EF per fraction, with cutoff ceil(N*f), including the full-library
        endpoint. Inputs are not sorted or retained, and output is independent.

    Raises:
        TypeError: Wrong iterable/row/column/label/fraction types.
        ValueError: Empty input, nonbinary labels, bad fractions, missing column,
            size limits or undefined EF when no active records exist.

    Column/fractions are validated before score rows, and all labels are copied
    before native calculation. Caller iteration/indexing exceptions propagate.
    This is an explicit binary-label profile, not RDKit's arbitrary truthiness.
    Ranking, score magnitudes, ties and statistical uncertainty are caller-owned.
    """
    column = _column(col)
    selected_fractions = _fractions(fractions)
    if isinstance(scores, (str, bytes, bytearray)):
        raise TypeError("scores must be an iterable of indexed records")
    labels: list[bool] = []
    for row in scores:
        if len(labels) >= _MAX_RECORDS:
            raise ValueError("scores exceed the 100000-record limit")
        if isinstance(row, (str, bytes, bytearray)):
            raise TypeError("score rows must be nontext indexed records")
        width = len(row)
        if not 1 <= width <= _MAX_COLUMNS:
            raise ValueError("score rows must contain 1..1024 columns")
        if not -width <= column < width:
            raise ValueError("col is outside the score row")
        labels.append(_active(row[column]))
    return _enrichment_factors(labels, selected_fractions)


class Scoring:
    """Bounded Scoring namespace; only CalcEnrichment is implemented."""

    CalcEnrichment = staticmethod(CalcEnrichment)
