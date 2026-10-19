"""Checked inputs for the native default Pattern fingerprint facade."""

from __future__ import annotations

from typing import Optional, SupportsIndex

import numpy as np

from ._layered_counts import (
    Counters,
    LayeredAtomCounts as PatternAtomCounts,
    LayeredBitMask as PatternBitMask,
    _mask,
    integer,
)


def prepare(
    atoms: int,
    bonds: int,
    fp_size: SupportsIndex,
    counts: Optional[PatternAtomCounts],
    mask: Optional[PatternBitMask],
    tautomer: bool,
) -> int:
    if not isinstance(tautomer, bool):
        raise TypeError("tautomerFingerprints must be bool")
    if tautomer:
        raise ValueError("tautomerFingerprints=True is not supported")
    width = integer(fp_size, "fpSize", 1, 2048)
    if width & (width - 1):
        raise ValueError("fpSize must be a power of two in 1..2048")
    if atoms > 1024 or bonds > 4096:
        raise ValueError(
            "PatternFingerprint supports at most 1024 atoms and 4096 bonds"
        )
    # RDKit's pinned Pattern algorithm validates but does not use these two
    # arguments. Do not manufacture counts or apply a mask absent upstream.
    if mask is not None and (
        not callable(getattr(mask, "GetNumBits", None))
        or not callable(getattr(mask, "GetOnBits", None))
    ):
        raise TypeError("setOnlyBits requires callable GetNumBits and GetOnBits")
    _mask(mask, width)
    empty = (type(counts) is list and not counts) or (
        type(counts) is np.ndarray and counts.ndim == 1 and counts.size == 0
    )
    Counters.prepare(counts, 0 if empty else atoms)
    return width
