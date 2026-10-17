"""Copy-preserving conformer alignment with atomic caller-owned RMS output."""

from __future__ import annotations

from collections.abc import Sequence as RuntimeSequence
from typing import Any, Optional, Sequence, Union

import numpy as np

from ._point_alignment import Points, Weights, _points, _weights, rdAlignment

Conformers = Union[Sequence[Points], np.ndarray[Any, np.dtype[Any]]]
MAX_CONFORMERS = 1024
MAX_CONFORMER_POINTS = 1_000_000


def _indices(value: Optional[Sequence[int]], size: int, name: str) -> list[int]:
    if value is None:
        return list(range(size))
    if not isinstance(value, RuntimeSequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a sequence of integer indices")
    if len(value) > size:
        raise ValueError(f"{name} cannot contain more than {size} indices")
    result = list(value)
    if any(type(i) is not int for i in result):
        raise TypeError(f"{name} requires integer indices, excluding booleans")
    if any(i < 0 or i >= size for i in result):
        raise ValueError(f"{name} indices must be in 0..{size - 1}")
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must not contain duplicate indices")
    return result or list(range(size))


def align_conformers_checked(
    conformers: Conformers,
    *,
    atom_ids: Optional[Sequence[int]] = None,
    conformer_ids: Optional[Sequence[int]] = None,
    weights: Optional[Weights] = None,
    rms_list: Optional[list[float]] = None,
) -> np.ndarray[Any, np.dtype[np.float64]]:
    """Align selected coordinate frames to the first selected frame.

    Args:
        conformers: Finite real coordinates with shape (C, N, 3), in common
            units and atom order. Arrays and sequences of point clouds work.
        atom_ids: Unique atom indices used for fitting; None/empty means all.
        conformer_ids: Unique frame indices in fit/output order; None/empty
            means all. The first selected frame is the fixed reference.
        weights: Positive weights for the selected atoms; None/empty means
            equal weights. RMS is sqrt(weighted SSD / selected atom count),
            matching RDKit rather than dividing by the sum of weights.
        rms_list: Ordinary list to append one float per selected probe frame.
            Existing entries are preserved. It changes only after every fit
            succeeds; unsupported containers or nonfinite entries raise.

    Returns:
        A new float64 (C, N, 3) array. Each fit transforms every atom of its
        frame, including atoms outside atom_ids. Reference and unselected
        frames retain their input coordinates. Input frames are never edited.

    Raises:
        TypeError: Unsupported containers, booleans or nonnumeric values.
        ValueError: Empty/ragged/nonfinite coordinates, repeated/out-of-range
            indices, invalid weights, limits, or an unrepresentable fit.

    This coordinate API does not mutate an RDKit Mol, generate conformers,
    match atoms or permit reflections. Limits: 1024 frames, 65536 atoms per
    frame, one million points in total, coordinate magnitude <= 1e100.
    """
    if rms_list is not None:
        if type(rms_list) is not list:
            raise TypeError("rms_list must be an ordinary list of floats or None")
        if any(type(v) is not float for v in rms_list):
            raise TypeError("rms_list existing entries must be floats")
        if not np.isfinite(rms_list).all():
            raise ValueError("rms_list existing entries must be finite")
    if isinstance(conformers, np.ndarray):
        if type(conformers) is not np.ndarray or conformers.dtype.kind not in "iuf":
            raise TypeError("conformers requires an ordinary real numeric ndarray")
        if conformers.ndim != 3 or conformers.shape[2] != 3:
            raise ValueError("conformers must have shape (C, N, 3)")
        count = conformers.shape[0]
    elif isinstance(conformers, RuntimeSequence) and not isinstance(
        conformers, (str, bytes)
    ):
        count = len(conformers)
    else:
        raise TypeError("conformers must be a coordinate sequence or numeric ndarray")
    if not 1 <= count <= MAX_CONFORMERS:
        raise ValueError(f"conformers requires 1..{MAX_CONFORMERS} frames")
    # Check sizes before allocating the output or copying any frame.
    n_atoms = 0
    for frame in conformers:
        if not isinstance(frame, (RuntimeSequence, np.ndarray)) or isinstance(
            frame, (str, bytes)
        ):
            raise TypeError("each conformer must be a point sequence or ndarray")
        if isinstance(frame, np.ndarray) and frame.ndim == 0:
            raise ValueError("each conformer requires a point dimension")
        if n_atoms == 0:
            n_atoms = len(frame)
            if not 1 <= n_atoms <= 65_536:
                raise ValueError("each conformer requires 1..65536 atoms")
        elif len(frame) != n_atoms:
            raise ValueError("all conformers must have the same atom count")
    if count * n_atoms > MAX_CONFORMER_POINTS:
        raise ValueError("conformers exceeds one million coordinate points")
    atoms = _indices(atom_ids, n_atoms, "atom_ids")
    frames = _indices(conformer_ids, count, "conformer_ids")
    w = _weights(weights, len(atoms))
    relative = w / w.max()
    if (relative / relative.sum() == 0).any():
        raise ValueError("weight ratios exceed float64 resolution")
    result: np.ndarray[Any, np.dtype[np.float64]] = np.empty(
        (count, n_atoms, 3), dtype=np.float64
    )
    for i, frame in enumerate(conformers):
        result[i] = _points(frame, f"conformers[{i}]")
    reference = result[frames[0], atoms].copy()
    rms_values: list[float] = []
    for i in frames[1:]:
        ssd, transform = rdAlignment.GetAlignmentTransform(
            reference, result[i, atoms], w
        )
        result[i] = result[i] @ transform[:3, :3].T + transform[:3, 3]
        if not np.isfinite(result[i]).all():
            raise ValueError(
                "aligned coordinates are not representable as finite float64"
            )
        rms_values.append(float(np.sqrt(ssd / len(atoms))))
    if rms_list is not None:
        rms_list.extend(rms_values)
    return result
