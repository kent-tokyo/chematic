"""Checked point-cloud alignment for the rdAlignment compatibility namespace."""

from __future__ import annotations

from collections.abc import Sequence as RuntimeSequence
from numbers import Real
from typing import Any, Optional, Protocol, Sequence, Union

import numpy as np

MAX_ALIGNMENT_POINTS = 65_536
MAX_ALIGNMENT_MAGNITUDE = 1e100


class Point3DLike(Protocol):
    """Read-only Cartesian point interface, including RDKit Point3D objects."""

    @property
    def x(self) -> float: ...

    @property
    def y(self) -> float: ...

    @property
    def z(self) -> float: ...


Points = Union[
    Sequence[Union[Sequence[float], Point3DLike]], np.ndarray[Any, np.dtype[Any]]
]
Weights = Union[Sequence[float], np.ndarray[Any, np.dtype[Any]]]


def _real(value: Any, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must contain real numbers, excluding booleans")
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} is outside the supported numeric range") from exc
    if not np.isfinite(result) or abs(result) > MAX_ALIGNMENT_MAGNITUDE:
        raise ValueError(f"{name} must be finite with magnitude <= 1e100")
    return result


def _numeric_array(
    value: np.ndarray[Any, np.dtype[Any]], name: str
) -> np.ndarray[Any, np.dtype[np.float64]]:
    if type(value) is not np.ndarray or value.dtype.kind not in "iuf":
        raise TypeError(f"{name} requires an ordinary real numeric ndarray")
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.array(value, dtype=np.float64, copy=True)
    if (
        not np.isfinite(result).all()
        or (np.abs(result) > MAX_ALIGNMENT_MAGNITUDE).any()
    ):
        raise ValueError(f"{name} must be finite with magnitude <= 1e100")
    return result


def _points(value: Points, name: str) -> np.ndarray[Any, np.dtype[np.float64]]:
    if isinstance(value, np.ndarray):
        if value.ndim != 2 or value.shape[1] != 3:
            raise ValueError(f"{name} must have shape (N, 3)")
        n = value.shape[0]
    elif isinstance(value, RuntimeSequence) and not isinstance(value, (str, bytes)):
        n = len(value)
    else:
        raise TypeError(f"{name} must be a point sequence or numeric ndarray")
    if not 1 <= n <= MAX_ALIGNMENT_POINTS:
        raise ValueError(f"{name} requires 1..{MAX_ALIGNMENT_POINTS} points")
    if isinstance(value, np.ndarray):
        return _numeric_array(value, name)
    result: np.ndarray[Any, np.dtype[np.float64]] = np.empty((n, 3), dtype=np.float64)
    for i in range(n):
        point = value[i]
        if isinstance(point, np.ndarray):
            if point.shape != (3,):
                raise ValueError(f"{name}[{i}] must have three coordinates")
            result[i] = _numeric_array(point, name)
        elif isinstance(point, RuntimeSequence) and not isinstance(point, (str, bytes)):
            if len(point) != 3:
                raise ValueError(f"{name}[{i}] must have three coordinates")
            result[i] = [_real(v, name) for v in point]
        else:
            try:
                xyz = (point.x, point.y, point.z)
            except AttributeError as exc:
                raise TypeError(
                    f"{name}[{i}] requires three coordinates or x/y/z"
                ) from exc
            result[i] = [_real(v, name) for v in xyz]
    return result


def _weights(value: Optional[Weights], n: int) -> np.ndarray[Any, np.dtype[np.float64]]:
    if value is None:
        return np.ones(n, dtype=np.float64)
    if isinstance(value, np.ndarray):
        if value.ndim != 1:
            raise ValueError("weights must be one-dimensional")
        count = value.size
    elif isinstance(value, RuntimeSequence) and not isinstance(value, (str, bytes)):
        count = len(value)
    else:
        raise TypeError("weights must be a numeric sequence or ndarray")
    if count not in (0, n):
        raise ValueError("weights must be empty or have one value per point pair")
    result: np.ndarray[Any, np.dtype[np.float64]]
    if isinstance(value, np.ndarray):
        result = _numeric_array(value, "weights")
    else:
        result = np.array([_real(v, "weights") for v in value], dtype=np.float64)
    if count == 0:
        return np.ones(n, dtype=np.float64)
    if np.any(result <= 0):
        raise ValueError("weights must be strictly positive")
    return result


class rdAlignment:
    """Typed subset of rdkit.Numerics.rdAlignment, using a NumPy SVD solver."""

    @staticmethod
    def GetAlignmentTransform(
        refPoints: Points,
        probePoints: Points,
        weights: Optional[Weights] = None,
        reflect: bool = False,
        maxIterations: int = 50,
    ) -> tuple[float, np.ndarray[Any, np.dtype[np.float64]]]:
        """Return weighted SSD and a 4×4 probe-to-reference affine transform.

        Args:
            refPoints: Nonempty N×3 real array, coordinate rows, or x/y/z points.
            probePoints: Paired points in the same order, units and shape.
            weights: Positive pair weights; None or empty means equal weights.
            reflect: Require an improper rotation when True, proper otherwise.
            maxIterations: Only the default 50 is accepted. The SVD backend
                has no iterative-eigensolver budget; other values raise.

        Returns:
            Weighted sum of squared deviations (not RMSD), and a float64 4×4
            matrix acting on homogeneous column vectors. Inputs are copied.

        Raises:
            TypeError: Unsupported containers, scalar types or option types.
            ValueError: Empty/mismatched points, invalid numeric data, limits,
                unsupported iteration budget, or a failed numerical solve.
        """
        if type(reflect) is not bool or type(maxIterations) is not int:
            raise TypeError("reflect must be bool and maxIterations must be int")
        if maxIterations != 50:
            raise ValueError("only maxIterations=50 is supported by the SVD backend")
        reference = _points(refPoints, "refPoints")
        probe = _points(probePoints, "probePoints")
        if reference.shape != probe.shape:
            raise ValueError("point counts must match")
        w = _weights(weights, len(reference))
        relative = w / w.max()
        normalized = relative / relative.sum()
        if (normalized == 0).any():
            raise ValueError("weight ratios exceed float64 resolution")
        ref_center = normalized @ reference
        probe_center = normalized @ probe
        ref = reference - ref_center
        moving = probe - probe_center
        scale = max(float(np.abs(ref).max()), float(np.abs(moving).max()))
        if scale == 0:
            rotation = np.diag([-1.0, -1.0, -1.0]) if reflect else np.eye(3)
        else:
            covariance = ((moving / scale) * normalized[:, None]).T @ (ref / scale)
            try:
                u, _, vt = np.linalg.svd(covariance)
            except np.linalg.LinAlgError as exc:
                raise ValueError("point alignment SVD did not converge") from exc
            sign = -1.0 if reflect else 1.0
            correction: np.ndarray[Any, np.dtype[np.float64]] = np.eye(3)
            correction[2, 2] = sign * (1.0 if np.linalg.det(vt.T @ u.T) >= 0 else -1.0)
            rotation = vt.T @ correction @ u.T
        transform: np.ndarray[Any, np.dtype[np.float64]] = np.eye(4, dtype=np.float64)
        transform[:3, :3] = rotation
        transform[:3, 3] = ref_center - rotation @ probe_center
        delta = moving @ rotation.T - ref
        ssd = float(np.sum(w * np.sum(delta * delta, axis=1)))
        if not np.isfinite(ssd) or not np.isfinite(transform).all():
            raise ValueError("alignment result is not representable as finite float64")
        return ssd, transform
