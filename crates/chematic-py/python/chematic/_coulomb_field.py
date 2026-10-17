"""Typed explicit-charge vacuum Coulomb field, with pinned RDKit MIF behavior.

Independent implementation of the Coulomb equation. Behavioral reference:
RDKit Release_2026_09_1 Code/GraphMol/MolInteractionFields/MIFDescriptors.cpp.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from numbers import Real
from typing import Protocol, Tuple, runtime_checkable

_MAX_SITES = 100_000
_MAX_COORDINATE = 1e100
_MAX_CHARGE = 1e6
# RDKit 2026.09.1's legacy conversion factor, for compatibility (kJ/mol, e, A).
_PREFACTOR = 1 / (4.0 * 3.141592 * 8.854188) * 1.602 * 1.602 * 6.02214129 * 10000


@runtime_checkable
class _Indexed(Protocol):
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> object: ...


class MifScalarSequence(Protocol):
    """Bounded indexed real scalars, including ordinary numeric NumPy vectors."""

    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> float: ...


class MifPositionSequence(Protocol):
    """Indexed rows of three real scalars, in Cartesian angstroms."""

    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> MifScalarSequence: ...


def _length(value: _Indexed, name: str, maximum: int) -> int:
    if isinstance(value, (str, bytes, bytearray, Mapping)) or not isinstance(
        value, _Indexed
    ):
        raise TypeError(f"{name} must be an indexed numeric sequence")
    size = len(value)
    if size > maximum:
        raise ValueError(f"{name} has more than {maximum} entries")
    return size


def _real(value: float, name: str, lower: float, upper: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number, excluding bool")
    try:
        number = float(value)
    except OverflowError as error:
        raise ValueError(f"{name} does not fit a finite float") from error
    if not math.isfinite(number) or not lower <= number <= upper:
        raise ValueError(f"{name} must be finite and in {lower}..{upper}")
    return number


def _snapshot(
    charges: MifScalarSequence, positions: MifPositionSequence
) -> tuple[Tuple[float, ...], Tuple[Tuple[float, float, float], ...]]:
    size = _length(charges, "charges", _MAX_SITES)
    if _length(positions, "positions", _MAX_SITES) != size:
        raise ValueError("charges and positions must have equal lengths")
    q = []
    xyz = []
    for i in range(size):
        q.append(_real(charges[i], "charge", -_MAX_CHARGE, _MAX_CHARGE))
        row = positions[i]
        if _length(row, "position row", 3) != 3:
            raise ValueError("position rows must contain exactly 3 values")
        xyz.append(
            tuple(
                _real(row[j], "coordinate", -_MAX_COORDINATE, _MAX_COORDINATE)
                for j in range(3)
            )
        )
    # Each row's length was checked before indexing/conversion.
    return tuple(q), tuple((row[0], row[1], row[2]) for row in xyz)


@dataclass(frozen=True, init=False, repr=False)
class Coulomb:
    """A copied, typed explicit-charge vacuum interaction field.

    Args:
        charges: Indexed per-site charges in elementary-charge units.
        positions: Indexed rows of exactly three Cartesian angstrom coordinates.
        probeCharge: Charge of the virtual probe, in elementary-charge units.
        absVal: Return the negative absolute value of the total interaction.
        softcoreParam: Zero for hard cutoff, or at least 1e-6 A^2 for soft core.
        cutoff: Minimum hard-core distance, 0.001 through 1e100 angstroms.

    Molecule/conformer construction, charge generation, dielectric fields and
    grid output are outside this profile. Non-sequence input raises TypeError;
    numeric bounds/length violations raise ValueError. Inputs are copied into
    immutable tuples; later source mutations do not change the field.
    """

    _charges: Tuple[float, ...]
    _positions: Tuple[Tuple[float, float, float], ...]
    _coefficient: float
    _abs: bool
    _alpha: float
    _cutoff2: float

    def __init__(
        self,
        charges: MifScalarSequence,
        positions: MifPositionSequence,
        probeCharge: float = 1.0,
        absVal: bool = False,
        softcoreParam: float = 0.0,
        cutoff: float = 1.0,
    ) -> None:
        probe = _real(probeCharge, "probeCharge", -_MAX_CHARGE, _MAX_CHARGE)
        alpha = _real(softcoreParam, "softcoreParam", 0.0, 1e200)
        floor = _real(cutoff, "cutoff", 0.001, _MAX_COORDINATE)
        if not isinstance(absVal, bool):
            raise TypeError("absVal must be bool")
        if 0.0 < alpha < 1e-6:
            raise ValueError("positive softcoreParam must be at least 1e-6 A^2")
        q, xyz = _snapshot(charges, positions)
        object.__setattr__(self, "_charges", q)
        object.__setattr__(self, "_positions", xyz)
        object.__setattr__(self, "_coefficient", _PREFACTOR * probe)
        object.__setattr__(self, "_abs", absVal)
        object.__setattr__(self, "_alpha", alpha)
        object.__setattr__(self, "_cutoff2", floor * floor)

    @property
    def site_count(self) -> int:
        return len(self._charges)

    def __repr__(self) -> str:
        return f"Coulomb(sites={self.site_count})"

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """Return kJ/mol for sites with squared distance strictly < threshold.

        `threshold` is in A^2, matching the pinned scalar RDKit call. It must be
        finite in 0..1e200. Probe coordinates obey the constructor's 1e100 bound.
        Rounded site-boundary inclusion can differ from compiled C++ arithmetic.
        """
        px = _real(x, "x", -_MAX_COORDINATE, _MAX_COORDINATE)
        py = _real(y, "y", -_MAX_COORDINATE, _MAX_COORDINATE)
        pz = _real(z, "z", -_MAX_COORDINATE, _MAX_COORDINATE)
        limit = _real(threshold, "threshold", 0.0, 1e200)
        result = 0.0
        for q, (ax, ay, az) in zip(self._charges, self._positions):
            dx, dy, dz = px - ax, py - ay, pz - az
            distance2 = dx * dx + dy * dy + dz * dz
            if distance2 < limit:
                denominator2 = (
                    self._alpha + distance2
                    if self._alpha
                    else max(distance2, self._cutoff2)
                )
                result += q * (1.0 / math.sqrt(denominator2))
        energy = result * self._coefficient
        if not math.isfinite(energy):
            raise ValueError("Coulomb energy is not representable as a finite float")
        return -abs(energy) if self._abs else energy
