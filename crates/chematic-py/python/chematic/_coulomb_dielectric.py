"""Typed local-solvation Coulomb model with RDKit 2026.09.1 conventions.

Independent Python implementation; reference conventions and the empirical
neighbor table come from Release_2026_09_1 MIFDescriptors.cpp (BSD-3-Clause).
See THIRD_PARTY_NOTICES.md. No molecule/conformer overload or charge inference.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Tuple

from ._coulomb_field import (
    Coulomb,
    MifPositionSequence,
    MifScalarSequence,
    _MAX_COORDINATE,
    _PREFACTOR,
    _length,
    _real,
)

_MAX_DIELECTRIC_SITES = 1024


def _solvation_depth(neighbors: int) -> float:
    if neighbors <= 6:
        return 0.0
    if neighbors >= 12:
        return 4.0
    return (0.4, 0.9, 1.4, 1.9, 2.6)[neighbors - 7]


def _site_depths(
    positions: Tuple[Tuple[float, float, float], ...]
) -> Tuple[float, ...]:
    counts = [0] * len(positions)
    for i, (x, y, z) in enumerate(positions):
        for j in range(i + 1, len(positions)):
            ax, ay, az = positions[j]
            dx, dy, dz = x - ax, y - ay, z - az
            # Explicit-charge reference constructor compares distance to 4 A;
            # unlike probe counting, it takes sqrt before the comparison.
            if math.sqrt(dx * dx + dy * dy + dz * dz) < 4.0:
                counts[i] += 1
                counts[j] += 1
    return tuple(_solvation_depth(count) for count in counts)


@dataclass(frozen=True, init=False, repr=False)
class CoulombDielectric:
    """Copied charges/positions with a bounded local-solvation dielectric model.

    Args:
        charges: Indexed charges in elementary-charge units; at most 1024 sites.
        positions: Indexed rows of exactly three Cartesian angstrom coordinates.
        probeCharge: Probe charge in elementary-charge units, in +/-1e6.
        absVal: Return negative absolute value of the total interaction.
        softcoreParam: Zero for hard cutoff, or 1e-6..1e200 A^2 for soft core.
        cutoff: Minimum hard-core distance, 0.001..1e100 angstroms.
        epsilon: Solvent relative permittivity, finite in 1..1e6 (default 80).
        xi: Solute relative permittivity, finite in 1..1e6 (default 4).

    Other scalar/input rules are shared with Coulomb. Construction copies data
    and counts all site pairs within 4 A, taking O(N^2) time. Evaluation uses
    O(N) local scratch, including all sites when counting probe neighbors,
    even zero-charge sites or those excluded by the energy threshold.
    """

    _base: Coulomb
    _depths: Tuple[float, ...]
    _contrast: float
    _coefficient: float

    def __init__(
        self,
        charges: MifScalarSequence,
        positions: MifPositionSequence,
        probeCharge: float = 1.0,
        absVal: bool = False,
        softcoreParam: float = 0.0,
        cutoff: float = 1.0,
        epsilon: float = 80.0,
        xi: float = 4.0,
    ) -> None:
        solvent = _real(epsilon, "epsilon", 1.0, 1e6)
        solute = _real(xi, "xi", 1.0, 1e6)
        # Enforce the quadratic-work cap before indexing/copying any site.
        _length(charges, "charges", _MAX_DIELECTRIC_SITES)
        _length(positions, "positions", _MAX_DIELECTRIC_SITES)
        probe = _real(probeCharge, "probeCharge", -1e6, 1e6)
        base = Coulomb(charges, positions, probe, absVal, softcoreParam, cutoff)
        object.__setattr__(self, "_base", base)
        object.__setattr__(self, "_depths", _site_depths(base._positions))
        object.__setattr__(self, "_contrast", (solute - solvent) / (solute + solvent))
        object.__setattr__(self, "_coefficient", _PREFACTOR * (1.0 / solute) * probe)

    @property
    def site_count(self) -> int:
        return self._base.site_count

    def __repr__(self) -> str:
        return f"CoulombDielectric(sites={self.site_count})"

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """Return kJ/mol; include sites with squared distance strictly < threshold.

        The threshold is finite in 0..1e200 A^2. Probe coordinates are finite
        in +/-1e100 A. Rounded threshold and 4 A neighbor boundaries can differ
        from compiled C++ arithmetic. No scratch data is retained on the field.
        """
        px = _real(x, "x", -_MAX_COORDINATE, _MAX_COORDINATE)
        py = _real(y, "y", -_MAX_COORDINATE, _MAX_COORDINATE)
        pz = _real(z, "z", -_MAX_COORDINATE, _MAX_COORDINATE)
        limit = _real(threshold, "threshold", 0.0, 1e200)
        distance2 = []
        for ax, ay, az in self._base._positions:
            dx, dy, dz = px - ax, py - ay, pz - az
            distance2.append(dx * dx + dy * dy + dz * dz)
        probe_depth = _solvation_depth(sum(d2 < 16.0 for d2 in distance2))
        total = 0.0
        for charge, d2, depth in zip(self._base._charges, distance2, self._depths):
            if d2 < limit:
                core2 = (
                    self._base._alpha + d2
                    if self._base._alpha
                    else max(d2, self._base._cutoff2)
                )
                inverse = 1.0 / math.sqrt(core2)
                screened = self._contrast / math.sqrt(core2 + 4.0 * depth * probe_depth)
                total += charge * (inverse + screened)
        energy = total * self._coefficient
        if not math.isfinite(energy):
            raise ValueError(
                "Coulomb dielectric energy is not representable as a finite float"
            )
        return -abs(energy) if self._base._abs else energy
