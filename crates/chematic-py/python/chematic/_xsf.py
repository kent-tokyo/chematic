"""Checked static 3D-periodic XSF geometry using existing native structures.

Independent implementation of the XCrySDen format specification. Primitive
keywords carry the supplied cell; no primitive/conventional reduction is done.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    from . import PeriodicStructure
else:
    from .chematic import PeriodicStructure

_MAX_SITES = 50_000
_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs "
    "Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db "
    "Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
_NUMBERS = {symbol: i + 1 for i, symbol in enumerate(_SYMBOLS)}


def _vector(row: Sequence[float]) -> str:
    if len(row) != 3:
        raise ValueError("XSF needs exactly three Cartesian components")
    xyz = tuple(float(value) for value in row)
    if not all(math.isfinite(value) for value in xyz):
        raise ValueError("XSF lattice and Cartesian positions must be finite")
    return " ".join(format(value, ".17g") for value in xyz)


def write_xsf_checked(
    structure: PeriodicStructure, *, allow_metadata_loss: bool = False
) -> str:
    """Return static CRYSTAL/PRIMVEC/PRIMCOORD geometry in angstroms.

    Args:
        structure: Native validated periodic structure, 1..50,000 ordered sites.
        allow_metadata_loss: Must be True: XSF geometry omits site labels,
            CIF provenance and any POSCAR constraints/velocities/restart data.

    Partial/mixed occupancy and known unexpanded CIF symmetry are always
    refused, even with metadata-loss consent. The source structure is unchanged;
    site order and the supplied lattice basis/Cartesian positions are retained.
    No wrapping, symmetry search, reduction, topology, forces or grids are added.
    Native geometry getters require NumPy, as in the existing crystal API.
    """
    if not isinstance(structure, PeriodicStructure):
        raise TypeError("structure must be a chematic PeriodicStructure")
    if not isinstance(allow_metadata_loss, bool):
        raise TypeError("allow_metadata_loss must be bool")
    if not allow_metadata_loss:
        raise ValueError("XSF exports geometry only; set allow_metadata_loss=True")
    count = structure.site_count()
    if not 1 <= count <= _MAX_SITES:
        raise ValueError("XSF requires 1..50000 sites")
    status = structure.symmetry_status
    if status is not None and not status.is_complete_cell:
        raise ValueError("XSF refuses known unexpanded CIF symmetry")
    # The native wrapper is immutable; check the cap before cloning sites or
    # allocating Cartesian arrays via its existing getters.
    sites = structure.sites
    numbers = []
    for site in sites:
        species = site.species
        if len(species) != 1 or species[0][1] != 1.0:
            raise ValueError("XSF requires exactly one fully occupied species per site")
        symbol = species[0][0]
        if symbol not in _NUMBERS:
            raise ValueError("XSF needs a canonical element symbol")
        numbers.append(_NUMBERS[symbol])
    cell = structure.lattice.matrix.tolist()
    coords = structure.cartesian_positions().tolist()
    if len(cell) != 3 or len(coords) != count or len(numbers) != count:
        raise ValueError("inconsistent native geometry dimensions")
    lines = ["CRYSTAL", "PRIMVEC"]
    lines.extend(_vector(row) for row in cell)
    lines.extend(("PRIMCOORD", f"{count} 1"))
    lines.extend(f"{number} {_vector(row)}" for number, row in zip(numbers, coords))
    return "\n".join(lines) + "\n"
