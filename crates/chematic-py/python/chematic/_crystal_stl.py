"""Bounded ASCII STL atom spheres with explicit millimeter display scaling."""
from __future__ import annotations

import io
import math
from dataclasses import dataclass
from numbers import Real
from typing import TYPE_CHECKING, Literal, Sequence, Tuple

from ._crystal_obj import Vector, _OMITTED, _prepare_spheres, _sphere_positions

if TYPE_CHECKING:
    from . import PeriodicStructure
else:
    from .chematic import PeriodicStructure

_MAX_TRIANGLES = 50000
_MAX_OUTPUT = 16 * 1024 * 1024
_STL_OMISSIONS = _OMITTED + (
    "site_and_sphere_component_identity",
    "surface_union",
    "printability_validation",
)


@dataclass(frozen=True)
class StlMeshExport:
    """ASCII geometry and explicit scale; no chemical or printable-solid claim."""

    text: str
    sphere_count: int
    triangle_count: int
    millimeters_per_angstrom: float
    units: Literal["millimeter"] = "millimeter"
    format: Literal["ascii_stl"] = "ascii_stl"
    omitted_features: Tuple[str, ...] = _STL_OMISSIONS


def _scale(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError("millimeters_per_angstrom must be a real number, not bool")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError("millimeters_per_angstrom is outside the finite range") from error
    if not math.isfinite(result) or result <= 0:
        raise ValueError("millimeters_per_angstrom must be finite and positive")
    return result


def _facet_normal(a: Vector, b: Vector, c: Vector) -> Vector:
    # Scale coordinates before subtraction/cross products to avoid overflow.
    scale = max(abs(x) for point in (a, b, c) for x in point)
    if scale == 0:
        raise ValueError("STL triangle collapses at the requested scale")
    u = tuple(b[i] / scale - a[i] / scale for i in range(3))
    v = tuple(c[i] / scale - a[i] / scale for i in range(3))
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    length = math.hypot(*cross)
    if not math.isfinite(length) or length == 0:
        raise ValueError("STL triangle has an unrepresentable normal")
    return (cross[0] / length, cross[1] / length, cross[2] / length)


def write_crystal_stl_checked(
    structure: PeriodicStructure,
    radii: Sequence[float],
    *,
    millimeters_per_angstrom: float,
    subdivisions: int = 2,
) -> StlMeshExport:
    """Render ordered periodic sites as finite atom spheres in ASCII STL.

    Args:
        structure: Native structure with 1..1000 ordered, fully occupied sites.
            Known unexpanded CIF symmetry is refused. No wrapping is performed.
        radii: One finite positive caller-selected radius per site, in angstroms.
        millimeters_per_angstrom: Required finite positive display scale, applied
            equally to Cartesian positions and radii. Configure STL import as mm;
            STL has no standardized unit field. No default scale is inferred.
        subdivisions: Integer 0..4, default 2; octahedral sphere subdivision.

    Returns:
        Frozen report with ASCII text, sphere/triangle counts, scale, millimeter
        unit convention and explicit omissions. Normals follow outward winding.

    Raises:
        TypeError: Invalid structure, radius, subdivision or scale type.
        ValueError: Invalid geometry/site/symmetry, nonfinite or collapsed scaled
            geometry, or limits (50000 triangles, 16 MiB, shared OBJ vertex cap).

    Retains site order; does not infer bonds, expand symmetry, combine surfaces,
    repair overlaps, or guarantee a printable solid. Chemical identity, lattice,
    labels, colors and per-sphere identity are omitted. The function performs no
    I/O, leaves its input unchanged and requires no mesh library at runtime.
    """
    factor = _scale(millimeters_per_angstrom)
    centers, sizes, _, unit_vertices, faces = _prepare_spheres(
        structure, radii, subdivisions
    )
    triangle_count = len(centers) * len(faces)
    if triangle_count > _MAX_TRIANGLES:
        raise ValueError("STL triangle budget exceeds 50000")
    output, size = io.StringIO(), 0

    def emit(line: str) -> None:
        nonlocal size
        size += len(line) + 1
        if size > _MAX_OUTPUT:
            raise ValueError("STL output exceeds 16 MiB")
        output.write(line + "\n")

    emit("solid chematic_atom_spheres_mm")
    for center, radius in zip(centers, sizes):
        scaled = (center[0] * factor, center[1] * factor, center[2] * factor)
        scaled_radius = radius * factor
        if (
            not all(math.isfinite(x) for x in scaled)
            or not math.isfinite(scaled_radius)
            or scaled_radius <= 0
        ):
            raise ValueError("scaled STL center/radius overflows or underflows")
        vertices = _sphere_positions(scaled, scaled_radius, unit_vertices)
        for i, j, k in faces:
            a, b, c = vertices[i], vertices[j], vertices[k]
            normal = _facet_normal(a, b, c)
            emit("  facet normal " + " ".join(format(x, ".17g") for x in normal))
            emit("    outer loop")
            for point in (a, b, c):
                emit("      vertex " + " ".join(format(x, ".17g") for x in point))
            emit("    endloop")
            emit("  endfacet")
    emit("endsolid chematic_atom_spheres_mm")
    return StlMeshExport(output.getvalue(), len(centers), triangle_count, factor)
