"""Checked atom-sphere OBJ geometry for native periodic structures."""

from __future__ import annotations

import io
import math
from collections.abc import Sequence as RuntimeSequence
from dataclasses import dataclass
from numbers import Real
from typing import TYPE_CHECKING, List, Literal, Sequence, Tuple

if TYPE_CHECKING:
    from . import PeriodicStructure
else:
    from .chematic import PeriodicStructure

Vector = Tuple[float, float, float]
Face = Tuple[int, int, int]
_MAX_SITES = 1000
_MAX_VERTICES = 100000
_MAX_FACES = 200000
_MAX_OUTPUT = 16 * 1024 * 1024
_OMITTED = (
    "lattice",
    "periodic_boundary_conditions",
    "bonds",
    "site_labels_and_properties",
    "occupancy_and_chemical_semantics",
    "materials_and_colors",
)


@dataclass(frozen=True)
class ObjMeshExport:
    """Finite display geometry, with one independently closed sphere per site."""

    text: str
    sphere_count: int
    vertex_count: int
    face_count: int
    units: Literal["angstrom"] = "angstrom"
    omitted_features: Tuple[str, ...] = _OMITTED


def _unit_sphere(subdivisions: int) -> Tuple[List[Vector], List[Face]]:
    vertices: List[Vector] = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, -1.0),
    ]
    faces: List[Face] = [
        (4, 0, 2),
        (4, 2, 1),
        (4, 1, 3),
        (4, 3, 0),
        (5, 2, 0),
        (5, 1, 2),
        (5, 3, 1),
        (5, 0, 3),
    ]
    for _ in range(subdivisions):
        edges = {}

        def midpoint(a: int, b: int) -> int:
            key = (min(a, b), max(a, b))
            if key not in edges:
                xyz = tuple(vertices[a][i] + vertices[b][i] for i in range(3))
                length = math.hypot(*xyz)
                edges[key] = len(vertices)
                vertices.append((xyz[0] / length, xyz[1] / length, xyz[2] / length))
            return edges[key]

        refined: List[Face] = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)))
        faces = refined
    return vertices, faces


def _checked_radii(values: Sequence[float], count: int) -> List[float]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values, RuntimeSequence
    ):
        raise TypeError("radii must be a sequence of real numbers")
    if len(values) != count:
        raise ValueError("radii must have one value per site")
    result: List[float] = []
    for index in range(count):
        value = values[index]
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError("each radius must be a real number, not bool")
        try:
            radius = float(value)
        except OverflowError as error:
            raise ValueError(
                "radius is outside the finite floating-point range"
            ) from error
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("radii must be finite and positive, in angstroms")
        result.append(radius)
    return result


def _sphere_positions(
    center: Vector, radius: float, normals: Sequence[Vector]
) -> List[Vector]:
    result = []
    for normal in normals:
        vertex = tuple(center[i] + radius * normal[i] for i in range(3))
        if not all(math.isfinite(x) for x in vertex):
            raise ValueError("mesh vertex overflow")
        # Refuse collapsed/distorted spheres when the radius is too small
        # relative to an enormous translated coordinate. 17g preserves these
        # checked floating-point positions on serialization.
        recovered = tuple((vertex[i] - center[i]) / radius for i in range(3))
        if math.dist(recovered, normal) > 1e-8:
            raise ValueError(
                "sphere radius cannot be represented at this coordinate scale"
            )
        result.append((vertex[0], vertex[1], vertex[2]))
    return result


def _prepare_spheres(
    structure: PeriodicStructure, radii: Sequence[float], subdivisions: int
) -> Tuple[List[Vector], List[float], List[str], List[Vector], List[Face]]:
    """Shared checked sites, caller radii and finite unit-sphere template."""
    if not isinstance(structure, PeriodicStructure):
        raise TypeError("structure must be a chematic PeriodicStructure")
    if type(subdivisions) is not int:
        raise TypeError("subdivisions must be an integer, not bool")
    if not 0 <= subdivisions <= 4:
        raise ValueError("subdivisions must be between 0 and 4")
    count = structure.site_count()
    if not 1 <= count <= _MAX_SITES:
        raise ValueError("OBJ mesh requires 1..1000 sites")
    vertex_count = count * (4 ** (subdivisions + 1) + 2)
    face_count = count * 8 * 4**subdivisions
    if vertex_count > _MAX_VERTICES or face_count > _MAX_FACES:
        raise ValueError("OBJ mesh vertex or face budget exceeded")
    sizes = _checked_radii(radii, count)
    status = structure.symmetry_status
    if status is not None and not status.is_complete_cell:
        raise ValueError("OBJ mesh refuses known unexpanded CIF symmetry")
    symbols = []
    for site in structure.sites:
        species = site.species
        if len(species) != 1 or species[0][1] != 1.0:
            raise ValueError("OBJ mesh needs one fully occupied species per site")
        symbol = species[0][0]
        if not symbol.isascii() or not symbol.isalpha():
            raise ValueError("OBJ mesh requires canonical element symbols")
        symbols.append(symbol)
    positions = structure.cartesian_positions().tolist()
    if len(positions) != count or len(symbols) != count:
        raise ValueError("inconsistent native site dimensions")
    centers: List[Vector] = []
    for row in positions:
        if len(row) != 3 or not all(math.isfinite(float(x)) for x in row):
            raise ValueError("Cartesian positions must have three finite components")
        centers.append((float(row[0]), float(row[1]), float(row[2])))
    normals, faces = _unit_sphere(subdivisions)
    return centers, sizes, symbols, normals, faces


def write_crystal_obj_checked(
    structure: PeriodicStructure,
    radii: Sequence[float],
    *,
    subdivisions: int = 2,
) -> ObjMeshExport:
    """Tessellate ordered periodic sites as atom spheres in Wavefront OBJ.

    Args:
        structure: Native periodic structure with 1..1000 ordered sites. Known
            unexpanded CIF symmetry and partial/mixed occupancy are refused.
        radii: One finite positive display radius per site, in angstroms. These
            are caller-selected display sizes; no chemical radii are inferred.
        subdivisions: Integer 0..4. Each step subdivides an octahedron's faces
            and projects the new vertices onto the sphere. Default 2.

    Returns:
        Frozen report containing OBJ text, sphere/vertex/triangle counts and
        omitted features. One coordinate unit means one angstrom; OBJ readers
        must be configured for that scale because OBJ has no standard unit field.

    Raises:
        TypeError: Wrong structure, radii or subdivision type.
        ValueError: Unsupported site/symmetry, invalid radius, unrepresentable
            vertex, or resource cap (100000 vertices, 200000 faces, 16 MiB).

    Uses supplied Cartesian positions without wrapping, symmetry expansion or
    lattice reduction. Each site has a named object, shared radial unit normals,
    and counterclockwise triangles. This is an atom display mesh: it omits
    periodic topology, bonds, labels, chemical properties, colors and materials.
    Spheres may overlap or be disconnected; no surface union or printable-solid
    guarantee is made. The input is unchanged and no files are written.
    """
    centers, sizes, symbols, normals, faces = _prepare_spheres(
        structure, radii, subdivisions
    )
    count = len(centers)
    vertex_count, face_count = count * len(normals), count * len(faces)
    output = io.StringIO()
    size = 0

    def emit(line: str) -> None:
        nonlocal size
        size += len(line) + 1  # Only ASCII tokens, all names internally generated.
        if size > _MAX_OUTPUT:
            raise ValueError("OBJ mesh output exceeds 16 MiB")
        output.write(line + "\n")

    emit("# chematic atom-sphere mesh; coordinate_unit angstrom")
    emit("# finite geometry; periodic boundaries and chemical semantics omitted")
    for normal in normals:
        emit("vn " + " ".join(format(x, ".17g") for x in normal))
    for index, (center, radius, symbol) in enumerate(zip(centers, sizes, symbols)):
        vertices = _sphere_positions(center, radius, normals)
        emit(f"o site_{index + 1}_{symbol}")
        emit(f"# site_index {index} radius_angstrom {radius:.17g}")
        for vertex in vertices:
            emit("v " + " ".join(format(x, ".17g") for x in vertex))
        offset = index * len(normals)
        for a, b, c in faces:
            emit("f " + " ".join(f"{offset + i + 1}//{i + 1}" for i in (a, b, c)))
    return ObjMeshExport(output.getvalue(), count, vertex_count, face_count)
