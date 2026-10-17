"""Checked static XBS scene output, independently based on the XBS manual.

XBS atom/spec/bonds records are distinct from Open Babel's Ball and Stick
format with the same .bs extension. No viewer or bond perception is invoked.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math
from numbers import Real
import re

_LABEL = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,7}\Z", re.ASCII)


def _label(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("XBS species labels must be strings")
    if not _LABEL.fullmatch(value):
        raise ValueError(
            "XBS labels must be 1..8 ASCII letters/digits/underscores, starting with a letter"
        )
    return value


def _number(value: float, name: str, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number, excluding bool")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(f"{name} is not a finite float") from error
    if not math.isfinite(result) or not low <= result <= high:
        raise ValueError(f"{name} must be finite in {low}..{high}")
    return result


def _size(value: Sequence[object], name: str, maximum: int) -> int:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence, such as a list or tuple")
    size = len(value)
    if size > maximum:
        raise ValueError(f"{name} exceeds {maximum} entries")
    return size


@dataclass(frozen=True)
class XbsSpecies:
    """Explicit species drawing style: label, radius and gray (0 black, 1 white)."""

    label: str
    radius: float
    gray: float

    def __post_init__(self) -> None:
        _label(self.label)
        radius = _number(self.radius, "radius", 0.0, 1e6)
        if radius == 0.0:
            raise ValueError("radius must be positive")
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "gray", _number(self.gray, "gray", 0.0, 1.0))


@dataclass(frozen=True)
class XbsBondRule:
    """Species-pair distance drawing rule, separate from chemical bond orders."""

    species1: str
    species2: str
    min_distance: float
    max_distance: float
    radius: float
    gray: float

    def __post_init__(self) -> None:
        _label(self.species1)
        _label(self.species2)
        lower = _number(self.min_distance, "min_distance", 0.0, 1e6)
        upper = _number(self.max_distance, "max_distance", 0.0, 1e6)
        if lower >= upper:
            raise ValueError("min_distance must be less than max_distance")
        radius = _number(self.radius, "radius", 0.0, 1e6)
        if radius == 0.0:
            raise ValueError("radius must be positive")
        object.__setattr__(self, "min_distance", lower)
        object.__setattr__(self, "max_distance", upper)
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "gray", _number(self.gray, "gray", 0.0, 1.0))


def write_xbs(
    labels: Sequence[str],
    positions: Sequence[Sequence[float]],
    *,
    species: Sequence[XbsSpecies],
    bonds: Sequence[XbsBondRule] = (),
) -> str:
    """Return checked atom/spec/bonds XBS text without filesystem writes.

    All distances/radii share the caller's coordinate unit (normally angstrom).
    Atom order and binary64 values are retained with 17 significant digits.
    Requires 1..10,000 atoms, <=128 styles and <=1,024 undirected pair rules.
    Every used label needs exactly one style; extra styles are allowed. Bounds
    on coordinates are +/-1e6; styles validate their own finite bounds.
    Generators, wildcards, color names/RGB, view commands and movies are outside
    this static gray-style profile. No graph, charge or topology is inferred.
    """
    size = _size(labels, "labels", 10_000)
    if not size:
        raise ValueError("XBS output needs at least one atom")
    if _size(positions, "positions", 10_000) != size:
        raise ValueError("labels and positions must have equal lengths")
    style_count = _size(species, "species", 128)
    rule_count = _size(bonds, "bonds", 1024)
    styles = tuple(species[i] for i in range(style_count))
    rules = tuple(bonds[i] for i in range(rule_count))
    defined = set()
    lines = []
    for style in styles:
        if not isinstance(style, XbsSpecies):
            raise TypeError("species entries must be XbsSpecies")
        if style.label in defined:
            raise ValueError("duplicate species style")
        defined.add(style.label)
    for i in range(size):
        label = _label(labels[i])
        if label not in defined:
            raise ValueError(f"missing style for species {label}")
        row = positions[i]
        if _size(row, "position row", 3) != 3:
            raise ValueError("each position must have exactly three coordinates")
        xyz = [_number(row[j], "coordinate", -1e6, 1e6) for j in range(3)]
        lines.append("atom " + label + " " + " ".join(format(x, ".17g") for x in xyz))
    for style in styles:
        lines.append(f"spec {style.label} {style.radius:.17g} {style.gray:.17g}")
    pairs = set()
    for rule in rules:
        if not isinstance(rule, XbsBondRule):
            raise TypeError("bonds entries must be XbsBondRule")
        if rule.species1 not in defined or rule.species2 not in defined:
            raise ValueError("bond rule refers to an undefined species")
        pair = tuple(sorted((rule.species1, rule.species2)))
        if pair in pairs:
            raise ValueError("duplicate undirected species-pair rule")
        pairs.add(pair)
        lines.append(
            f"bonds {rule.species1} {rule.species2} {rule.min_distance:.17g} "
            f"{rule.max_distance:.17g} {rule.radius:.17g} {rule.gray:.17g}"
        )
    return "\n".join(lines) + "\n"
