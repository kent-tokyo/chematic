"""Bounded static gray-style XBS input, independently based on the XBS manual."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re

from ._xbs import XbsBondRule, XbsSpecies, _label, _number, write_xbs

_MAX_BYTES = 4 << 20
_MAX_LINES = 20_000
_MAX_LINE_BYTES = 4096
_REAL = re.compile(
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z", re.ASCII
)


def _fail(line: int, message: str) -> ValueError:
    return ValueError(f"XBS line {line}: {message}")


def _decode(text: str | bytes) -> str:
    if not isinstance(text, (str, bytes)):
        raise TypeError("XBS input must be str or immutable ASCII bytes")
    if len(text) > _MAX_BYTES:
        raise ValueError("XBS input exceeds 4 MiB")
    try:
        source = text.decode("ascii") if isinstance(text, bytes) else text
        source.encode("ascii")
    except UnicodeError as error:
        raise ValueError("XBS input must be ASCII") from error
    return source


def _real(token: str) -> float:
    if not _REAL.fullmatch(token):
        raise ValueError("expected a decimal real number")
    value = float(token)
    if not math.isfinite(value):
        raise ValueError("nonfinite numeric field")
    if value == 0.0 and any(c in "123456789" for c in token.lower().split("e")[0]):
        raise ValueError("nonzero numeric field underflows binary64")
    return value


@dataclass
class _Rows:
    labels: list[str] = field(default_factory=list)
    positions: list[tuple[float, float, float]] = field(default_factory=list)
    species: list[XbsSpecies] = field(default_factory=list)
    bonds: list[XbsBondRule] = field(default_factory=list)
    atom_lines: list[int] = field(default_factory=list)
    species_lines: list[int] = field(default_factory=list)
    bond_lines: list[int] = field(default_factory=list)
    definitions: dict[str, int] = field(default_factory=dict)
    pairs: set[tuple[str, ...]] = field(default_factory=set)

    def add(self, fields: list[str], line: int) -> None:
        command = fields[0]
        if command == "atom":
            if len(fields) != 5:
                raise ValueError("atom requires exactly species x y z")
            if len(self.labels) == 10_000:
                raise ValueError("atom count exceeds 10000")
            label = _label(fields[1])
            xyz = tuple(_number(_real(s), "coordinate", -1e6, 1e6) for s in fields[2:])
            self.labels.append(label)
            self.positions.append((xyz[0], xyz[1], xyz[2]))
            self.atom_lines.append(line)
        elif command == "spec":
            if len(fields) != 4:
                raise ValueError(
                    "spec requires exactly name radius gray; RGB/named colors unsupported"
                )
            if len(self.species) == 128:
                raise ValueError("species count exceeds 128")
            style = XbsSpecies(fields[1], _real(fields[2]), _real(fields[3]))
            if style.label in self.definitions:
                raise ValueError(f"duplicate species style {style.label}")
            self.definitions[style.label] = line
            self.species.append(style)
            self.species_lines.append(line)
        elif command == "bonds":
            if len(fields) != 7:
                raise ValueError(
                    "bonds requires exactly name1 name2 min max radius gray"
                )
            if len(self.bonds) == 1024:
                raise ValueError("bond-rule count exceeds 1024")
            rule = XbsBondRule(fields[1], fields[2], *(_real(s) for s in fields[3:]))
            pair = tuple(sorted((rule.species1, rule.species2)))
            if pair in self.pairs:
                raise ValueError("duplicate undirected species-pair rule")
            self.pairs.add(pair)
            self.bonds.append(rule)
            self.bond_lines.append(line)
        else:
            raise ValueError(
                f"unsupported record {command!r}; expected atom, spec or bonds"
            )

    def validate_references(self) -> None:
        if not self.labels:
            raise _fail(1, "scene requires at least one atom")
        for label, line in zip(self.labels, self.atom_lines):
            if label not in self.definitions:
                raise _fail(line, f"missing style for species {label}")
        for rule, line in zip(self.bonds, self.bond_lines):
            for label in (rule.species1, rule.species2):
                if label not in self.definitions:
                    raise _fail(line, f"bond rule refers to undefined species {label}")


def _parse(source: str) -> _Rows:
    normalized = source.replace("\r\n", "\n")
    if (
        normalized.count("\n") + int(bool(normalized) and not normalized.endswith("\n"))
        > _MAX_LINES
    ):
        raise ValueError("XBS input exceeds 20000 physical lines")
    rows = _Rows()
    for number, line in enumerate(normalized.split("\n"), 1):
        if len(line) > _MAX_LINE_BYTES:
            raise _fail(number, "physical line exceeds 4096 bytes")
        if any((ord(c) < 32 and c != "\t") or ord(c) == 127 for c in line):
            raise _fail(number, "unsupported control character or line ending")
        stripped = line.strip(" \t")
        if not stripped or stripped.startswith("*"):
            continue
        try:
            rows.add(stripped.split(), number)
        except ValueError as error:
            raise _fail(number, str(error)) from error
    rows.validate_references()
    return rows


@dataclass(frozen=True, init=False, repr=False)
class XbsScene:
    """Validated immutable static XBS scene; construction parses serialized text.

    Atom/species/rule order, source text and 1-based physical line numbers are
    retained. Labels are drawing species; distance rules are not chemical bonds.
    No element lookup, units, topology inference or viewer execution is implied.
    """

    source_text: str
    labels: tuple[str, ...]
    positions: tuple[tuple[float, float, float], ...]
    species: tuple[XbsSpecies, ...]
    bonds: tuple[XbsBondRule, ...]
    atom_lines: tuple[int, ...]
    species_lines: tuple[int, ...]
    bond_lines: tuple[int, ...]

    def __init__(self, text: str | bytes) -> None:
        source = _decode(text)
        rows = _parse(source)
        object.__setattr__(self, "source_text", source)
        object.__setattr__(self, "labels", tuple(rows.labels))
        object.__setattr__(self, "positions", tuple(rows.positions))
        object.__setattr__(self, "species", tuple(rows.species))
        object.__setattr__(self, "bonds", tuple(rows.bonds))
        object.__setattr__(self, "atom_lines", tuple(rows.atom_lines))
        object.__setattr__(self, "species_lines", tuple(rows.species_lines))
        object.__setattr__(self, "bond_lines", tuple(rows.bond_lines))

    def __repr__(self) -> str:
        return f"XbsScene(atoms={len(self.labels)}, species={len(self.species)}, bond_rules={len(self.bonds)})"

    def to_xbs(self, *, preserve_source: bool = True) -> str:
        """Return exact decoded source, or normalized writer output when False."""
        if not isinstance(preserve_source, bool):
            raise TypeError("preserve_source must be bool")
        if preserve_source:
            return self.source_text
        return write_xbs(
            self.labels, self.positions, species=self.species, bonds=self.bonds
        )


def read_xbs(text: str | bytes) -> XbsScene:
    """Read the checked static gray-style XBS profile, without filesystem I/O.

    Accepts ASCII str/bytes up to 4 MiB, 20000 physical lines, 4096 bytes/line;
    1..10000 atoms, <=128 styles and <=1024 undirected rules. LF/CRLF, blank
    lines and whole-line '*' comments are accepted. Records may interleave or
    refer to later styles; duplicates and unresolved references are errors.
    Numbers are finite decimal binary64, with coordinates in +/-1e6; a nonzero
    literal that underflows to zero is rejected. Style/rule bounds match write_xbs.
    RGB/named colors, wildcards, viewer commands, movies, periodic metadata and
    the other 'bs' dialect are outside this profile and are not silently ignored.
    Wrong input types raise TypeError; malformed/unsupported input ValueError,
    including a physical line number for record errors. No partial scene returns.
    """
    return XbsScene(text)
