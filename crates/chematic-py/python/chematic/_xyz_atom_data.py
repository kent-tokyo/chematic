"""Checked headerless XYZ + per-atom real data (Open Babel issue #1703).

Independent data-only adapter. No topology, charge, vector units or tensor
layout is inferred; named extxyz output uses the caller's property name.
"""

from __future__ import annotations

from dataclasses import dataclass
import io
import json
import math
import re
from typing import Literal, Tuple, Union

_MAX_BYTES = 16 << 20
_MAX_ATOMS = 50_000
_MAX_LINE_BYTES = 4096
_REAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
_PROPERTY = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,31}\Z", re.ASCII)
_CONTROL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f\r\u0085\u2028\u2029]")
_ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs "
    "Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db "
    "Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)


def _decode(text: Union[str, bytes]) -> str:
    if not isinstance(text, (str, bytes)):
        raise TypeError("XYZ atom data must be str or immutable UTF-8 bytes")
    if len(text) > _MAX_BYTES:
        raise ValueError("XYZ atom data exceeds 16 MiB")
    try:
        source = text.decode("utf-8") if isinstance(text, bytes) else text
        if len(source.encode("utf-8")) > _MAX_BYTES:
            raise ValueError("XYZ atom data exceeds 16 MiB")
    except UnicodeError as error:
        raise ValueError("XYZ atom data must be valid UTF-8") from error
    return source


def _line(stream: io.StringIO, name: str) -> str:
    line = stream.readline(_MAX_LINE_BYTES + 2)
    if not line:
        raise ValueError(f"missing {name}")
    line = line.removesuffix("\n")
    if len(line.encode("utf-8")) > _MAX_LINE_BYTES:
        raise ValueError(f"{name} exceeds 4096 bytes")
    return line


def _number(token: str, line: int) -> float:
    if not _REAL.fullmatch(token):
        raise ValueError(f"line {line}: expected decimal real data")
    value = float(token)
    if not math.isfinite(value):
        raise ValueError(f"line {line}: nonfinite real data")
    if value == 0.0 and any(c in "123456789" for c in token.lower().split("e")[0]):
        raise ValueError(f"line {line}: nonzero real data underflows binary64")
    return value


@dataclass(frozen=True, init=False, repr=False)
class XyzAtomDataFrame:
    """Immutable one-frame XYZ record with 1, 3 or 9 real values per atom.

    Construct from serialized text; default components=3 reads xyz+vec.
    Element symbols are canonical H..Og (no isotope/dummy/atomic-number aliases).
    Positions and payload values retain input order and units without inference.
    source_text preserves the exact decoded input, including CRLF/spacing.
    """

    elements: Tuple[str, ...]
    positions: Tuple[Tuple[float, float, float], ...]
    values: Tuple[Tuple[float, ...], ...]
    components: int
    comment: str
    source_text: str

    def __init__(
        self, text: Union[str, bytes], *, components: Literal[1, 3, 9] = 3
    ) -> None:
        if isinstance(components, bool) or not isinstance(components, int):
            raise TypeError("components must be an integer: 1, 3 or 9")
        if components not in (1, 3, 9):
            raise ValueError("components must be 1, 3 or 9")
        source = _decode(text)
        normalized = source.replace("\r\n", "\n")
        if _CONTROL.search(normalized):
            raise ValueError("unsupported control character or line ending")
        stream = io.StringIO(normalized)
        count_line = _line(stream, "atom count").strip(" \t")
        if not re.fullmatch(r"[0-9]{1,6}", count_line):
            raise ValueError("atom count must be decimal digits")
        count = int(count_line)
        if not 1 <= count <= _MAX_ATOMS:
            raise ValueError("atom count must be 1..50000")
        comment = _line(stream, "comment line")
        if re.search(r"(?:^|\s)(?:Properties|Lattice)\s*=", comment):
            raise ValueError("named/periodic extxyz headers require from_extxyz")
        elements = []
        positions = []
        values = []
        columns = 4 + components
        for i in range(count):
            fields = _line(stream, f"atom row {i + 1}").split(None, columns)
            if len(fields) != columns:
                raise ValueError(f"line {i + 3}: expected exactly {columns} columns")
            if fields[0] not in _ELEMENTS:
                raise ValueError(f"line {i + 3}: unsupported element symbol")
            numbers = tuple(_number(token, i + 3) for token in fields[1:])
            elements.append(fields[0])
            positions.append((numbers[0], numbers[1], numbers[2]))
            values.append(numbers[3:])
        if stream.read().strip(" \t\n"):
            raise ValueError("extra records/data after the single XYZ frame")
        object.__setattr__(self, "elements", tuple(elements))
        object.__setattr__(self, "positions", tuple(positions))
        object.__setattr__(self, "values", tuple(values))
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "comment", comment)
        object.__setattr__(self, "source_text", source)

    def __repr__(self) -> str:
        return f"XyzAtomDataFrame(atoms={len(self.elements)}, components={self.components})"

    def _rows(self) -> str:
        lines = []
        for symbol, xyz, data in zip(self.elements, self.positions, self.values):
            lines.append(symbol + " " + " ".join(format(x, ".17g") for x in xyz + data))
        return "\n".join(lines) + "\n"

    def to_xyz(self) -> str:
        """Return normalized XYZ+data with LF and 17-digit binary64 values."""
        return f"{len(self.elements)}\n{self.comment}\n" + self._rows()

    def to_extxyz(self, property_name: str = "atom_data") -> str:
        """Name the payload column explicitly for existing extxyz readers.

        Nine components remain flat in input order; no tensor reshape is done.
        The original comment is JSON-tagged quoted metadata under `comment`.
        ASE decodes the marker; readers retaining raw strings need JSON decoding.
        """
        if not isinstance(property_name, str):
            raise TypeError("property_name must be str")
        if not _PROPERTY.fullmatch(property_name) or property_name in {
            "species",
            "pos",
            "positions",
            "numbers",
            "Z",
        }:
            raise ValueError("invalid/reserved extxyz property name")
        metadata = "_JSON " + json.dumps(self.comment, ensure_ascii=False)
        escaped = metadata.replace("\\", "\\\\").replace('"', '\\"')
        header = (
            f"Properties=species:S:1:pos:R:3:{property_name}:R:{self.components} "
            f'comment="{escaped}"'
        )
        return f"{len(self.elements)}\n{header}\n" + self._rows()


def parse_xyz_atom_data(
    text: Union[str, bytes], *, components: Literal[1, 3, 9] = 3
) -> XyzAtomDataFrame:
    """Parse one bounded XYZ+scalar/vector/flat-tensor frame without guessing."""
    return XyzAtomDataFrame(text, components=components)
