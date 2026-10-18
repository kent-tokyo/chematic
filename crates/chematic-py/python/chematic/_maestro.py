"""Checked Maestro geometry output for a concrete, fully explicit-H profile."""
from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from numbers import Real
import re
from typing import TYPE_CHECKING, List, Optional, Sequence, Tuple, TypedDict, cast

import numpy as np

if TYPE_CHECKING:
    from . import Mol, from_mol_block, write_cdjson_graph
else:
    from .chematic import Mol, from_mol_block, write_cdjson_graph

_MAX_ATOMS = 4096
_MAX_BONDS = 8192
_MAX_OUTPUT = 2 << 20
_MAX_INPUT = 16 << 20
_REAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
_INTEGER = re.compile(r"[+-]?[0-9]+\Z")
_NUMBERS = dict(H=1, Li=3, Be=4, C=6, N=7, O=8, F=9, Na=11, Mg=12,
                P=15, S=16, Cl=17, K=19, Ca=20, Br=35, I=53)
_VALENCES = {
    ("H", 0): (1,), ("H", 1): (0,), ("H", -1): (0,), ("C", 0): (4,),
    ("N", 0): (3,), ("N", 1): (4,), ("N", -1): (2,),
    ("O", 0): (2,), ("O", 1): (3,), ("O", -1): (1,),
    ("P", 0): (3, 5), ("S", 0): (2, 4, 6),
    **{(s, 0): (1,) for s in ("F", "Cl", "Br", "I")},
    **{(s, -1): (0,) for s in ("F", "Cl", "Br", "I")},
    **{(s, 1): (0,) for s in ("Li", "Na", "K")},
    **{(s, 2): (0,) for s in ("Be", "Mg", "Ca")},
}


class _AtomBase(TypedDict):
    l: str
    c: int
    h: int


class _Atom(_AtomBase, total=False):
    m: int


class _Bond(TypedDict):
    b: int
    e: int
    o: int


class _Graph(TypedDict):
    a: List[_Atom]
    b: List[_Bond]


@dataclass(frozen=True)
class MaeExport:
    """One ASCII f_m_ct block with explicit geometry and declared omissions."""

    text: str
    title: str
    atom_count: int
    bond_count: int
    source_text: Optional[str] = None
    omitted_features: Tuple[str, ...] = (
        "source_records", "sd_properties", "atom_and_bond_metadata",
        "residue_information", "force_field_types",
    )


def _sequence(value: object, name: str) -> None:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be an indexed sequence")


def _coordinates(values: Sequence[Sequence[float]], count: int) -> List[List[float]]:
    _sequence(values, "coordinates")
    if len(values) != count:
        raise ValueError("coordinates must have one 3D row per atom")
    copied: List[List[float]] = []
    for row in values:
        _sequence(row, "coordinate row")
        if len(row) != 3:
            raise ValueError("coordinate rows must have exactly three values")
        point: List[float] = []
        for value in row:
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
                raise TypeError("coordinates must be real numbers, excluding bool")
            try:
                number = float(value)
            except OverflowError as error:
                raise ValueError("coordinate exceeds finite range") from error
            if not math.isfinite(number) or abs(number) > 1e8:
                raise ValueError("coordinates must be finite and within +/-1e8 angstrom")
            if number != 0 and abs(number) < 1e-100:
                raise ValueError("nonzero coordinates below 1e-100 are unsupported")
            point.append(number)
        copied.append(point)
    return copied


def write_mae_checked(
    mol: Mol, coordinates: Sequence[Sequence[float]], *, title: Optional[str] = None,
    allow_metadata_loss: bool = False,
) -> MaeExport:
    """Write a single Maestro structure with explicit atoms, bonds and angstrom XYZ.

    Requires metadata-loss consent. Supports fully explicit H, selected ordinary
    closed-shell valences/ions, plain 1/2/3 bonds, no aromatic/query/isotope/map or
    declared stereo. Reuses the checked native CDJSON graph projection. Does not
    generate atoms/coordinates, assign force-field types or write files. Readers
    may infer stereochemistry from geometry; no stereo-label parity is promised.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic Mol")
    if not isinstance(allow_metadata_loss, bool):
        raise TypeError("allow_metadata_loss must be bool")
    if not allow_metadata_loss:
        raise ValueError("Maestro graph export requires explicit metadata-loss consent")
    if not 1 <= mol.atom_count <= _MAX_ATOMS:
        raise ValueError("Maestro atom count must be 1..4096")
    name = mol.GetName() if title is None else title
    if not isinstance(name, str):
        raise TypeError("title must be str or None")
    if len(name) > 256 or any(not " " <= c <= "~" for c in name):
        raise ValueError("title must be at most 256 printable ASCII characters")
    points = _coordinates(coordinates, mol.atom_count)
    smiles = mol.smiles
    if any(c in smiles for c in ("@", "/", "\\")):
        raise ValueError("declared stereo is outside the Maestro output profile")
    graph = cast(List[_Graph], json.loads(write_cdjson_graph(mol, points))["m"])[0]
    atoms, bonds = graph["a"], graph["b"]
    if len(bonds) > _MAX_BONDS:
        raise ValueError("Maestro supports at most 8192 bonds")
    valence = [0] * len(atoms)
    pairs = set()
    for bond in bonds:
        pair = tuple(sorted((bond["b"], bond["e"])))
        if pair in pairs:
            raise ValueError("duplicate bonds are unsupported")
        pairs.add(pair)
        valence[bond["b"]] += bond["o"]
        valence[bond["e"]] += bond["o"]
    for atom, drawn in zip(atoms, valence):
        if "m" in atom or atom["h"] != 0:
            raise ValueError("Maestro requires fully explicit H and no isotope labels")
        if drawn not in _VALENCES.get((atom["l"], atom["c"]), ()):
            raise ValueError("element/charge/valence is outside the closed-shell profile")
    escaped = name.replace("\\", "\\\\").replace('"', '\\"')
    lines = ["{", " s_m_m2io_version", " :::", " 2.0.0", "}", "f_m_ct {",
             " s_m_title", " :::", f' "{escaped}"', f" m_atom[{len(atoms)}] {{",
             " r_m_x_coord", " r_m_y_coord", " r_m_z_coord", " i_m_atomic_number",
             " i_m_formal_charge", " :::"]
    for index, (atom, point) in enumerate(zip(atoms, points), 1):
        xyz = " ".join(format(value, ".17g") for value in point)
        lines.append(f" {index} {xyz} {_NUMBERS[atom['l']]} {atom['c']}")
    lines += [" :::", " }", f" m_bond[{len(bonds)}] {{", " i_m_from", " i_m_to",
              " i_m_order", " :::"]
    for index, bond in enumerate(bonds, 1):
        lines.append(f" {index} {bond['b'] + 1} {bond['e'] + 1} {bond['o']}")
    lines += [" :::", " }", "}", ""]
    text = "\n".join(lines)
    if len(text) > _MAX_OUTPUT:
        raise ValueError("Maestro output exceeds 2 MiB")
    return MaeExport(text, name, len(atoms), len(bonds))


def _v2000_point(row: str) -> List[float]:
    point = []
    for start in (0, 10, 20):
        token = row[start:start+10].strip()
        if not _REAL.fullmatch(token):
            raise ValueError("invalid V2000 coordinate token")
        value = float(token)
        if value == 0 and any(c in "123456789" for c in token.lower().split("e")[0]):
            raise ValueError("nonzero V2000 coordinate underflows binary64")
        point.append(value)
    return point


def _charge_records(rows: Sequence[str], count: int) -> None:
    seen = set()
    for row in rows:
        fields = row.split()
        try:
            size = _v2000_integer(fields[2])
            if (not row.startswith("M  CHG") or fields[:2] != ["M", "CHG"]
                    or not 1 <= size <= 8 or len(fields) != 3+2*size):
                raise ValueError("unsupported ancillary record")
            for offset in range(size):
                atom = _v2000_integer(fields[3+2*offset])
                charge = _v2000_integer(fields[4+2*offset])
                if not 1 <= atom <= count or atom in seen or not -3 <= charge <= 3:
                    raise ValueError("invalid/duplicate charge reference or charge")
                seen.add(atom)
        except (ValueError, IndexError) as error:
            raise ValueError("only checked, unique M CHG records are supported") from error


def _v2000_integer(token: str) -> int:
    token = token.strip()
    if not token:
        return 0
    if not _INTEGER.fullmatch(token):
        raise ValueError("invalid V2000 integer token")
    return int(token)


def v2000_to_mae_checked(text: str, *, allow_metadata_loss: bool = False) -> MaeExport:
    """Convert one V2000 MOL/SD record, retaining XYZ and original text.

    Restricts raw stereo/query/property records before native graph validation.
    SD values are opaque omitted metadata; multiple records/V3000 refuse. No
    invalid record is skipped. This is not a general SDF schema validator.
    """
    if not isinstance(text, str):
        raise TypeError("text must be str")
    if len(text) > _MAX_INPUT or len(text.encode("utf-8")) > _MAX_INPUT:
        raise ValueError("V2000 input exceeds 16 MiB")
    if any(c not in "\n\r\t" and not " " <= c <= "~" for c in text):
        raise ValueError("V2000 requires ASCII text")
    rows = text.splitlines()
    if any(len(row) > 4096 for row in rows) or len(rows) > 20000:
        raise ValueError("V2000 line/record limit exceeded")
    if len(rows) < 5 or not rows[3].rstrip().endswith("V2000"):
        raise ValueError("one V2000 MOL/SD record is required")
    try:
        count, bonds = _v2000_integer(rows[3][:3]), _v2000_integer(rows[3][3:6])
        if not 1 <= count <= 999 or not 0 <= bonds <= 999:
            raise ValueError("invalid V2000 counts")
        if len(rows) < 5 + count + bonds or _v2000_integer(rows[3][12:15]) != 0:
            raise ValueError("truncated or chiral V2000 record")
        points = [_v2000_point(row) for row in rows[4:4+count]]
        for row in rows[4:4+count]:
            if (_v2000_integer(row[34:36])
                    or _v2000_integer(row[36:39]) not in (0, 1, 2, 3, 5, 6, 7)):
                raise ValueError("atom isotope/radical flags are unsupported")
            if (any(_v2000_integer(value) for value in row[39:48].split())
                    or _v2000_integer(row[48:51]) not in (0, 15)
                    or any(_v2000_integer(value) for value in row[51:].split())):
                raise ValueError("atom stereo/query/valence/map flags are unsupported")
        for row in rows[4+count:4+count+bonds]:
            a, b = _v2000_integer(row[:3]), _v2000_integer(row[3:6])
            if a == b or not 1 <= a <= count or not 1 <= b <= count:
                raise ValueError("invalid bond reference")
            if (_v2000_integer(row[6:9]) not in (1, 2, 3)
                    or any(_v2000_integer(v) for v in row[9:].split())):
                raise ValueError("bond order/stereo/query flags are unsupported")
        end = rows.index("M  END", 4 + count + bonds)
    except (ValueError, IndexError) as error:
        raise ValueError("invalid or unsupported V2000 structure") from error
    _charge_records(rows[4+count+bonds:end], count)
    tail = rows[end+1:]
    if "$$$$" in tail:
        separator = tail.index("$$$$")
        if any(row.strip() for row in tail[separator+1:]):
            raise ValueError("multiple SD records or data after delimiter")
        tail = tail[:separator]
    # Opaque SD fields have a header followed by values terminated by a blank.
    in_value = False
    for row in tail:
        if not row.strip():
            in_value = False
        elif not in_value:
            if not row.startswith(">") or "<" not in row or ">" not in row[1:]:
                raise ValueError("unsupported text after M END")
            in_value = True
    mol = from_mol_block("\n".join(rows[:end+1]) + "\n")
    result = write_mae_checked(mol, points, title=rows[0], allow_metadata_loss=allow_metadata_loss)
    if result.atom_count != count or result.bond_count != bonds:
        raise ValueError("native projection changed source atom/bond counts")
    return replace(result, source_text=text)
