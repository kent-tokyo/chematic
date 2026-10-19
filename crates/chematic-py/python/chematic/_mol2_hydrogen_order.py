"""Group explicit MOL2 H rows with their bonded parents, updating ID references."""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from typing import TYPE_CHECKING, Any, Dict, List, Tuple

if TYPE_CHECKING:
    from . import parse_mol2_record_json
else:
    from .chematic import parse_mol2_record_json

_MAX_BYTES = 16 << 20
_MAX_ATOMS = 16_384
_MAX_BONDS = 32_768
_SECTIONS = {"MOLECULE", "ATOM", "BOND", "SUBSTRUCTURE", "UNITY_ATOM_ATTR"}


@dataclass(frozen=True)
class Mol2HydrogenOrder:
    """Source/output text and exact atom correspondence; indices are zero-based.

    atom_id_map is (old ID, output ID) in source-row order. hydrogen_parents
    contains (source H ID, source parent ID) in source H-row order. Source IDs
    remain distinct from row indices. No file is written.
    """

    source_text: str
    mol2: str
    new_to_old: Tuple[int, ...]
    old_to_new: Tuple[int, ...]
    atom_id_map: Tuple[Tuple[int, int], ...]
    hydrogen_parents: Tuple[Tuple[int, int], ...]


def _fields(line: str, values: Dict[int, int]) -> str:
    spans = list(re.finditer(r"\S+", line))
    for field, value in sorted(values.items(), reverse=True):
        start, end = spans[field].span()
        line = line[:start] + str(value) + line[end:]
    return line


def _rows(lines: List[str]) -> Dict[str, List[int]]:
    result: Dict[str, List[int]] = {}
    section = ""
    for number, line in enumerate(lines):
        stripped = line.strip()
        if stripped.upper().startswith("@<TRIPOS>"):
            section = stripped[9:].upper()
            if section not in _SECTIONS or section in result:
                raise ValueError("unsupported or duplicate MOL2 section: " + section)
            result[section] = []
        elif stripped and not stripped.startswith("#"):
            if not section:
                raise ValueError("data outside a MOL2 section")
            result[section].append(number)
    if not {"MOLECULE", "ATOM"}.issubset(result):
        raise ValueError("MOLECULE and ATOM sections are required")
    return result


def _references(
    lines: List[str], rows: Dict[str, List[int]], record: Dict[str, Any]
) -> Dict[int, List[int]]:
    """Validate source-only references; return token positions to remap."""
    ids = {atom["source_id"] for atom in record["atoms"]}
    references: Dict[int, List[int]] = {}
    for line in rows.get("BOND", []):
        references[line] = [1, 2]
    substructures = set()
    for line in rows.get("SUBSTRUCTURE", []):
        fields = lines[line].split()
        if len(fields) < 3:
            raise ValueError("SUBSTRUCTURE requires ID, name and root atom ID")
        try:
            sub_id, root = int(fields[0]), int(fields[2])
        except ValueError as error:
            raise ValueError("invalid SUBSTRUCTURE ID/root") from error
        if sub_id <= 0 or sub_id in substructures or (root != 0 and root not in ids):
            raise ValueError("invalid/duplicate SUBSTRUCTURE ID or missing root atom")
        substructures.add(sub_id)
        if root != 0:
            references[line] = [2]
    if "SUBSTRUCTURE" in rows:
        for atom in record["atoms"]:
            sub_id = atom["substructure_id"]
            if sub_id is not None and sub_id not in substructures:
                raise ValueError("atom references an undeclared SUBSTRUCTURE")
    attributes = rows.get("UNITY_ATOM_ATTR", [])
    if len(attributes) % 2:
        raise ValueError("UNITY_ATOM_ATTR supports one charge per group")
    for offset in range(0, len(attributes), 2):
        header, value = attributes[offset : offset + 2]
        fields = lines[header].split()
        charge = lines[value].split()
        if (
            len(fields) != 2
            or fields[1] != "1"
            or len(charge) != 2
            or charge[0] != "charge"
        ):
            raise ValueError("only single charge UNITY_ATOM_ATTR groups are supported")
        if int(fields[0]) not in ids:
            raise ValueError("UNITY_ATOM_ATTR references a missing atom")
        references[header] = [0]
    return references


def _verify(
    before: Dict[str, Any],
    after: Dict[str, Any],
    order: List[int],
    mapping: Dict[int, int],
) -> None:
    expected_atoms = []
    for index in order:
        atom = dict(before["atoms"][index])
        atom["source_id"] = mapping[atom["source_id"]]
        expected_atoms.append(atom)
    expected_bonds = []
    for item in before["bonds"]:
        bond = dict(item)
        bond["atom_ids"] = [mapping[index] for index in bond["atom_ids"]]
        expected_bonds.append(bond)
    expected_coordinates = [before["coordinates"][i] for i in order]
    expected_charges = [before["formal_charges"][i] for i in order]
    if (
        after["atoms"] != expected_atoms
        or after["bonds"] != expected_bonds
        or after["coordinates"] != expected_coordinates
        or after["formal_charges"] != expected_charges
        or any(
            before[key] != after[key]
            for key in ("name", "molecule_type", "charge_type", "opaque_sections")
        )
    ):
        raise RuntimeError("MOL2 hydrogen ordering failed mapped readback verification")


def order_mol2_hydrogens_checked(text: str) -> Mol2HydrogenOrder:
    """Place explicit H rows immediately after their bonded non-H parent.

    Preserve non-H order and source order among each parent's H atoms. Only
    explicit H types (H or H.*) with exactly one ordinary single bond to a
    non-H atom are accepted. No H addition, distance inference or renaming.

    One source record is required: at most 16 MiB, 16384 atoms and 32768 bonds.
    Supported sections are MOLECULE/ATOM/BOND, SUBSTRUCTURE (ID roots), and
    charge-only UNITY_ATOM_ATTR. Unsupported sections/attributes are refused
    because opaque data may contain atom references.

    If row order changes, atom IDs become 1..N and references in BOND,
    SUBSTRUCTURE roots and UNITY headers are updated. Other row fields,
    coordinates, types, partial charges, bond IDs and source text are retained.
    With unchanged row order, the exact original text and IDs are returned.
    Atom rows require a common LF or CRLF terminator, including the final row.
    Blank/comment lines keep their original physical slots; they are not
    interpreted as atom annotations.

    The native source reader validates input and mapped output metadata,
    coordinates, charges and bonds. This does not validate force-field typing,
    valence or broader MOL2 chemistry. Wrong types raise TypeError; unsupported,
    malformed and oversized sources raise ValueError. Verification failures
    raise RuntimeError. No partial result or file write occurs.
    """
    if type(text) is not str:
        raise TypeError("text must be a plain str")
    if len(text) > _MAX_BYTES or len(text.encode("utf-8")) > _MAX_BYTES:
        raise ValueError("MOL2 hydrogen ordering input exceeds 16 MiB")
    if any(c in text for c in ("\x00", "\x0b", "\x0c", "\x85", "\u2028", "\u2029")):
        raise ValueError("unsupported MOL2 control/line separator")
    lines = text.splitlines(keepends=True)
    if len(lines) > 131_072 or any(
        len(line.encode("utf-8")) > 16_384 for line in lines
    ):
        raise ValueError("MOL2 physical line limits exceeded")
    rows = _rows(lines)
    if (
        not 1 <= len(rows["ATOM"]) <= _MAX_ATOMS
        or len(rows.get("BOND", [])) > _MAX_BONDS
    ):
        raise ValueError("MOL2 requires 1..16384 atoms and at most 32768 bonds")
    endings = {
        "\r\n" if lines[i].endswith("\r\n") else "\n" if lines[i].endswith("\n") else ""
        for i in rows["ATOM"]
    }
    if len(endings) != 1 or "" in endings:
        raise ValueError("ATOM rows require a common LF or CRLF terminator")
    record: Dict[str, Any] = json.loads(parse_mol2_record_json(text))
    atoms = record["atoms"]
    if len(atoms) != len(rows["ATOM"]) or len(record["bonds"]) != len(
        rows.get("BOND", [])
    ):
        raise ValueError("native record and physical row counts disagree")
    references = _references(lines, rows, record)
    indices = {atom["source_id"]: i for i, atom in enumerate(atoms)}
    hydrogens = {
        i for i, atom in enumerate(atoms) if atom["atom_type"].split(".")[0] == "H"
    }
    neighbors: List[List[Tuple[int, str]]] = [[] for _ in atoms]
    for bond in record["bonds"]:
        left, right = (indices[i] for i in bond["atom_ids"])
        neighbors[left].append((right, bond["bond_type"]))
        neighbors[right].append((left, bond["bond_type"]))
    children: Dict[int, List[int]] = {}
    parents = []
    for hydrogen in sorted(hydrogens):
        attached = neighbors[hydrogen]
        if len(attached) != 1 or attached[0][0] in hydrogens or attached[0][1] != "1":
            raise ValueError("each H requires one ordinary single bond to a non-H atom")
        parent = attached[0][0]
        children.setdefault(parent, []).append(hydrogen)
        parents.append((atoms[hydrogen]["source_id"], atoms[parent]["source_id"]))
    order = [
        index
        for parent in range(len(atoms))
        if parent not in hydrogens
        for index in [parent, *children.get(parent, [])]
    ]
    inverse = [0] * len(atoms)
    for new, old in enumerate(order):
        inverse[old] = new
    changed = order != list(range(len(atoms)))
    mapping = {
        atom["source_id"]: inverse[i] + 1 if changed else atom["source_id"]
        for i, atom in enumerate(atoms)
    }
    if changed:
        output = lines.copy()
        for destination, old in zip(rows["ATOM"], order):
            output[destination] = _fields(
                lines[rows["ATOM"][old]], {0: mapping[atoms[old]["source_id"]]}
            )
        for line, fields in references.items():
            tokens = lines[line].split()
            output[line] = _fields(
                lines[line], {field: mapping[int(tokens[field])] for field in fields}
            )
        result = "".join(output)
        if len(result.encode("utf-8")) > _MAX_BYTES:
            raise ValueError("MOL2 reordered output exceeds 16 MiB")
        _verify(record, json.loads(parse_mol2_record_json(result)), order, mapping)
    else:
        result = text
    return Mol2HydrogenOrder(
        text,
        result,
        tuple(order),
        tuple(inverse),
        tuple(mapping.items()),
        tuple(parents),
    )
