"""Bounded, source-preserving pymatgen Structure JSON reader (no Monty imports)."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from typing import TYPE_CHECKING, Any, Union

if TYPE_CHECKING:
    from . import Lattice, PeriodicStructure, Site
else:
    from .chematic import Lattice, PeriodicStructure, Site

MAX_INPUT_BYTES = 16 * 1024 * 1024
MAX_RECORDS = 256
MAX_SITES = 50_000
MAX_SPECIES = 100_000


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _json_float(token: str) -> float:
    value = float(token)
    mantissa = token.lower().split("e")[0]
    if not math.isfinite(value) or (
        value == 0 and any(c in "123456789" for c in mantissa)
    ):
        raise ValueError("JSON real number overflows or underflows float64")
    return value


def _bad_constant(token: str) -> Any:
    raise ValueError(f"non-JSON numeric constant: {token}")


def _text(value: Union[str, bytes]) -> str:
    if type(value) not in (str, bytes):
        raise TypeError("input must be UTF-8 str or bytes")
    if len(value) > MAX_INPUT_BYTES:
        raise ValueError("input exceeds 16 MiB")
    try:
        text = value.decode("utf-8") if isinstance(value, bytes) else value
        if len(text.encode("utf-8")) > MAX_INPUT_BYTES:
            raise ValueError("input exceeds 16 MiB")
    except UnicodeError as exc:
        raise ValueError("input must be valid UTF-8") from exc
    depth = 0
    quoted = escaped = False
    for char in text:
        if quoted:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                quoted = False
        elif char == '"':
            quoted = True
        elif char in "[{":
            depth += 1
            if depth > 64:
                raise ValueError("JSON nesting exceeds 64 levels")
        elif char in "]}":
            depth -= 1
    return text


def _object(value: Any, path: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise ValueError(f"{path} must be a JSON object")
    return value


def _number(value: Any, path: str) -> float:
    if type(value) not in (int, float):
        raise ValueError(f"{path} must be a real number, excluding booleans")
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError(f"{path} exceeds float64 range") from exc
    if not math.isfinite(result) or abs(result) > 1e100:
        raise ValueError(f"{path} requires a finite magnitude <= 1e100")
    return result


def _triple(value: Any, path: str) -> tuple[float, float, float]:
    if type(value) is not list or len(value) != 3:
        raise ValueError(f"{path} requires exactly three numbers")
    return (_number(value[0], path), _number(value[1], path), _number(value[2], path))


def _extras(data: dict[str, Any], known: set[str], path: str) -> list[str]:
    return [f"{path}.{key}" for key in data if key not in known]


@dataclass(frozen=True)
class PymatgenStructureRecord:
    """One validated structure plus the original complete JSON response.

    Obtain records with parse_pymatgen_structures. source_text retains the
    complete input; record_index addresses its data array (zero for a bare
    Structure). structure_dict returns a fresh ordinary JSON dictionary,
    including charge, oxidation states, site properties and unknown fields.
    No serialized module/class is ever imported or instantiated.
    """

    source_text: str = field(repr=False)
    record_index: int
    projection_losses: tuple[str, ...]
    _structure_json: str = field(repr=False)
    _structure: PeriodicStructure = field(repr=False, compare=False)

    def structure_dict(self) -> dict[str, Any]:
        """Return a separate JSON-value dictionary for this structure."""
        result: dict[str, Any] = json.loads(self._structure_json)
        return result

    def to_periodic_structure(
        self, *, allow_metadata_loss: bool = False
    ) -> PeriodicStructure:
        """Project to native geometry, refusing unrepresented metadata by default.

        Lattice, fractional coordinates, site/species order, occupancies and
        labels are retained. Charge, oxidation/spin, properties and unknown
        fields remain in this record, not in the native projection. Inspect
        projection_losses and explicitly opt in when these fields are present.
        Response-level metadata always remains only in source_text.
        """
        if type(allow_metadata_loss) is not bool:
            raise TypeError("allow_metadata_loss must be bool")
        if self.projection_losses and not allow_metadata_loss:
            raise ValueError(
                "native projection loses metadata: "
                + ", ".join(self.projection_losses[:8])
            )
        return self._structure


def _structure(data: dict[str, Any]) -> tuple[PeriodicStructure, tuple[str, ...]]:
    if data.get("@class") not in (None, "Structure", "IStructure") or data.get(
        "@module"
    ) not in (None, "pymatgen.core.structure"):
        raise ValueError("only pymatgen Structure/IStructure documents are supported")
    lattice_data = _object(data.get("lattice"), "lattice")
    pbc = lattice_data.get("pbc", [True, True, True])
    if (
        type(pbc) is not list
        or len(pbc) != 3
        or any(type(v) is not bool for v in pbc)
        or not all(pbc)
    ):
        raise ValueError("only fully periodic pbc=[true,true,true] is supported")
    matrix_data = lattice_data.get("matrix")
    if type(matrix_data) is not list or len(matrix_data) != 3:
        raise ValueError("lattice.matrix requires three row vectors in angstrom")
    matrix = [list(_triple(row, "lattice.matrix")) for row in matrix_data]
    lattice = Lattice.from_matrix(matrix)
    derived = dict(zip(("a", "b", "c"), lattice.lengths))
    derived.update(zip(("alpha", "beta", "gamma"), lattice.angles_degrees))
    derived["volume"] = lattice.volume
    for key, expected in derived.items():
        if key in lattice_data and not math.isclose(
            _number(lattice_data[key], "lattice." + key),
            expected,
            rel_tol=1e-7,
            abs_tol=1e-7,
        ):
            raise ValueError(f"lattice.{key} disagrees with lattice.matrix")
    losses = _extras(
        data,
        {"@module", "@class", "lattice", "sites", "charge", "properties"},
        "structure",
    )
    losses += _extras(
        lattice_data,
        {"matrix", "pbc", "a", "b", "c", "alpha", "beta", "gamma", "volume"},
        "lattice",
    )
    if "charge" in data and data["charge"] is not None:
        _number(data["charge"], "charge")
        losses.append("structure.charge")
    properties = _object(data.get("properties", {}), "properties")
    if properties:
        losses.append("structure.properties")
    sites: list[Site] = []
    for i, value in enumerate(data["sites"]):
        path = f"sites[{i}]"
        site_data = _object(value, path)
        fractional = _triple(site_data.get("abc"), path + ".abc")
        cartesian = [
            sum(fractional[j] * matrix[j][k] for j in range(3)) for k in range(3)
        ]
        if not all(math.isfinite(v) for v in cartesian):
            raise ValueError(f"{path} Cartesian coordinates overflow")
        if "xyz" in site_data:
            declared = _triple(site_data["xyz"], path + ".xyz")
            if not all(
                math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-7)
                for a, b in zip(cartesian, declared)
            ):
                raise ValueError(f"{path}.xyz disagrees with abc @ lattice.matrix")
        label = site_data.get("label")
        if label is not None:
            if type(label) is not str:
                raise ValueError(f"{path}.label must be a string or null")
            try:
                label.encode("utf-8")
            except UnicodeError as exc:
                raise ValueError(f"{path}.label must be valid Unicode") from exc
        props = _object(site_data.get("properties", {}), path + ".properties")
        if props:
            losses.append(path + ".properties")
        losses += _extras(
            site_data, {"species", "abc", "xyz", "label", "properties"}, path
        )
        species: list[tuple[str, float]] = []
        for j, entry in enumerate(site_data["species"]):
            spec = _object(entry, path + ".species")
            element = spec.get("element")
            if type(element) is not str:
                raise ValueError(f"{path}.species element must be a symbol")
            try:
                element.encode("utf-8")
            except UnicodeError as exc:
                raise ValueError(
                    f"{path}.species element must be valid Unicode"
                ) from exc
            occupancy = _number(spec.get("occu"), path + ".species.occu")
            species.append((element, occupancy))
            losses += _extras(spec, {"element", "occu"}, f"{path}.species[{j}]")
        sites.append(Site(species, fractional, label))
    return PeriodicStructure(lattice, sites), tuple(losses)


def parse_pymatgen_structures(text: Union[str, bytes]) -> list[PymatgenStructureRecord]:
    """Read bare Structure JSON or a Materials Project data/structure response.

    UTF-8 text is parsed as ordinary JSON without pymatgen, Monty or network
    access. Each data entry must contain one supported structure; failures
    abort the entire call instead of skipping records. Matrix rows are lattice
    vectors in angstrom, abc fractional coordinates; xyz, when supplied, must
    agree within atol=1e-7 angstrom and rtol=1e-8. No wrapping or symmetry
    expansion is performed; all supplied sites are taken as the explicit cell.

    Limits: 16 MiB UTF-8, 64 JSON nesting levels, 256 structures, 50000 total
    sites, 32 species per site, 100000 total species. Duplicate keys, nonfinite
    numbers, float64 numeric overflow/underflow, malformed/empty structures,
    dummy species, partial PBC and non-Structure classes raise ValueError.
    Wrong input types raise TypeError. Native lattice/occupancy validation also
    applies. Optional redundant lengths/angles/volume must agree with the
    matrix within atol=rtol=1e-7. See projection_losses before native conversion.
    """
    source = _text(text)
    try:
        document = json.loads(
            source,
            object_pairs_hook=_pairs,
            parse_float=_json_float,
            parse_constant=_bad_constant,
        )
    except (json.JSONDecodeError, RecursionError) as exc:
        raise ValueError("invalid or excessively nested JSON input") from exc
    root = _object(document, "input")
    response_losses: list[str] = []
    entry_losses: list[list[str]] = []
    if "lattice" in root or "sites" in root:
        if "data" in root:
            raise ValueError("ambiguous Structure and response.data document")
        entries = [root]
        entry_losses = [[]]
    else:
        response_entries = root.get("data")
        if (
            type(response_entries) is not list
            or not 1 <= len(response_entries) <= MAX_RECORDS
        ):
            raise ValueError("response.data requires 1..256 entries")
        response_losses = ["response"] if any(key != "data" for key in root) else []
        entry_losses = [
            _extras(_object(entry, "data entry"), {"structure"}, f"data[{i}]")
            for i, entry in enumerate(response_entries)
        ]
        entries = [
            _object(_object(entry, "data entry").get("structure"), "data.structure")
            for entry in response_entries
        ]
    n_sites = n_species = 0
    # Bound all record/site/species collections before constructing native objects.
    for data in entries:
        sites = data.get("sites")
        if type(sites) is not list or not 1 <= len(sites) <= MAX_SITES:
            raise ValueError("structure.sites requires 1..50000 entries")
        n_sites += len(sites)
        for value in sites:
            species = _object(value, "site").get("species")
            if type(species) is not list or not 1 <= len(species) <= 32:
                raise ValueError("site.species requires 1..32 entries")
            n_species += len(species)
        if n_sites > MAX_SITES or n_species > MAX_SPECIES:
            raise ValueError("total site/species limit exceeded")
    records = []
    for i, data in enumerate(entries):
        geometry, losses = _structure(data)
        losses = losses + tuple(response_losses + entry_losses[i])
        records.append(
            PymatgenStructureRecord(
                source,
                i,
                losses,
                json.dumps(data, ensure_ascii=True, allow_nan=False),
                geometry,
            )
        )
    return records
