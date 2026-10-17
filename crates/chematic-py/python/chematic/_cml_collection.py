"""Bounded source-subtree CML collection framing, without graph projection."""

from __future__ import annotations

from collections.abc import Sequence as RuntimeSequence
from dataclasses import dataclass, field
import re
from typing import Any, Optional, Sequence, Union
from xml.parsers import expat
from xml.sax.saxutils import quoteattr

CML_NS = "http://www.xml-cml.org/schema"
MAX_INPUT_BYTES = 16 * 1024 * 1024
MAX_RECORD_BYTES = 4 * 1024 * 1024
MAX_MOLECULES = 1000
MAX_ELEMENTS = 200_000
_ID = re.compile(r"[A-Za-z][A-Za-z0-9_.-]*\Z")


def _text(value: Union[str, bytes], limit: int) -> str:
    if type(value) not in (str, bytes):
        raise TypeError("CML input must be UTF-8 str or bytes")
    if len(value) > limit:
        raise ValueError("CML input exceeds byte limit")
    try:
        text = value.decode("utf-8") if isinstance(value, bytes) else value
        if len(text.encode("utf-8")) > limit:
            raise ValueError("CML input exceeds byte limit")
    except UnicodeError as exc:
        raise ValueError("CML input must be valid UTF-8") from exc
    return text


def _end_tag(raw: bytes, start: int) -> int:
    quote = 0
    for i in range(start, len(raw)):
        c = raw[i]
        if quote:
            if c == quote:
                quote = 0
        elif c in (34, 39):
            quote = c
        elif c == 62:
            return i + 1
    raise ValueError("unterminated XML tag")


def _expanded(name: str) -> tuple[str, str]:
    namespace, sep, local = name.rpartition("}")
    return (namespace, local) if sep else ("", name)


def _insert(raw: bytes, attributes: str) -> bytes:
    # Start-element offsets come from Expat, never an XML-looking regex match.
    end = _end_tag(raw, 0)
    at = end - 2 if raw[end - 2 : end] == b"/>" else end - 1
    return raw[:at] + attributes.encode("utf-8") + raw[at:]


@dataclass(frozen=True)
class CmlMoleculeRecord:
    """One source-ordered molecular subtree, with its complete input text.

    molecule_xml is a standalone molecule document: inherited namespace
    declarations are materialized on its opening tag. No atom, bond, stereo,
    coordinate, property or other payload inside the molecule is reserialized.
    source_text retains the original complete document, including XML preamble
    and collection-level comments/PIs. Records do not interpret molecular graphs.
    """

    source_text: str = field(repr=False)
    index: int
    molecule_id: Optional[str]
    molecule_xml: str = field(repr=False)


def _parse(text: str) -> tuple[bool, list[CmlMoleculeRecord]]:
    raw = text.encode("utf-8")
    parser = expat.ParserCreate(namespace_separator="}")
    parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_NEVER)
    depth = elements = 0
    collection = False
    root_namespaces: dict[Optional[str], str] = {}
    pending: dict[Optional[str], str] = {}
    spans: list[tuple[int, int, dict[Optional[str], str], Optional[str]]] = []
    active: Optional[tuple[int, int, bool, dict[Optional[str], str], Optional[str]]] = (
        None
    )
    ids: set[str] = set()

    def forbidden(*args: Any) -> Any:
        raise ValueError("CML DTD and entity declarations are not supported")

    def declaration(version: str, encoding: Optional[str], standalone: int) -> None:
        if version != "1.0" or (
            encoding is not None
            and encoding.lower() not in ("utf-8", "utf8", "us-ascii", "ascii")
        ):
            raise ValueError("CML requires XML 1.0 and UTF-8/ASCII encoding")

    def namespace(prefix: Optional[str], uri: Optional[str]) -> None:
        pending[prefix] = uri or ""

    def start(name: str, attrs: dict[str, str]) -> None:
        nonlocal depth, elements, collection, active
        depth += 1
        elements += 1
        if depth > 64 or elements > MAX_ELEMENTS:
            raise ValueError("CML XML depth/element limit exceeded")
        ns, local = _expanded(name)
        if depth == 1:
            if ns != CML_NS or local not in ("cml", "molecule"):
                raise ValueError("CML root must be a namespaced cml or molecule")
            collection = local == "cml"
            root_namespaces.update(pending)
            if collection and any(k != "version" or v != "3" for k, v in attrs.items()):
                raise ValueError("collection root supports only optional version='3'")
        molecule_depth = 2 if collection else 1
        if collection and depth == 2 and (ns, local) != (CML_NS, "molecule"):
            raise ValueError("CML collection children must be molecule elements")
        if depth > molecule_depth and (ns, local) == (CML_NS, "molecule"):
            raise ValueError("nested molecule hierarchies are outside this profile")
        if depth == molecule_depth:
            if len(spans) >= MAX_MOLECULES:
                raise ValueError("CML collection exceeds 1000 molecules")
            mol_id = attrs.get("id")
            if mol_id is not None:
                if not _ID.fullmatch(mol_id):
                    raise ValueError("invalid CML molecule id")
                if mol_id in ids:
                    raise ValueError("duplicate CML molecule id")
                ids.add(mol_id)
            begin = parser.CurrentByteIndex
            end = _end_tag(raw, begin)
            inherited = (
                {k: v for k, v in root_namespaces.items() if k not in pending}
                if collection
                else {}
            )
            active = (begin, end, raw[end - 2 : end] == b"/>", inherited, mol_id)
        pending.clear()

    def end(name: str) -> None:
        nonlocal depth, active
        molecule_depth = 2 if collection else 1
        if depth == molecule_depth:
            assert active is not None
            begin, open_end, self_closed, inherited, mol_id = active
            finish = open_end if self_closed else _end_tag(raw, parser.CurrentByteIndex)
            if finish - begin > MAX_RECORD_BYTES:
                raise ValueError("CML molecule exceeds 4 MiB")
            spans.append((begin, finish, inherited, mol_id))
            active = None
        depth -= 1

    def characters(data: str) -> None:
        if collection and depth <= 1 and data.strip():
            raise ValueError("CML collection cannot contain non-whitespace text")

    parser.StartDoctypeDeclHandler = forbidden
    parser.EntityDeclHandler = forbidden
    parser.ExternalEntityRefHandler = forbidden
    parser.XmlDeclHandler = declaration
    parser.StartNamespaceDeclHandler = namespace
    parser.StartElementHandler = start
    parser.EndElementHandler = end
    parser.CharacterDataHandler = characters
    try:
        parser.Parse(raw, True)
    except expat.ExpatError as exc:
        raise ValueError(f"invalid CML XML: {exc}") from exc
    finally:
        # Both callbacks close over parser for its byte position. Break those
        # references on success and failure instead of retaining input buffers
        # until cyclic garbage collection during a large writer batch.
        parser.StartElementHandler = None
        parser.EndElementHandler = None
    if not spans:
        raise ValueError("CML requires at least one molecule")
    records = []
    materialized_size = 0
    for i, (begin, finish, inherited, mol_id) in enumerate(spans):
        fragment = raw[begin:finish]
        additions = "".join(
            " xmlns" + (":" + prefix if prefix else "") + "=" + quoteattr(uri)
            for prefix, uri in inherited.items()
        )
        standalone_size = len(fragment) + len(additions.encode("utf-8"))
        materialized_size += standalone_size
        if standalone_size > MAX_RECORD_BYTES or materialized_size > MAX_INPUT_BYTES:
            raise ValueError("materialized CML namespace/byte limit exceeded")
        if additions:
            fragment = _insert(fragment, additions)
        records.append(CmlMoleculeRecord(text, i, mol_id, fragment.decode("utf-8")))
    return collection, records


def read_cml_collection(text: Union[str, bytes]) -> list[CmlMoleculeRecord]:
    """Read every molecular subtree from one CML XML document in input order.

    Accepts a CML-namespaced molecule, or cml root containing only direct
    molecules plus whitespace/comments/PIs. Collection root attributes are
    limited to optional version='3'; nested molecules and global CML metadata
    are rejected rather than flattened. Each record retains the complete input
    and a standalone subtree with inherited namespaces materialized. Molecular
    payload is opaque: this is XML framing, not chemical/schema validation.

    Limits: UTF-8 input and total standalone fragments each 16 MiB,
    1..1000 molecules, source/standalone subtree each 4 MiB, XML depth 64
    and 200000 elements. DTD/entities and non-UTF-8 declarations
    raise ValueError. Wrong input types raise TypeError. No partial list is
    returned after a malformed later molecule. Supplied molecule IDs must be
    unique; atom IDs may repeat across sibling molecule scopes.
    """
    return _parse(_text(text, MAX_INPUT_BYTES))[1]


def write_cml_collection(molecule_blocks: Sequence[Union[str, bytes]]) -> str:
    """Frame standalone CML molecule documents as one complete XML collection.

    Supply serialized molecule blocks, for example [m.to_cml() for m in mols].
    Each block must have one CML-namespaced molecule root, not a collection.
    Molecular subtrees are retained verbatim except for generated missing
    molecule IDs and an empty default-namespace shield for prefixed roots.
    Existing IDs and atom/bond references are never renamed. Native to_cml's
    chemistry/projection limitations are unchanged; this API adds framing.

    Returns XML 1.0 UTF-8 text with one cml root and all molecules in order.
    Missing molecule IDs receive unused m1/m2/... IDs; duplicate supplied IDs
    raise. Atom IDs are molecule-scoped and may repeat across siblings.
    Per-block preamble/comments/PIs outside the molecule root are omitted.
    No path/stream is written and no output is returned on a later failure.
    Limits follow read_cml_collection; input and final output are each <=16 MiB.
    """
    if not isinstance(molecule_blocks, RuntimeSequence) or isinstance(
        molecule_blocks, (str, bytes)
    ):
        raise TypeError("molecule_blocks must be a sequence of str/bytes CML documents")
    if not 1 <= len(molecule_blocks) <= MAX_MOLECULES:
        raise ValueError("CML collection requires 1..1000 molecule blocks")
    blocks: list[str] = []
    size = 0
    for block in molecule_blocks:
        text = _text(block, MAX_RECORD_BYTES)
        size += len(text.encode("utf-8"))
        if size > MAX_INPUT_BYTES:
            raise ValueError("total CML inputs exceed 16 MiB")
        blocks.append(text)
    records: list[CmlMoleculeRecord] = []
    supplied: set[str] = set()
    for text in blocks:
        is_collection, parsed = _parse(text)
        if is_collection:
            raise ValueError("each block requires a molecule root, not cml collection")
        record = parsed[0]
        if record.molecule_id is not None:
            if record.molecule_id in supplied:
                raise ValueError("duplicate supplied molecule id")
            supplied.add(record.molecule_id)
        records.append(record)
    fragments = []
    for i, record in enumerate(records):
        fragment = record.molecule_xml.encode("utf-8")
        attrs = ""
        if record.molecule_id is None:
            mol_id = f"m{i+1}"
            while mol_id in supplied:
                mol_id += "_"
            supplied.add(mol_id)
            attrs += ' id="' + mol_id + '"'
        # The collection introduces a default CML namespace. A prefixed
        # molecule with no default declaration must keep unqualified payload
        # elements in no namespace. Expat records declarations on that root.
        root_parser = expat.ParserCreate()
        root_attrs: Optional[dict[str, str]] = None

        def capture_root(name: str, attributes: dict[str, str]) -> None:
            nonlocal root_attrs
            if root_attrs is None:
                root_attrs = attributes

        root_parser.StartElementHandler = capture_root
        # Parse the complete valid fragment so namespace declarations can be
        # inspected without regex interpretation of quoted attribute values.
        root_parser.Parse(fragment, True)
        assert root_attrs is not None
        if "xmlns" not in root_attrs:
            attrs += ' xmlns=""'
        if attrs:
            fragment = _insert(fragment, attrs)
        fragments.append(fragment.decode("utf-8"))
    output = (
        '<?xml version="1.0" encoding="UTF-8"?>\n<cml xmlns="'
        + CML_NS
        + '">\n'
        + "\n".join(fragments)
        + "\n</cml>\n"
    )
    if len(output.encode("utf-8")) > MAX_INPUT_BYTES:
        raise ValueError("CML output exceeds 16 MiB")
    # Combined output must also satisfy the total element budget.
    _parse(output)
    return output
