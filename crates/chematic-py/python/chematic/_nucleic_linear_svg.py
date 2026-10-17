"""Linear sequence diagrams for validated explicit nucleic-acid documents."""

from __future__ import annotations

from collections.abc import Sequence as RuntimeSequence
import json
from typing import TYPE_CHECKING, Any, Sequence, Union
from xml.etree.ElementTree import Element, SubElement, tostring

if TYPE_CHECKING:
    from . import nucleic_acid_validate_json
else:
    from .chematic import nucleic_acid_validate_json

MAX_RENDER_RESIDUES = 512
MAX_RENDER_STRANDS = 8
COLORS = {
    "A": "#25754a",
    "C": "#2466a6",
    "G": "#a86612",
    "U": "#ae3d48",
    "T": "#814a9d",
    "other": "#555555",
}


def _unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    data: dict[str, Any] = {}
    for key, value in pairs:
        if key in data:
            raise ValueError(f"duplicate JSON key: {key!r}")
        data[key] = value
    return data


def _source(value: Union[str, bytes]) -> str:
    if type(value) not in (str, bytes):
        raise TypeError("document_json must be UTF-8 str or bytes")
    if len(value) > 1_048_576:
        raise ValueError("document exceeds 1 MiB")
    try:
        text = value.decode("utf-8") if isinstance(value, bytes) else value
        if len(text.encode("utf-8")) > 1_048_576:
            raise ValueError("document exceeds 1 MiB")
    except UnicodeError as exc:
        raise ValueError("document must be valid UTF-8") from exc
    try:
        json.loads(text, object_pairs_hook=_unique)
    except (json.JSONDecodeError, RecursionError) as exc:
        raise ValueError("invalid or excessively nested document JSON") from exc
    return text


def _xml_text(value: str) -> str:
    if any(
        not (
            c in "\t\r\n"
            or 0x20 <= ord(c) <= 0xD7FF
            or 0xE000 <= ord(c) <= 0xFFFD
            or 0x10000 <= ord(c) <= 0x10FFFF
        )
        for c in value
    ):
        raise ValueError("SVG text must contain valid XML 1.0 characters")
    return value


def _pair_layout(
    strands: list[dict[str, Any]],
    base_pairs: Sequence[tuple[str, str]],
    antiparallel: bool,
) -> tuple[int, list[tuple[str, str, int]]]:
    if not isinstance(base_pairs, RuntimeSequence) or isinstance(
        base_pairs, (str, bytes)
    ):
        raise TypeError("base_pairs must be a sequence of residue-ID pairs")
    if len(base_pairs) > MAX_RENDER_RESIDUES // 2:
        raise ValueError("too many base pairs")
    if base_pairs and not antiparallel:
        raise ValueError("base_pairs require antiparallel=True and two strands")
    index = {
        r["id"]: (s, i, r)
        for s, strand in enumerate(strands)
        for i, r in enumerate(strand["residues"])
    }
    used: set[str] = set()
    offset = None
    result: list[tuple[str, str, int]] = []
    for pair in base_pairs:
        if (
            not isinstance(pair, (list, tuple))
            or len(pair) != 2
            or any(type(v) is not str for v in pair)
        ):
            raise TypeError("each base pair requires two residue-ID strings")
        a, b = pair
        if a not in index or b not in index:
            raise ValueError("base pair refers to an unknown residue")
        sa, ia, ra = index[a]
        sb, ib, rb = index[b]
        if sa == sb:
            raise ValueError("base pairs must connect the two distinct strands")
        if a in used or b in used:
            raise ValueError("a residue cannot occur in more than one base pair")
        if ra.get("modification") or rb.get("modification"):
            raise ValueError(
                "canonical pair marks are not supported for modified residues"
            )
        bases = frozenset((ra["base"], rb["base"]))
        if bases in (frozenset(("A", "U")), frozenset(("A", "T"))):
            count = 2
        elif bases == frozenset(("C", "G")):
            count = 3
        else:
            raise ValueError("only canonical A-U, A-T and C-G pair marks are supported")
        if sa == 1:
            a, b, ia, ib = b, a, ib, ia
        required = ia - (len(strands[1]["residues"]) - 1 - ib)
        if offset is not None and required != offset:
            raise ValueError("pairs require a single ungapped antiparallel register")
        offset = required
        used.update((a, b))
        result.append((a, b, count))
    return offset or 0, result


def _line(parent: Element, x1: int, y1: int, x2: int, y2: int, **attrs: str) -> Element:
    return SubElement(
        parent,
        "line",
        {"x1": str(x1), "y1": str(y1), "x2": str(x2), "y2": str(y2), **attrs},
    )


def _label(parent: Element, x: int, y: int, text: str, **attrs: str) -> None:
    node = SubElement(
        parent,
        "text",
        {
            "x": str(x),
            "y": str(y),
            "text-anchor": "middle",
            "dominant-baseline": "middle",
            **attrs,
        },
    )
    node.text = text


def nucleic_acid_linear_svg(
    document_json: Union[str, bytes],
    *,
    antiparallel: bool = False,
    base_pairs: Sequence[tuple[str, str]] = (),
) -> str:
    """Draw a linear sequence diagram from explicit nucleic-acid metadata.

    Args:
        document_json: Existing chematic.nucleic-acid.v1 JSON, validated by
            the native document API. Full original text is embedded in SVG.
        antiparallel: For exactly two strands, draw the second right-to-left
            and make bases face inward. Otherwise every strand runs left-to-right.
        base_pairs: Explicit residue-ID pairs. Requires antiparallel mode.
            Unmodified canonical A-U/A-T/C-G pairs in one ungapped register
            receive two/three dashed symbolic H-bond marks. Unpaired overhangs
            are allowed. No complementarity or secondary structure is inferred.

    Returns:
        A self-contained SVG with fixed 64-unit residue spacing, upright base
        glyphs on one side of each straight backbone, 5′/3′ labels, stable source
        IDs, explicit phosphodiester links and the original JSON in metadata.
        Modifications use an asterisk and a full tooltip. This is a sequence
        diagram, not full atom/bond depiction or a change to Mol.svg().

    Raises:
        TypeError: Unsupported input/option/pair types.
        ValueError: Invalid document, unsupported pairs/register, XML control
            characters, or size limits (1 MiB JSON, 8 strands, 512 residues).
    """
    if type(antiparallel) is not bool:
        raise TypeError("antiparallel must be bool")
    source = _source(document_json)
    envelope = json.loads(nucleic_acid_validate_json(source))
    if not envelope["ok"]:
        error = envelope["error"]
        raise ValueError(f"{error['code']} at {error['path']}")
    document = envelope["document"]
    strands = document["strands"]
    if (
        len(strands) > MAX_RENDER_STRANDS
        or sum(len(s["residues"]) for s in strands) > MAX_RENDER_RESIDUES
    ):
        raise ValueError("linear diagram supports at most 8 strands and 512 residues")
    if antiparallel and len(strands) != 2:
        raise ValueError("antiparallel=True requires exactly two strands")
    offset, pairs = _pair_layout(strands, base_pairs, antiparallel)
    # XML-invalid IDs fail instead of corrupting an otherwise valid document.
    for strand in strands:
        _xml_text(strand["id"])
        for residue in strand["residues"]:
            _xml_text(residue["id"])
            _xml_text(residue.get("modification") or "")
            for atom_ref in residue["atom_refs"]:
                _xml_text(atom_ref)
    for link in document["linkages"]:
        _xml_text(link["id"])
    columns: dict[str, int] = {}
    for s, strand in enumerate(strands):
        n = len(strand["residues"])
        for i, residue in enumerate(strand["residues"]):
            columns[residue["id"]] = (
                n - 1 - i + offset if antiparallel and s == 1 else i
            )
    lo, hi = min(columns.values()), max(columns.values())
    width = 192 + 64 * (hi - lo)
    height = 320 if antiparallel else 170 * len(strands)
    root = Element(
        "svg",
        {
            "xmlns": "http://www.w3.org/2000/svg",
            "viewBox": f"0 0 {width} {height}",
            "width": str(width),
            "height": str(height),
            "role": "img",
            "aria-labelledby": "nucleic-title nucleic-description",
        },
    )
    SubElement(root, "title", {"id": "nucleic-title"}).text = (
        "Linear nucleic-acid sequence diagram"
    )
    SubElement(root, "desc", {"id": "nucleic-description"}).text = (
        "Explicit residue order and phosphodiester links; letters represent bases. Pair marks, when present, describe caller-declared canonical pairs, not measured hydrogen bonds."
    )
    SubElement(
        root,
        "metadata",
        {"id": "nucleic-source", "data-format": "chematic.nucleic-acid.v1"},
    ).text = _xml_text(source)
    SubElement(root, "metadata", {"id": "nucleic-pairs"}).text = json.dumps(
        {"antiparallel": antiparallel, "base_pairs": [[a, b] for a, b, _ in pairs]},
        ensure_ascii=True,
    )
    SubElement(
        root, "rect", {"width": str(width), "height": str(height), "fill": "white"}
    )
    _label(
        root,
        width // 2,
        20,
        "Nucleic-acid sequence diagram",
        fill="#333333",
        **{"font-size": "14", "font-family": "sans-serif"},
    )
    locations: dict[str, tuple[int, int, int]] = {}
    for s, strand in enumerate(strands):
        y = (80 if s == 0 else 260) if antiparallel else 110 + 170 * s
        base_y = y + 50 if antiparallel and s == 0 else y - 50
        group = SubElement(
            root,
            "g",
            {
                "data-strand-id": strand["id"],
                "data-kind": strand["kind"],
                "data-direction": "reverse" if antiparallel and s == 1 else "forward",
                "font-family": "sans-serif",
            },
        )
        SubElement(group, "title").text = strand["id"]
        for i, residue in enumerate(strand["residues"]):
            rid = residue["id"]
            x = 96 + 64 * (columns[rid] - lo)
            locations[rid] = (x, y, base_y)
            node = SubElement(
                group,
                "g",
                {
                    "data-residue-id": rid,
                    "data-position": str(i + 1),
                    "data-atom-refs": json.dumps(
                        residue["atom_refs"], ensure_ascii=True
                    ),
                },
            )
            modification = residue.get("modification")
            SubElement(node, "title").text = f"{rid}: {residue['base']}" + (
                f" ({modification})" if modification else ""
            )
            _line(
                node,
                x,
                y,
                x,
                base_y,
                stroke="#666666",
                **{"stroke-width": "2", "data-kind": "base-stem"},
            )
            SubElement(
                node,
                "circle",
                {"cx": str(x), "cy": str(y), "r": "4", "fill": "#555555"},
            )
            SubElement(
                node,
                "rect",
                {
                    "x": str(x - 21),
                    "y": str(base_y - 18),
                    "width": "42",
                    "height": "36",
                    "rx": "6",
                    "fill": COLORS[residue["base"]],
                },
            )
            glyph = ("?" if residue["base"] == "other" else residue["base"]) + (
                "*" if modification else ""
            )
            _label(
                node,
                x,
                base_y,
                glyph,
                fill="white",
                **{"font-size": "18", "font-weight": "bold"},
            )
        first, last = strand["residues"][0]["id"], strand["residues"][-1]["id"]
        direction = -1 if antiparallel and s == 1 else 1
        _label(
            group,
            locations[first][0] - 38 * direction,
            y,
            "5′",
            fill="#222222",
            **{"font-size": "14"},
        )
        _label(
            group,
            locations[last][0] + 38 * direction,
            y,
            "3′",
            fill="#222222",
            **{"font-size": "14"},
        )
        _label(
            group,
            width // 2,
            y + (29 if base_y < y else -25),
            f"{strand['kind'].upper()} {s + 1}",
            fill="#555555",
            **{"font-size": "12"},
        )
    for link in document["linkages"]:
        start_location, end_location = (
            locations[link["three_prime_residue_id"]],
            locations[link["five_prime_residue_id"]],
        )
        _line(
            root,
            start_location[0],
            start_location[1],
            end_location[0],
            end_location[1],
            stroke="#555555",
            **{"stroke-width": "2", "data-linkage-id": link["id"]},
        )
    for a, b, count in pairs:
        x, _, ya = locations[a]
        _, _, yb = locations[b]
        group = SubElement(
            root,
            "g",
            {
                "data-pair": json.dumps([a, b], ensure_ascii=True),
                "data-symbolic-hbonds": str(count),
            },
        )
        for k in range(count):
            dx = (2 * k - count + 1) * 3
            _line(
                group,
                x + dx,
                ya + 22,
                x + dx,
                yb - 22,
                stroke="#777777",
                **{"stroke-dasharray": "3 3", "stroke-width": "1.25"},
            )
    # XML line-ending normalization would otherwise change source CR/CRLF.
    return tostring(root, encoding="unicode").replace("\r", "&#13;")
