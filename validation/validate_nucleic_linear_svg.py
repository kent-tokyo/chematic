"""Exhaustive short-sequence geometry and source-preservation checks."""

import hashlib
import itertools
import json
from pathlib import Path
from xml.etree import ElementTree as ET

from chematic import nucleic_acid_linear_svg

NS = {"s": "http://www.w3.org/2000/svg"}


def document(sequences, kind):
    # Independent authoring of explicit topology; no molecule/SMILES inference.
    result = {
        "schema": "chematic.nucleic-acid.v1",
        "atom_ids": [],
        "strands": [],
        "linkages": [],
        "annotations": {"synthetic": True},
    }
    for s, sequence in enumerate(sequences):
        strand = {"id": f"strand{s}", "kind": kind, "residues": [], "annotations": {}}
        for i, base in enumerate(sequence):
            rid = f"{s}:{i}"
            refs = [rid + ":5", rid + ":3"]
            result["atom_ids"].extend(refs)
            residue = {
                "id": rid,
                "base": base,
                "sugar": "ribose" if kind == "rna" else "deoxyribose",
                "atom_refs": refs,
                "annotations": {},
            }
            if i:
                residue["five_prime_linkage"] = f"{s}:link:{i-1}"
            if i + 1 < len(sequence):
                lid = f"{s}:link:{i}"
                residue["three_prime_linkage"] = lid
                result["linkages"].append(
                    {
                        "id": lid,
                        "kind": "phosphodiester",
                        "three_prime_residue_id": rid,
                        "five_prime_residue_id": f"{s}:{i+1}",
                        "three_prime_atom_ref": refs[1],
                        "five_prime_atom_ref": f"{s}:{i+1}:5",
                        "annotations": {},
                    }
                )
            strand["residues"].append(residue)
        result["strands"].append(strand)
    return result


def check(source, sequences, paired):
    options = (
        {
            "antiparallel": True,
            "base_pairs": [
                (f"0:{i}", f"1:{len(sequences[0])-1-i}")
                for i in range(len(sequences[0]))
            ],
        }
        if paired
        else {}
    )
    svg = nucleic_acid_linear_svg(source, **options)
    root = ET.fromstring(svg)
    assert root.find("s:metadata[@id='nucleic-source']", NS).text == source
    width, height = float(root.attrib["width"]), float(root.attrib["height"])
    rectangles = []
    for s, sequence in enumerate(sequences):
        group = root.findall("s:g[@data-strand-id]", NS)[s]
        nodes = group.findall("s:g[@data-residue-id]", NS)
        assert [node.find("s:text", NS).text for node in nodes] == list(sequence)
        x_values = []
        backbone_y = set()
        side = set()
        for node in nodes:
            line = node.find("s:line", NS)
            x, y1, y2 = map(
                float, [line.attrib["x1"], line.attrib["y1"], line.attrib["y2"]]
            )
            x_values.append(x)
            backbone_y.add(y1)
            side.add(y2 > y1)
            rect = node.find("s:rect", NS)
            left, top, w, h = [
                float(rect.attrib[k]) for k in ["x", "y", "width", "height"]
            ]
            assert 0 <= left < left + w <= width and 0 <= top < top + h <= height
            rectangles.append((left, top, left + w, top + h))
        direction = -1 if paired and s == 1 else 1
        assert all((b - a) * direction == 64 for a, b in zip(x_values, x_values[1:]))
        assert len(backbone_y) == len(side) == 1
        assert side == ({True} if paired and s == 0 else {False})
    for i, a in enumerate(rectangles):
        for b in rectangles[i + 1 :]:
            assert a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1]
    assert len(root.findall("s:line[@data-linkage-id]", NS)) == sum(
        len(s) - 1 for s in sequences
    )
    groups = root.findall("s:g[@data-symbolic-hbonds]", NS)
    assert len(groups) == (len(sequences[0]) if paired else 0)
    for i, group in enumerate(groups):
        expected = 3 if sequences[0][i] in "CG" else 2
        assert int(group.attrib["data-symbolic-hbonds"]) == expected
        assert len(group.findall("s:line", NS)) == expected
        for line in group.findall("s:line", NS):
            assert line.attrib["x1"] == line.attrib["x2"]
            assert float(line.attrib["y1"]) < float(line.attrib["y2"])
    return len(rectangles), len(groups)


def main():
    cases = residues = pairs = 0
    checksum = hashlib.sha256()
    for kind, alphabet, complement in [
        ("rna", "ACGU", str.maketrans("ACGU", "UGCA")),
        ("dna", "ACGT", str.maketrans("ACGT", "TGCA")),
    ]:
        for n in range(1, 5):
            for letters in itertools.product(alphabet, repeat=n):
                sequence = "".join(letters)
                for paired in [False, True]:
                    sequences = (
                        [sequence, sequence.translate(complement)[::-1]]
                        if paired
                        else [sequence]
                    )
                    source = json.dumps(document(sequences, kind), indent=2).replace(
                        "\n", "\r\n" if n % 2 == 0 else "\n"
                    )
                    checksum.update(source.encode())
                    r, p = check(source, sequences, paired)
                    cases += 1
                    residues += r
                    pairs += p
    print(
        json.dumps(
            {
                "cases": cases,
                "residues_checked": residues,
                "explicit_pairs_checked": pairs,
                "input_sha256": checksum.hexdigest(),
                "source_roundtrip": "exact decoded JSON text, including CRLF",
                "geometry": "straight backbones, 64-unit spacing, uniform base side, opposite duplex direction, non-overlapping in-bounds glyph boxes",
                "scope": "all length-1..4 canonical RNA/DNA sequences and their explicit reverse complements; sequence diagrams, not atomic chemical depictions or measured secondary structure",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
