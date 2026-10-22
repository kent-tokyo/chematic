"""Linear document diagrams motivated by RDKit #9617 (not atom depiction)."""

import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest
from chematic import nucleic_acid_linear_svg

NS = {"s": "http://www.w3.org/2000/svg"}


def make_document(sequences=("ACGU",), kinds=None):
    kinds = kinds or ["rna"] * len(sequences)
    document = {
        "schema": "chematic.nucleic-acid.v1",
        "atom_ids": [],
        "strands": [],
        "linkages": [],
        "annotations": {"source": "original synthetic fixture"},
    }
    for s, (sequence, kind) in enumerate(zip(sequences, kinds)):
        strand = {"id": f"strand-{s}", "kind": kind, "residues": [], "annotations": {}}
        for i, base in enumerate(sequence):
            rid = f"s{s}r{i}"
            refs = [rid + "-5", rid + "-3"]
            document["atom_ids"].extend(refs)
            residue = {
                "id": rid,
                "base": base,
                "sugar": "ribose" if kind == "rna" else "deoxyribose",
                "atom_refs": refs,
                "annotations": {},
            }
            if i:
                residue["five_prime_linkage"] = f"s{s}l{i-1}"
            if i + 1 < len(sequence):
                lid = f"s{s}l{i}"
                residue["three_prime_linkage"] = lid
                document["linkages"].append(
                    {
                        "id": lid,
                        "kind": "phosphodiester",
                        "three_prime_residue_id": rid,
                        "five_prime_residue_id": f"s{s}r{i+1}",
                        "three_prime_atom_ref": refs[1],
                        "five_prime_atom_ref": f"s{s}r{i+1}-5",
                        "annotations": {},
                    }
                )
            strand["residues"].append(residue)
        document["strands"].append(strand)
    return document


def render(doc, **kwargs):
    return ET.fromstring(nucleic_acid_linear_svg(json.dumps(doc), **kwargs))


def residues(root):
    return {
        node.attrib["data-residue-id"]: node
        for node in root.findall(".//s:g[@data-residue-id]", NS)
    }


@pytest.mark.parametrize(
    "sequence,kind",
    [("A", "rna"), ("ACGU", "rna"), ("ACGTACGT", "dna"), ("A" * 64, "rna")],
)
def test_straight_backbone_uniform_base_side_and_source_metadata(sequence, kind):
    doc = make_document([sequence], [kind])
    source = json.dumps(doc, indent=2)
    text = nucleic_acid_linear_svg(source)
    assert text == nucleic_acid_linear_svg(source.encode())
    root = ET.fromstring(text)
    nodes = list(residues(root).values())
    stems = [node.find("s:line[@data-kind='base-stem']", NS) for node in nodes]
    assert [int(stem.attrib["x1"]) for stem in stems] == [
        96 + 64 * i for i in range(len(sequence))
    ]
    assert all(
        stem.attrib["y1"] == "110" and stem.attrib["y2"] == "60" for stem in stems
    )
    assert [node.find("s:text", NS).text for node in nodes] == list(sequence)
    assert len(root.findall("s:line[@data-linkage-id]", NS)) == len(sequence) - 1
    assert root.find("s:metadata[@id='nucleic-source']", NS).text == source
    assert root.findall(".//s:text", NS)[-3].text == "5′"
    assert root.findall(".//s:text", NS)[-2].text == "3′"


def test_antiparallel_pairs_are_explicit_in_register_and_symbolic():
    doc = make_document(["AGC", "GCU"])
    pairs = [("s0r0", "s1r2"), ("s0r1", "s1r1"), ("s0r2", "s1r0")]
    root = render(doc, antiparallel=True, base_pairs=pairs)
    nodes = residues(root)
    for a, b in pairs:
        sa = nodes[a].find("s:line", NS)
        sb = nodes[b].find("s:line", NS)
        assert sa.attrib["x1"] == sb.attrib["x1"]
        assert int(sa.attrib["y2"]) > int(sa.attrib["y1"])
        assert int(sb.attrib["y2"]) < int(sb.attrib["y1"])
    groups = root.findall("s:g[@data-symbolic-hbonds]", NS)
    assert [g.attrib["data-symbolic-hbonds"] for g in groups] == ["2", "3", "3"]
    assert [len(g.findall("s:line", NS)) for g in groups] == [2, 3, 3]
    assert len(root.findall("s:g[@data-direction='reverse']", NS)) == 1
    # Complementary-looking sequences alone never imply base pairs.
    assert render(doc, antiparallel=True).findall("s:g[@data-pair]", NS) == []


def test_dna_canonical_pairs_and_pair_argument_orientation():
    doc = make_document(["ATCG", "CGAT"], ["dna", "dna"])
    pairs = [(f"s1r{3-i}", f"s0r{i}") for i in range(4)]
    root = render(doc, antiparallel=True, base_pairs=pairs)
    assert [
        g.attrib["data-symbolic-hbonds"] for g in root.findall("s:g[@data-pair]", NS)
    ] == ["2", "2", "3", "3"]


def test_duplex_overhang_register_offset_and_multiple_independent_strands():
    doc = make_document(["AAGC", "GCU"])
    root = render(
        doc,
        antiparallel=True,
        base_pairs=[("s0r1", "s1r2"), ("s0r2", "s1r1"), ("s0r3", "s1r0")],
    )
    nodes = residues(root)
    assert (
        nodes["s0r1"].find("s:line", NS).attrib["x1"]
        == nodes["s1r2"].find("s:line", NS).attrib["x1"]
    )
    many = render(make_document(["AU", "CG", "GC"]))
    assert len(many.findall("s:g[@data-strand-id]", NS)) == 3
    assert all(
        g.attrib["data-direction"] == "forward"
        for g in many.findall("s:g[@data-strand-id]", NS)
    )


def test_modifications_and_full_document_metadata_are_not_flattened():
    root_path = Path(__file__).resolve().parents[3]
    fixture = json.loads(
        (root_path / "validation/nucleic_acid_document_contract.json").read_text(encoding="utf-8")
    )
    for case in fixture["valid_cases"]:
        source = json.dumps(case["document"], indent=2)
        root = ET.fromstring(nucleic_acid_linear_svg(source))
        assert root.find("s:metadata[@id='nucleic-source']", NS).text == source
        if case["id"] == "explicit_modification":
            modified = residues(root)["m1"]
            assert modified.find("s:text", NS).text == "C*"
            assert "5mC" in modified.find("s:title", NS).text
            assert json.loads(modified.attrib["data-atom-refs"]) == ["m1-5", "m1-3"]


def test_source_identifiers_and_annotations_are_xml_escaped():
    doc = make_document(["A"])
    residue = doc["strands"][0]["residues"][0]
    residue["id"] = 'r" onload="alert(1)<&'
    doc["strands"][0]["id"] = "<script>not code</script>"
    doc["annotations"] = {
        "xml": "</metadata><script>alert(1)</script>",
        "unicode": "α😀",
    }
    root = render(doc)
    assert root.findall(".//s:script", NS) == []
    assert list(residues(root)) == [residue["id"]]
    assert all("onload" not in n.attrib for n in root.iter())
    assert json.loads(root.find("s:metadata[@id='nucleic-source']", NS).text) == doc


@pytest.mark.parametrize(
    "pairs,exception",
    [
        ("ids", TypeError),
        (None, TypeError),
        ((x for x in []), TypeError),
        ([1], TypeError),
        ([("s0r0",)], TypeError),
        ([("s0r0", 1)], TypeError),
        ([("missing", "s1r0")], ValueError),
        ([("s0r0", "s0r1")], ValueError),
        ([("s0r0", "s1r2"), ("s0r0", "s1r2")], ValueError),
        ([("s0r0", "s1r0")], ValueError),
        ([("s0r0", "s1r2"), ("s0r2", "s1r1")], ValueError),
    ],
)
def test_invalid_pairs_are_rejected(pairs, exception):
    with pytest.raises(exception):
        render(make_document(["AGC", "GCU"]), antiparallel=True, base_pairs=pairs)


def test_pairs_require_duplex_mode_and_modifications_are_not_assumed_canonical():
    doc = make_document(["CG", "CG"])
    with pytest.raises(ValueError, match="antiparallel"):
        render(doc, base_pairs=[("s0r0", "s1r1")])
    doc["strands"][0]["residues"][0]["modification"] = "5mC"
    with pytest.raises(ValueError, match="modified"):
        render(doc, antiparallel=True, base_pairs=[("s0r0", "s1r1")])


@pytest.mark.parametrize("count", [1, 3])
def test_antiparallel_requires_exactly_two_strands(count):
    with pytest.raises(ValueError, match="exactly two"):
        render(make_document(["AU"] * count), antiparallel=True)


@pytest.mark.parametrize("value", [None, {}, [], 1, bytearray(b"{}")])
def test_wrong_input_types(value):
    with pytest.raises(TypeError):
        nucleic_acid_linear_svg(value)


@pytest.mark.parametrize("value", [None, 0, 1, "yes"])
def test_wrong_option_types(value):
    with pytest.raises(TypeError):
        render(make_document(), antiparallel=value)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "[]",
        "{}",
        "{",
        '{"a":1,"a":2}',
        '{"a":NaN}',
        b"\xff",
        "\ud800",
        "[" * 2000 + "0" + "]" * 2000,
    ],
)
def test_invalid_json_schema_and_encoding(text):
    with pytest.raises(ValueError):
        nucleic_acid_linear_svg(text)


def test_existing_document_invalid_cases_remain_invalid():
    root_path = Path(__file__).resolve().parents[3]
    fixture = json.loads(
        (root_path / "validation/nucleic_acid_document_contract.json").read_text(encoding="utf-8")
    )
    for case in fixture["invalid_cases"]:
        with pytest.raises(ValueError, match=case["expected_code"]):
            render(case["document"])


@pytest.mark.parametrize("control", ["\x01", "\x0b", "\uffff"])
def test_xml_invalid_identifiers_raise(control):
    doc = make_document(["A"])
    doc["strands"][0]["id"] += control
    with pytest.raises(ValueError, match="XML"):
        render(doc)


def test_resource_limits():
    for doc in [make_document(["A"] * 9), make_document(["A" * 513])]:
        with pytest.raises(ValueError, match="512"):
            render(doc)
    with pytest.raises(ValueError, match="1 MiB"):
        nucleic_acid_linear_svg(" " * (1_048_576 + 1))
    with pytest.raises(ValueError, match="1 MiB"):
        nucleic_acid_linear_svg("α" * (524_288 + 1))
    with pytest.raises(ValueError, match="too many"):
        render(
            make_document(["A", "U"]),
            antiparallel=True,
            base_pairs=[("s0r0", "s1r0")] * 257,
        )


@pytest.mark.parametrize("ending", ["\r\n", "\r", "\n"])
def test_original_json_line_endings_survive_xml_parsing(ending):
    source = json.dumps(make_document(), indent=2).replace("\n", ending)
    root = ET.fromstring(nucleic_acid_linear_svg(source))
    assert root.find("s:metadata[@id='nucleic-source']", NS).text == source


def test_negative_overhang_offset_keeps_registered_pairs_inside_viewbox():
    doc = make_document(["AGC", "GCUU"])
    root = render(
        doc,
        antiparallel=True,
        base_pairs=[("s0r0", "s1r2"), ("s0r1", "s1r1"), ("s0r2", "s1r0")],
    )
    nodes = residues(root)
    assert (
        nodes["s0r0"].find("s:line", NS).attrib["x1"]
        == nodes["s1r2"].find("s:line", NS).attrib["x1"]
    )
    assert nodes["s1r3"].find("s:line", NS).attrib["x1"] == "96"


def test_render_boundary_accepts_512_residues_and_eight_strands():
    for sequences in [["A" * 512], ["ACGU" * 16] * 8]:
        root = render(make_document(sequences))
        assert len(residues(root)) == 512
        assert len(root.findall("s:g[@data-strand-id]", NS)) == len(sequences)


def test_xml_invalid_literal_source_metadata_is_rejected_without_losing_it():
    doc = make_document(["A"])
    doc["annotations"]["noncharacter"] = "\uffff"
    with pytest.raises(ValueError, match="XML"):
        nucleic_acid_linear_svg(json.dumps(doc, ensure_ascii=False))
    source = json.dumps(doc, ensure_ascii=True)
    root = ET.fromstring(nucleic_acid_linear_svg(source))
    assert root.find("s:metadata[@id='nucleic-source']", NS).text == source
