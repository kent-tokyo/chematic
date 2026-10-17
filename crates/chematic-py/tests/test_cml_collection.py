"""CML collection framing and subtree preservation for Open Babel #2339."""

from dataclasses import FrozenInstanceError
from xml.etree import ElementTree as ET

import pytest
import chematic
from chematic import read_cml_collection, write_cml_collection

NS = "http://www.xml-cml.org/schema"


def molecule(mid=None, prefix=""):
    tag = prefix + ":" if prefix else ""
    declaration = f"xmlns:{prefix}" if prefix else "xmlns"
    identity = f' id="{mid}"' if mid is not None else ""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<{tag}molecule {declaration}="{NS}"{identity} formalCharge="1">
  <!-- exact internal comment -->
  <{tag}atomArray>
    <{tag}atom id="a1" elementType="N" formalCharge="1" hydrogenCount="4" x3="1.1234567890123457" y3="-2.5" z3="3"/>
    <{tag}atom id="a2" elementType="H" isotopeNumber="2" hydrogenCount="0"/>
  </{tag}atomArray>
  <{tag}bondArray><{tag}bond id="b1" atomRefs2="a1 a2" order="1"><{tag}bondStereo atomRefs2="a1 a2">W</{tag}bondStereo></{tag}bond></{tag}bondArray>
  <{tag}propertyList><{tag}property dictRef="custom:energy" xmlns:custom="urn:custom"><{tag}scalar units="custom:unit">-1.25</{tag}scalar></{tag}property></{tag}propertyList>
</{tag}molecule>"""


def payload(text):
    node = ET.fromstring(text)

    def tree(element):
        return (
            element.tag,
            tuple(
                sorted(
                    (k, v)
                    for k, v in element.attrib.items()
                    if k != "id" or element is not node
                )
            ),
            element.text,
            tuple((tree(c), c.tail) for c in element),
        )

    return tree(node)


@pytest.mark.parametrize("kind", ["list", "tuple", "bytes", "mixed"])
def test_collection_has_one_root_all_records_and_local_ids(kind):
    blocks = [molecule(), molecule("explicit")]
    if kind == "tuple":
        blocks = tuple(blocks)
    if kind == "bytes":
        blocks = [b.encode() for b in blocks]
    if kind == "mixed":
        blocks[1] = blocks[1].encode()
    result = write_cml_collection(blocks)
    root = ET.fromstring(result)
    assert root.tag == f"{{{NS}}}cml" and len(root) == 2
    records = read_cml_collection(result.encode())
    assert [r.molecule_id for r in records] == ["m1", "explicit"]
    assert [r.index for r in records] == [0, 1]
    assert all(r.source_text == result for r in records)
    for i, record in enumerate(records):
        original = blocks[i].decode() if isinstance(blocks[i], bytes) else blocks[i]
        assert payload(record.molecule_xml) == payload(original)
        assert "<!-- exact internal comment -->" in record.molecule_xml
        assert (
            'id="a1"' in record.molecule_xml
            and 'atomRefs2="a1 a2"' in record.molecule_xml
        )
        assert "1.1234567890123457" in record.molecule_xml
    with pytest.raises(FrozenInstanceError):
        records[0].index = 9


def test_missing_ids_avoid_supplied_names_without_renaming_atoms():
    result = read_cml_collection(
        write_cml_collection([molecule(), molecule("m1"), molecule("m1_")])
    )
    assert [r.molecule_id for r in result] == ["m1__", "m1", "m1_"]
    assert [r.molecule_xml.count('id="a1"') for r in result] == [1, 1, 1]


@pytest.mark.parametrize("prefix", ["c", "mol", "cml"])
def test_prefixed_molecules_do_not_change_unqualified_payload_namespace(prefix):
    block = molecule("source", prefix).replace(
        f"</{prefix}:molecule>",
        '<opaque xmlns=""><item>value</item></opaque>' + f"</{prefix}:molecule>",
    )
    record = read_cml_collection(write_cml_collection([block]))[0]
    assert payload(record.molecule_xml) == payload(block)
    assert ET.fromstring(record.molecule_xml).find("opaque/item").text == "value"
    assert 'xmlns=""' in record.molecule_xml


def test_inherited_namespaces_and_child_overrides_are_materialized():
    source = f"""<c:cml xmlns:c="{NS}" xmlns:p="urn:parent">
<c:molecule id="one"><p:property p:name="p:value"/></c:molecule>
<c:molecule id="two" xmlns:p="urn:child"><p:property/></c:molecule>
</c:cml>"""
    records = read_cml_collection(source)
    one, two = [ET.fromstring(r.molecule_xml) for r in records]
    assert one.find("{urn:parent}property").attrib == {"{urn:parent}name": "p:value"}
    assert two.find("{urn:child}property") is not None
    assert "xmlns:c=" in records[0].molecule_xml
    assert records[0].source_text == records[1].source_text == source
    assembled = write_cml_collection([r.molecule_xml for r in records])
    assert [payload(r.molecule_xml) for r in read_cml_collection(assembled)] == [
        payload(r.molecule_xml) for r in records
    ]


@pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"])
def test_source_spelling_line_endings_unicode_comments_and_quoted_gt(ending):
    block = (
        molecule("unicode")
        .replace(' formalCharge="1"', ' title="α &amp; β &gt; γ" formalCharge="1"', 1)
        .replace("\n", ending)
    )
    record = read_cml_collection(block)[0]
    assert record.source_text == block
    fragment = block[block.index("<molecule") :]
    assert record.molecule_xml == fragment
    output = write_cml_collection([block])
    assert fragment in output


def test_self_closing_molecules_and_compact_adjacent_records():
    source = f'<cml xmlns="{NS}"><molecule id="one"/><molecule id="two"><atomArray/></molecule></cml>'
    records = read_cml_collection(source)
    assert len(records) == 2
    assert ET.fromstring(records[0].molecule_xml).attrib["id"] == "one"
    assert records[0].molecule_xml.endswith("/>")
    assert records[1].molecule_xml.endswith("</molecule>")


def test_single_molecule_bom_preamble_and_outside_comments():
    source = "\ufeff" + molecule() + "\n<!-- trailing source comment -->"
    record = read_cml_collection(source)[0]
    assert record.source_text == source and record.molecule_id is None
    assert record.molecule_xml.startswith("<molecule")
    assert "trailing source comment" not in write_cml_collection([source])


@pytest.mark.parametrize(
    "source",
    [
        "",
        "<",
        "<molecule/>",
        "<wrong/>",
        f'<cml xmlns="{NS}"/>',
        molecule() + molecule(),
        molecule().replace("</molecule>", ""),
        f'<cml xmlns="{NS}"><scalar>global metadata</scalar></cml>',
        f'<cml xmlns="{NS}" title="collection"><molecule/></cml>',
        f'<cml xmlns="{NS}">extra<molecule/></cml>',
        f'<cml xmlns="{NS}"><molecule id="same"/><molecule id="same"/></cml>',
        f'<molecule xmlns="{NS}" id="1invalid"/>',
        f'<molecule xmlns="{NS}"><molecule/></molecule>',
        '<!DOCTYPE molecule [<!ENTITY x "text">]>'
        + f'<molecule xmlns="{NS}">&x;</molecule>',
        '<!DOCTYPE molecule SYSTEM "https://example.invalid/not-read">'
        + f'<molecule xmlns="{NS}"/>',
        molecule().replace("UTF-8", "ISO-8859-1"),
        molecule().replace('version="1.0"', 'version="1.1"'),
    ],
)
def test_invalid_xml_profile_and_entities_raise(source):
    with pytest.raises(ValueError):
        read_cml_collection(source)


@pytest.mark.parametrize("value", [None, {}, [], 1, bytearray(b"xml")])
def test_reader_wrong_types(value):
    with pytest.raises(TypeError):
        read_cml_collection(value)


@pytest.mark.parametrize("value", [b"\xff", "\ud800"])
def test_invalid_utf8(value):
    with pytest.raises(ValueError):
        read_cml_collection(value)


@pytest.mark.parametrize("value", [None, {}, "xml", b"xml", (x for x in [])])
def test_writer_wrong_collection_types(value):
    with pytest.raises(TypeError):
        write_cml_collection(value)


def test_writer_rejects_bad_later_records_duplicate_ids_and_collection_input():
    for blocks in [
        [molecule(), "<bad>"],
        [molecule("same"), molecule("same")],
        [write_cml_collection([molecule()])],
    ]:
        with pytest.raises(ValueError):
            write_cml_collection(blocks)
    with pytest.raises(TypeError):
        write_cml_collection([molecule(), None])
    with pytest.raises(ValueError):
        write_cml_collection([])


def test_limits_without_native_graph_projection(monkeypatch):
    from chematic import _cml_collection as module

    with pytest.raises(ValueError, match="1000"):
        write_cml_collection([molecule()] * 1001)
    with pytest.raises(ValueError, match="byte limit"):
        read_cml_collection(" " * (16 * 1024 * 1024 + 1))
    with pytest.raises(ValueError, match="byte limit"):
        read_cml_collection("α" * (8 * 1024 * 1024 + 1))
    source = (
        f'<molecule xmlns="{NS}">'
        + "<property>" * 64
        + "</property>" * 64
        + "</molecule>"
    )
    with pytest.raises(ValueError, match="depth"):
        read_cml_collection(source)
    monkeypatch.setattr(module, "MAX_ELEMENTS", 7)
    with pytest.raises(ValueError, match="element"):
        read_cml_collection(molecule())


def test_native_blocks_get_complete_collection_instead_of_concatenated_xml():
    smiles = ["CCO", "[NH4+]", "[13CH4]", "[Na+].[Cl-]"]
    molecules = [chematic.from_smiles(s) for s in smiles]
    blocks = [m.to_cml() for m in molecules]
    with pytest.raises(ET.ParseError):
        ET.fromstring("".join(blocks))
    records = read_cml_collection(write_cml_collection(blocks))
    assert len(records) == 4
    assert [chematic.from_cml_strict(r.molecule_xml).formula for r in records] == [
        m.formula for m in molecules
    ]
    assert [
        chematic.from_cml_strict(r.molecule_xml).formal_charge for r in records
    ] == [m.formal_charge for m in molecules]


def test_materialized_namespaces_are_bounded_without_multiplicative_growth(monkeypatch):
    from chematic import _cml_collection as module

    source = (
        f'<cml xmlns="{NS}" xmlns:p="urn:'
        + "x" * 300
        + '">'
        + "".join(f'<molecule id="m{i}"/>' for i in range(10))
        + "</cml>"
    )
    assert len(source.encode()) < 1000
    monkeypatch.setattr(module, "MAX_INPUT_BYTES", 1000)
    with pytest.raises(ValueError, match="materialized"):
        read_cml_collection(source)
    monkeypatch.setattr(module, "MAX_RECORD_BYTES", 100)
    with pytest.raises(ValueError, match="materialized"):
        read_cml_collection(source)


def test_unicode_preamble_offsets_and_exact_source_subtree_boundaries():
    source = "<!-- α β 😀 -->\n" + molecule("source").split("?>", 1)[1]
    record = read_cml_collection(source)[0]
    assert record.source_text == source
    assert record.molecule_xml == source[source.index("<molecule") :]


def test_cdata_and_processing_instruction_do_not_create_false_molecule_boundaries():
    extra = "<scalar><![CDATA[</molecule><molecule>literal & value]]></scalar><?inside preserve?>"
    block = molecule().replace("</molecule>", extra + "</molecule>")
    output = write_cml_collection([block, molecule()])
    records = read_cml_collection(output)
    assert len(records) == 2
    assert extra in records[0].molecule_xml
    assert records[0].molecule_xml.endswith("</molecule>")
    assert payload(records[0].molecule_xml) == payload(block)
