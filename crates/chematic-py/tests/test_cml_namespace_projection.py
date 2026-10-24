"""CML namespace identity, graph scope and XML validity precede projection."""
import json
from pathlib import Path
import struct
import pytest
import chematic

PACKET=json.loads(
    (Path(__file__).parents[3]/"validation/cml_namespace_projection/cases.json")
    .read_text(encoding="utf-8")
)

def layout(text,variant):
    if variant=="compact":return text
    lf=text.replace("><",">\n<")
    return lf.replace("\n","\r\n") if variant=="crlf" else lf

def bits(rows):return [struct.pack(">d",float(v)) for row in rows for v in row]

@pytest.mark.parametrize("case",PACKET["valid"],ids=lambda c:c["id"])
@pytest.mark.parametrize("variant",PACKET["variants"])
def test_projection_selects_cml_graph_by_namespace_and_scope(case,variant):
    text=layout(case["input"],variant)
    assert text.encode("utf-8").decode("utf-8")==text
    mol,xyz=chematic.from_cml_with_3d(text.encode())
    expected=case["expected"]
    assert mol.atom_count==len(expected["numbers"])
    assert bits(xyz)==bits(expected["coordinates"])
    if case['id'] in {'foreign-atom','foreign-bond','foreign-attributes','foreign-ancestor',
                       'namespace-reset-note','rebound-prefix','empty-namespace-atom',
                       'formula-is-not-graph','comment-cdata'}:
        # The reader still selects the graph, while checked conversion refuses
        # annotations that its CJSON target cannot represent.
        with pytest.raises(ValueError):chematic.convert_format(text,'cml','cjson')
        return
    data=json.loads(chematic.convert_format(text,"cml","cjson"))
    a,b=data["atoms"],data["bonds"]
    assert bits([a["coords"]["3d"]])==bits(expected["coordinates"])
    assert a["elements"]["number"]==expected["numbers"]
    assert a.get("formalCharges",[0]*mol.atom_count)==expected["formal_charges"]
    assert a.get("isotopes",[0]*mol.atom_count)==expected["isotopes"]
    assert [[b["connections"]["index"][2*i],b["connections"]["index"][2*i+1],o] for i,o in enumerate(b["order"])]==expected["bonds"]

def test_unicode_namespace_prefix_has_identical_str_and_bytes_projection():
    case=next(c for c in PACKET["valid"] if c["id"]=="unicode-prefix")
    text=case["input"]
    from_text=chematic.from_cml_with_3d(text)
    from_bytes=chematic.from_cml_with_3d(text.encode("utf-8"))
    assert from_text[0].smiles==from_bytes[0].smiles
    assert bits(from_text[1])==bits(from_bytes[1])

@pytest.mark.parametrize("case",PACKET["invalid"],ids=lambda c:c["id"])
@pytest.mark.parametrize("variant",PACKET["variants"])
def test_ambiguous_or_invalid_input_never_returns_projection(case,variant):
    text=layout(case["input"],variant)
    with pytest.raises(ValueError):chematic.from_cml_with_3d(text)
    with pytest.raises(ValueError):chematic.convert_format(text,"cml","cjson")
