import json
from pathlib import Path
import struct
import xml.etree.ElementTree as ET
import pytest
import chematic

PACKET=json.loads((Path(__file__).resolve().parents[3]/'validation/cjson_cml_projection/cases.json').read_text())
NS='{http://www.xml-cml.org/schema}'


@pytest.mark.parametrize('case',PACKET['valid'],ids=lambda c:c['id'])
def test_common_and_direct_preserve_graph_xyz_and_name(case):
    for source in [case['input'],case['input'].replace('\n','\r\n'),json.dumps(json.loads(case['input']),ensure_ascii=False,separators=(',',':'))]:
        output=chematic.project_cjson_to_cml(source)
        assert chematic.project_cjson_to_cml(source.encode())==output
        assert chematic.convert_format(source,'.CJSON','cml')==output
        root=ET.fromstring(output);atoms=root.findall('./'+NS+'atomArray/'+NS+'atom');expected=case['expected']
        assert [a.attrib['elementType'] for a in atoms]==[{1:'H',6:'C',7:'N',8:'O'}[n] for n in expected['numbers']]
        assert [int(a.get('formalCharge','0')) for a in atoms]==expected['charges']
        assert [int(a.get('isotopeNumber','0')) for a in atoms]==expected['isotopes']
        names=root.findall(NS+'name')
        assert [n.text or '' for n in names]==([] if expected['name'] is None else [expected['name']])
        if expected['coordinates'] is None:
            assert all(not any(k in a.attrib for k in ['x3','y3','z3']) for a in atoms)
        else:
            actual=[float(a.attrib[k]) for a in atoms for k in ['x3','y3','z3']]
            assert [struct.pack('>d',v) for v in actual]==[struct.pack('>d',v) for row in expected['coordinates'] for v in row]
        ids={a.attrib['id']:i for i,a in enumerate(atoms)}
        orders={'1':1,'2':2,'3':3,'A':1.5}
        bonds=[[*(ids[r] for r in b.attrib['atomRefs2'].split()),orders[b.attrib['order']]] for b in root.findall('./'+NS+'bondArray/'+NS+'bond')]
        assert bonds==expected['bonds']


@pytest.mark.parametrize('case',PACKET['invalid'],ids=lambda c:c['id'])
def test_unsupported_or_malformed_source_does_not_return_a_candidate(case):
    for call in [chematic.project_cjson_to_cml,lambda s:chematic.convert_format(s,'cjson','cml')]:
        with pytest.raises(ValueError):call(case['input'])


@pytest.mark.parametrize('override',[dict(coords=[]),dict(charges=[]),dict(name='new'),dict(comment='new')])
def test_overrides_cannot_silently_replace_source_information(override):
    with pytest.raises(ValueError,match='does not accept overrides'):
        chematic.convert_format(PACKET['valid'][0]['input'],'cjson','cml',**override)


def test_direct_input_utf8_and_atom_limits():
    with pytest.raises(UnicodeDecodeError):chematic.project_cjson_to_cml(b'\xff')
    with pytest.raises(ValueError,match='input bytes'):chematic.project_cjson_to_cml(' '*(16*1024*1024+1))
    for count in [10000,10001]:
        source=json.dumps({'atoms':{'elements':{'number':[6]*count}}})
        if count==10000:assert len(ET.fromstring(chematic.project_cjson_to_cml(source)).findall('./'+NS+'atomArray/'+NS+'atom'))==count
        else:
            with pytest.raises(ValueError,match='atom records'):chematic.project_cjson_to_cml(source)
