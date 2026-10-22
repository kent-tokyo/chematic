import json
from pathlib import Path
import pytest
import chematic
PACKET=json.loads((Path(__file__).parents[3]/'validation/sdf_source_document/cases.json').read_text(encoding='utf-8'))
@pytest.mark.parametrize('case',PACKET['cases'],ids=lambda c:c['id'])
@pytest.mark.parametrize('layout',PACKET['variants'])
def test_exact_document_and_single_record_gate(case,layout):
    text=case['input'].replace('\n','\r\n') if layout=='crlf' else case['input'].rstrip('\n') if layout=='no-final-newline' else case['input']
    assert chematic.roundtrip_sdf_source(text.encode())==text
    assert chematic.convert_format(text,'.SDF','sdf')==text
    if case['records']>1:
        for source in ['mol','sdf']:
            for target in ['smiles','cjson','mol','mol_v3000']:
                with pytest.raises(ValueError,match='records.*limit 1'):
                    chematic.convert_format(text,source,target)
        with pytest.raises(ValueError):chematic.single_sdf_mol_block(text)
    else:
        block=chematic.single_sdf_mol_block(text)
        assert '<NOTE>' not in block
        assert chematic.convert_format(text,'sdf','smiles')=='O'
@pytest.mark.parametrize('kwargs',[dict(coords=[]),dict(charges=[]),dict(name='new'),dict(comment='new')])
def test_same_format_overrides_are_refused(kwargs):
    with pytest.raises(ValueError,match='does not accept overrides'):
        chematic.convert_format(PACKET['cases'][0]['input'],'sdf','sdf',**kwargs)
def test_sdf_output_is_framed_and_single_input_fields_do_not_confuse_the_parser():
    output=chematic.convert_format('CO','smiles','sdf')
    assert output.endswith('$$$$\n')
    assert chematic.convert_format(output,'sdf','smiles')
    source=PACKET['cases'][0]['input'].split('$$$$\n',1)[0]+'$$$$\n'
    assert '<ID>' not in chematic.single_sdf_mol_block(source)
    assert chematic.convert_format(source,'sdf','smiles')=='C'
    assert chematic.convert_format(source,'mol','smiles')=='C'
def test_no_valid_prefix_is_returned_on_later_failure():
    source=PACKET['cases'][0]['input']
    for bad in ['', ' \n', '$$$$\n', '$$$$\n'+source, source+'junk\n', source.replace('M  END','M  END suffix')]:
        for call in [chematic.roundtrip_sdf_source,lambda text:chematic.convert_format(text,'sdf','sdf')]:
            with pytest.raises(ValueError):call(bad)
def test_utf8_and_input_cap():
    with pytest.raises(UnicodeDecodeError):chematic.roundtrip_sdf_source(b'\xff')
    with pytest.raises(ValueError,match='input bytes'):chematic.roundtrip_sdf_source(' '*(16*1024*1024+1))
def test_zero_atom_records_are_retained_and_counted():
    text='zero\n\n\n  0  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n>  <NOTE>\nmetadata only\n\n$$$$\n'
    assert chematic.convert_format(text,'sdf','sdf')==text
    with pytest.raises(ValueError,match='records.*limit 1'):
        chematic.convert_format(text+text,'sdf','smiles')
def test_version_detection_uses_physical_lines_in_utf8_header_text():
    text=PACKET['cases'][2]['input'].replace('original high precision','original\u2028high precision')
    assert chematic.convert_format(text,'sdf','smiles')=='O'
