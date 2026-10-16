import json
from pathlib import Path
import pytest
import chematic

DATA=json.loads((Path(__file__).resolve().parents[3]/'validation/rdkit_issues/chemdraw-block/cases.json').read_text(encoding='utf-8'))

@pytest.mark.parametrize('case',DATA['cases'],ids=lambda c:c['id'])
def test_explicit_binary_and_xml_modes(case):
    source=bytes.fromhex(case['cdx_hex'])
    for mode in ('AUTO','CDX','cdx'):
        result=chematic.from_chemdraw_block(source,format=mode,allow_metadata_loss=case['allow_metadata_loss'])
        assert result['format']=='CDX' and result['source_bytes']==source
        assert len(result['omissions'])==int(case['allow_metadata_loss'])
        assert all(x['tag']==0x0805 for x in result['omissions'])
        assert len(result['molecules'])==len(result['coordinates'])
        assert '<CDXML' in result['projected_cdxml']
        if case['id'].startswith('aromatic'):
            assert result['molecules'][0].smiles=='c1ccccc1'
        if case['id'].startswith('ammonium'):
            assert result['molecules'][0].formula=='H4N+'
        xml_result=chematic.from_chemdraw_block(case['cdxml'].encode('utf-8'),format='CDXML')
        assert result['coordinates']==xml_result['coordinates']
        assert [m.atom_count for m in result['molecules']]==[m.atom_count for m in xml_result['molecules']]
        assert [m.formal_charge_per_atom() for m in result['molecules']]==[m.formal_charge_per_atom() for m in xml_result['molecules']]

def test_metadata_requires_explicit_consent():
    c=next(c for c in DATA['cases'] if c['allow_metadata_loss'])
    with pytest.raises(ValueError,match='requires opt-in'):
        chematic.from_chemdraw_block(bytes.fromhex(c['cdx_hex']))

@pytest.mark.parametrize('data,mode',[(b'not xml','AUTO'),(b'<CDXML/>','CDX'),(b'\xff','CDXML'),(b'','other'),(b'VjCD0100','AUTO')])
def test_wrong_format_and_malformed_input_raise(data,mode):
    with pytest.raises(ValueError):chematic.from_chemdraw_block(data,format=mode)

def test_every_truncation_and_trailing_data_raise():
    data=bytes.fromhex(DATA['cases'][0]['cdx_hex'])
    for i in range(len(data)):
        with pytest.raises(ValueError):chematic.from_chemdraw_block(data[:i],format='CDX')
    with pytest.raises(ValueError):chematic.from_chemdraw_block(data+b'\0')

def test_bytes_and_keyword_contract():
    with pytest.raises(TypeError):chematic.from_chemdraw_block('<CDXML/>')
    with pytest.raises(TypeError):chematic.from_chemdraw_block(b'<CDXML/>','AUTO')
    assert chematic.from_chemdraw_block(b'<CDXML><page/></CDXML>')['molecules']==[]
    bom=b'\xef\xbb\xbf<CDXML><page/></CDXML>'
    assert chematic.from_chemdraw_block(bom)['source_bytes']==bom

def test_stub_declares_new_binary_entrypoint():
    text=(Path(__file__).resolve().parents[1]/'python/chematic/__init__.pyi').read_text(encoding='utf-8')
    assert 'def from_chemdraw_block(data: bytes,' in text
