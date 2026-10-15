"""Original MIT OR Apache-2.0 RDFile record, source and boundary contracts."""
from pathlib import Path
import pytest
import chematic

ROOT=Path(__file__).resolve().parents[3]
TEXT=(ROOT/'validation/openbabel_issues/rdfile/mixed-records.rdf').read_text(encoding="utf-8")


def test_mixed_records_metadata_and_projection():
    d=chematic.parse_rdfile(TEXT)
    assert [r['kind'] for r in d['records']]==['molecule','molecule','reaction','molecule','reaction']
    assert d['records'][0]['internal_registry']=='0007'
    assert d['records'][0]['data']==[{'name':'note','value':'first'},{'name':'note','value':'second\nwith  two spaces\n\nlast'}]
    smiles=chematic.rdfile_to_smiles(TEXT)
    assert len(smiles)==5 and smiles[-2:]==[None,None]
    assert '13C' in smiles[1] and '>>' in smiles[2]
    assert ':1' in smiles[2] and ':2' in smiles[2]
    assert d['source_text']==TEXT


def test_read_write_read_retains_order_fields_and_structures():
    d=chematic.parse_rdfile(TEXT);back=chematic.parse_rdfile(chematic.write_rdfile(d))
    assert back['records']==d['records'] and back['date_comment']==d['date_comment']


def test_writer_uses_current_fields_and_preserves_original_sidecar():
    d=chematic.parse_rdfile(TEXT);d['records'][0]['data'][0]['value']='edited'
    assert chematic.parse_rdfile(chematic.write_rdfile(d))['records'][0]['data'][0]['value']=='edited'
    assert d['source_text']==TEXT


@pytest.mark.parametrize('n',[0,1,72,73,74,80,81,160,1000])
def test_long_ascii_data_roundtrip(n):
    d=chematic.parse_rdfile(TEXT);value='x'*n+'\n  next\n'
    d['records'][0]['data']=[{'name':'long','value':value}]
    text=chematic.write_rdfile(d);assert chematic.parse_rdfile(text)['records'][0]['data'][0]['value']==value


@pytest.mark.parametrize('suffix',['\n','\r\n',''])
def test_line_endings(suffix):
    text=TEXT.replace('\n','\r\n') if suffix=='\r\n' else TEXT.rstrip('\n')+suffix
    d=chematic.parse_rdfile(text);assert d['source_text']==text
    assert chematic.parse_rdfile(chematic.write_rdfile(d))['records']==d['records']


@pytest.mark.parametrize('kwargs',[{'max_input_bytes':1},{'max_records':4},{'max_fields_per_record':1},{'max_line_len':4},{'max_reaction_components':1}])
def test_limits(kwargs):
    with pytest.raises(ValueError):chematic.parse_rdfile(TEXT,**kwargs)


@pytest.mark.parametrize('text',['','$RDFILE 2\n$DATM\n','$RDFILE 1\n','$RDFILE 1\nwrong\n','$RDFILE 1\n$DATM\n\n','$RDFILE 1\n$DATM\n$MIREG -1\n',TEXT.replace('M  END','M  BAD',1),TEXT.replace('$RXN','$RXN V3000',1),TEXT.replace('$MOL','$UNKNOWN'),TEXT.replace('$MFMT $MIREG 0007','$MFMT $MIREG')])
def test_invalid_documents(text):
    with pytest.raises(ValueError):chematic.parse_rdfile(text)


@pytest.mark.parametrize('value',['x\n$MFMT','  leading whitespace','unicode 🧪','x\rbad','$MFMT'])
def test_writer_rejects_unsupported_data_and_injection(value):
    d=chematic.parse_rdfile(TEXT);d['records'][0]['data'][0]['value']=value
    with pytest.raises(ValueError):chematic.write_rdfile(d)


def test_writer_refuses_unknown_keys_and_ambiguous_registry():
    d=chematic.parse_rdfile(TEXT);d['records'][0]['unexpected']='discard me'
    with pytest.raises(ValueError):chematic.write_rdfile(d)
    d=chematic.parse_rdfile(TEXT);d['records'][0]['external_registry']='second'
    with pytest.raises(ValueError):chematic.write_rdfile(d)


def test_registry_only_records_keep_slots_without_lookup():
    text='$RDFILE 1\n$DATM\n$MIREG 1\n$MEREG name\n$RIREG 2\n$REREG rxn\n'
    assert chematic.rdfile_to_smiles(text)==[None]*4
    d=chematic.parse_rdfile(text);assert all(r['structure'] is None for r in d['records'])


def test_empty_library():
    assert chematic.rdfile_to_smiles('$RDFILE 1\n$DATM\n')==[]


def test_molecule_v3000_embedded_payload():
    mol=chematic.from_smiles('[NH4+]');block=mol.to_mol_v3000([[0,0]])
    text='$RDFILE 1\n$DATM\n$MFMT\n'+block
    d=chematic.parse_rdfile(text)
    assert d['records'][0]['structure']==block
    assert chematic.rdfile_to_smiles(text)==['[NH4+]']


def test_zero_agent_field_in_rxn_counts_is_supported():
    d=chematic.parse_rdfile(TEXT);r=d['records'][2];lines=r['structure'].splitlines();lines[4]+='  0';r['structure']='\n'.join(lines)+'\n'
    assert '>>' in chematic.rdfile_to_smiles(chematic.write_rdfile(d))[2]


def test_nonzero_rxn_agents_are_refused():
    d=chematic.parse_rdfile(TEXT);r=d['records'][2];lines=r['structure'].splitlines();lines[4]+='  1';r['structure']='\n'.join(lines)+'\n'
    with pytest.raises(ValueError):chematic.write_rdfile(d)


def test_unwrapped_interchange_profile_and_plain_long_input():
    d=chematic.parse_rdfile(TEXT);value='a'*320+'\nnext'
    d['records'][0]['data']=[{'name':'long','value':value}]
    text=chematic.write_rdfile(d,wrap_data=False)
    assert any(len(s)>81 for s in text.splitlines())
    assert chematic.parse_rdfile(text)['records'][0]['data'][0]['value']==value
