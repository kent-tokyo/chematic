"""Represented atom fields never silently default or come from quoted text."""
import json
from pathlib import Path
import chematic
import pytest

PACKET=json.loads((Path(__file__).resolve().parents[3]/'validation/v3000_atom_scalars/cases.json').read_text(encoding="utf-8"))


@pytest.mark.parametrize('case',PACKET['cases'],ids=lambda c:c['id'])
def test_scalar_identity_and_source_retention(case):
    for source in [case['input'],case['input'].encode()]:
        mol=chematic.from_mol_v3000(source)
        atom=json.loads(mol.to_moljson())['atoms'][0]
        e=case['expected']
        assert atom['element']=={1:'H',6:'C',7:'N',8:'O'}[e['number']]
        assert atom['charge']==e['charge']
        assert (atom['isotope'] or 0)==e['isotope']
        if e['atom_map']:
            assert f":{e['atom_map']}]" in mol.smiles
        else:
            assert ':' not in mol.smiles
        assert chematic.roundtrip_mol_v3000_source(source)==case['input']
    assert chematic.convert_format(case['input'],'v3000','smiles')==mol.smiles
    assert chematic.convert_format(case['input'],'v3000','v3000')==case['input']


@pytest.mark.parametrize('case',PACKET['invalid'],ids=lambda c:c['id'])
def test_bad_fields_fail_on_direct_common_and_source_routes(case):
    for call in [chematic.from_mol_v3000,chematic.from_mol_v3000_with_coords_3d,
                 chematic.roundtrip_mol_v3000_source,
                 lambda s:chematic.convert_format(s,'v3000','smiles'),
                 lambda s:chematic.convert_format(s,'v3000','v3000')]:
        with pytest.raises(ValueError,match=case['field']):call(case['input'])
