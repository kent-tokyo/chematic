"""The spec-authored oracle replaces each complete legacy atom-block channel."""
import io
import json
from pathlib import Path
import re

import chematic
import pytest

PACKET = json.loads((Path(__file__).parents[3] / 'validation/v2000_property_precedence/cases.json').read_text(encoding="utf-8"))
NUMBERS = dict(C=6, N=7, O=8)


def assert_fields(mol, expected):
    data = json.loads(mol.to_moljson())
    atoms = data['atoms']
    actual = sorted((NUMBERS[a['element']], a['charge'], a['isotope'] or 0) for a in atoms)
    wanted = sorted(zip(expected['numbers'], expected['charges'], expected['isotopes']))
    assert actual == wanted
    ids = {a['id']: NUMBERS[a['element']] for a in atoms}
    bonds = sorted([*sorted([ids[b['source_id']], ids[b['target_id']]]), b['order']] for b in data['bonds'])
    assert bonds == sorted([*sorted([expected['numbers'][a], expected['numbers'][b]]), order] for a,b,order in expected['bonds'])
    assert sorted(map(int, re.findall(r':(\d+)\]', mol.smiles))) == sorted(m for m in expected['atom_maps'] if m)


@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c: c['id'])
def test_precedence_in_mol_common_sdf_and_fast_supplier(case, tmp_path):
    for source in [case['input'], case['input'].encode()]:
        assert_fields(chematic.from_mol_block(source), case['expected'])
        assert_fields(chematic.from_mol_block_with_coords(source)[0], case['expected'])
    smiles = chematic.convert_format(case['input'], 'mol', 'smiles')
    assert_fields(chematic.from_smiles(smiles), case['expected'])
    sdf = case['input'] + '$$$$\n'
    records = list(chematic.iter_sdf_str(sdf))
    assert len(records) == 1
    assert_fields(records[0].mol, case['expected'])
    path = tmp_path / 'input.sdf'
    path.write_text(sdf)
    fast = list(chematic.SDMolSupplier(path))
    assert len(fast) == 1
    assert_fields(fast[0], case['expected'])
    forward = list(chematic.ForwardSDMolSupplier(io.BytesIO(sdf.encode())))
    assert len(forward) == 1
    assert_fields(forward[0], case['expected'])


def test_property_state_resets_per_record_and_sd_values_are_opaque(tmp_path):
    cases = PACKET['cases'] + list(reversed(PACKET['cases']))
    sdf = ''.join(c['input'] + '> <note>\nM  CHG  1   1   3\nM  ISO  1   1  99\n\n$$$$\n' for c in cases)
    path = tmp_path / 'many.sdf'
    path.write_text(sdf)
    for records in [list(chematic.SDMolSupplier(path)),
                    list(chematic.ForwardSDMolSupplier(io.BytesIO(sdf.encode())))]:
        assert len(records) == len(cases)
        for mol, case in zip(records, cases):
            assert_fields(mol, case['expected'])
