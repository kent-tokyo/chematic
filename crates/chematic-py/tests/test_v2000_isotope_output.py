"""Explicit isotope labels survive regenerated MOL/SDF, including mass difference 0."""
import json
from pathlib import Path
import re

import chematic
import pytest

PACKET = json.loads((Path(__file__).parents[3] / 'validation/v2000_isotope_output/cases.json').read_text(encoding="utf-8"))
NUMBERS = dict(H=1, C=6, N=7, O=8, Cl=17, Br=35, U=92)


def assert_fields(mol, expected):
    value = json.loads(mol.to_moljson())
    atoms = value['atoms']
    assert [NUMBERS[a['element']] for a in atoms] == expected['numbers']
    assert [a['charge'] for a in atoms] == expected['charges']
    assert [a['isotope'] or 0 for a in atoms] == expected['isotopes']
    ids = {a['id']: i for i, a in enumerate(atoms)}
    bonds = sorted([*sorted([ids[b['source_id']], ids[b['target_id']]]), b['order']] for b in value['bonds'])
    assert bonds == expected['bonds']
    assert sorted(map(int, re.findall(r':(\d+)\]', mol.smiles))) == sorted(m for m in expected['atom_maps'] if m)


@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c: c['id'])
def test_isotope_output_in_direct_coordinate_common_and_sdf_routes(case):
    mol = chematic.from_smiles(case['smiles'])
    coords = [[float(i), 0.0] for i in range(len(case['expected']['numbers']))]
    outputs = [mol.to_mol_block(), mol.to_mol_block(strict=True),
               mol.to_mol_block_2d(coords),
               chematic.convert_format(case['smiles'], 'smiles', 'mol'),
               chematic.convert_format(case['smiles'], 'smiles', 'sdf'),
               mol.to_sdf_with_charges([xy + [0.5] for xy in coords], [0.0] * len(coords))]
    for output in outputs:
        assert_fields(chematic.from_mol_block(output), case['expected'])
        if '$$$$' in output:
            records = list(chematic.iter_sdf_str(output))
            assert len(records) == 1
            assert_fields(records[0].mol, case['expected'])


def test_fixed_width_and_historical_spaced_iso_records_are_accepted():
    block = chematic.from_smiles('[9999CH4]').to_mol_block()
    assert 'M  ISO  1   19999\n' in block
    for text in [block, block.replace('M  ISO  1   19999', 'M  ISO 1 1 9999')]:
        assert json.loads(chematic.from_mol_block(text).to_moljson())['atoms'][0]['isotope'] == 9999


def test_partial_charge_sdf_keeps_supplied_geometry_and_data():
    mol = chematic.from_smiles('[14NH3+:7][12CH3:9]')
    coords = [[0.125, -0.25, 0.5], [1.0, 2.0, 3.0]]
    output = mol.to_sdf_with_charges(coords, [-0.125, 0.125])
    assert output.splitlines()[1][20:22] == '3D'
    for line, point in zip(output.splitlines()[4:6], coords):
        assert [float(line[i:i + 10]) for i in [0, 10, 20]] == point
    records = list(chematic.iter_sdf_str(output))
    assert len(records) == 1
    assert records[0].get('PARTIAL_CHARGES') == '-0.1250 0.1250'
    assert output.endswith('$$$$\n')
