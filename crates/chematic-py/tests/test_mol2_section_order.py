"""MOL2 core sections precede charge metadata without suppressing charges."""
import json
from pathlib import Path
import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
FOLDER = ROOT / 'validation/openbabel_issues/mol2-section-order'
CASES = json.loads((FOLDER / 'cases.json').read_text())['cases']

@pytest.mark.parametrize('case', [c for c in CASES if 'smiles' in c], ids=lambda c: c['id'])
def test_atom_bond_then_charge_metadata_and_roundtrip(case):
    mol = chematic.from_smiles(case['smiles'])
    output = mol.to_mol2(atom_bond_adjacent=True)
    default = mol.to_mol2()
    assert default == mol.to_mol2(atom_bond_adjacent=False)
    if any(case['formal_charges']):
        assert default.index('@<TRIPOS>UNITY_ATOM_ATTR') < default.index('@<TRIPOS>BOND')
    markers = [line for line in output.splitlines() if line.startswith('@<TRIPOS>')]
    expected = ['@<TRIPOS>MOLECULE','@<TRIPOS>ATOM','@<TRIPOS>BOND']
    if any(case['formal_charges']):
        expected.append('@<TRIPOS>UNITY_ATOM_ATTR')
    assert markers == expected
    record = json.loads(chematic.parse_mol2_record_json(output))
    assert record['formal_charges'] == case['formal_charges']
    assert len(record['atoms']) == mol.atom_count
    assert all(atom['partial_charge'] == 0 for atom in record['atoms'])

@pytest.mark.parametrize('crlf', [False, True])
def test_source_preserving_apis_keep_original_section_order(crlf):
    source = (FOLDER / 'charged-source.mol2').read_text()
    if crlf:
        source = source.replace('\n','\r\n')
    assert source.index('@<TRIPOS>UNITY_ATOM_ATTR') < source.index('@<TRIPOS>BOND')
    assert chematic.roundtrip_mol2_record(source) == source
    assert chematic.convert_format(source,'mol2','mol2') == source


def test_neutral_output_bytes_remain_the_baseline():
    assert chematic.from_smiles('CO').to_mol2() == (FOLDER / 'neutral-expected.mol2').read_text()
