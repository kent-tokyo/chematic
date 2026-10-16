"""Native explicit graph reader: no geometry perception or H removal."""
from pathlib import Path

import pytest
import chematic

TEXT = (Path(__file__).resolve().parents[3] / 'validation/rdkit_issues/pdb-conect/sulfur.pdb').read_text(encoding='utf-8')


def bond_rows(mol):
    block = mol.to_mol_v3000([[0.0, 0.0]] * mol.atom_count)
    section = block.split('M  V30 BEGIN BOND\n')[1].split('M  V30 END BOND')[0]
    return [tuple(map(int, row.split()[3:6])) for row in section.splitlines()]


def test_graph_coordinates_serials_and_charges():
    for text in (TEXT, TEXT.replace('\n', '\r\n')):
        frame = chematic.from_pdb_conect(text)
        assert set(frame) == {'mol', 'coords', 'serials'}
        assert frame['mol'].atom_count == 6
        assert frame['serials'] == [10, 20, 30, 40, 50, 90]
        assert frame['coords'] == [[0, 0, 0], [1.4, 0, 0], [-1.4, 0, 0], [0, 1.4, 0], [0, -1.4, 0], [0, 0, .8]]
        assert frame['mol'].formal_charge_per_atom() == [0, 0, 0, -1, -1, 0]
        assert bond_rows(frame['mol']) == [(2, 1, 2), (2, 1, 3), (1, 1, 4), (1, 1, 5)]


def test_extended_slots_require_keyword_opt_in():
    extended = TEXT.replace('CONECT   10   20   20   30   30\nCONECT   10   40   50', 'CONECT   10   20   20   30   30   40   50')
    with pytest.raises(ValueError, match='extended CONECT'):
        chematic.from_pdb_conect(extended)
    assert bond_rows(chematic.from_pdb_conect(extended, allow_extended_targets=True)['mol']) == bond_rows(chematic.from_pdb_conect(TEXT)['mol'])
    with pytest.raises(TypeError):
        chematic.from_pdb_conect(extended, True)


@pytest.mark.parametrize('text', [
    TEXT.replace('HETATM   20', 'HETATM   10', 1),
    TEXT.replace('CONECT   20   10   10', 'CONECT   20   10'),
    TEXT.replace('CONECT   50   10', 'CONECT   50   99'),
    TEXT.replace('CONECT   50   10', 'CONECT   50   50'),
    'MODEL        1\n' + TEXT, 'END\n' + TEXT,
    TEXT[:16] + 'A' + TEXT[17:], TEXT[:76] + '  ' + TEXT[78:],
    TEXT[:78] + '+1' + TEXT[80:], TEXT[:30] + '     NaN' + TEXT[38:],
    TEXT + 'CONECT   10   20\n',
    'ATOM\n',
])
def test_unsupported_or_malformed_inputs_raise(text):
    with pytest.raises(ValueError):
        chematic.from_pdb_conect(text)


def test_no_conect_and_empty_preserve_absence_of_edges():
    frame = chematic.from_pdb_conect('\n'.join(s for s in TEXT.splitlines() if s.startswith('HETATM')))
    assert frame['mol'].atom_count == 6
    assert 'M  V30 COUNTS 6 0' in frame['mol'].to_mol_v3000([[0, 0]] * 6)
    empty = chematic.from_pdb_conect('')
    assert empty['mol'].atom_count == 0 and empty['coords'] == empty['serials'] == []
