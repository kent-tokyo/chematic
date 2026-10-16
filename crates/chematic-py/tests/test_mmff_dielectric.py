"""Prepared MMFF fields update energy/gradient without preparing again."""
import json
from pathlib import Path

import pytest
import chematic

DATA = json.loads((Path(__file__).resolve().parents[3] / 'validation/openbabel_issues/mmff-dielectric/cases.json').read_text(encoding='utf-8'))


@pytest.mark.parametrize('row', DATA['cases'], ids=lambda r:r['smiles'])
@pytest.mark.parametrize('variant', ['MMFF94', 'MMFF94s'])
def test_live_updates_match_pinned_rdkit(row, variant):
    mol = chematic.from_smiles(row['smiles']).add_hydrogens()
    field = mol.mmff_force_field(row['coords'], variant, ignore_interfrag_interactions=False)
    assert field.dielectric_constant == 1.0
    before = field.energy_gradient()
    for epsilon in (1.0, 4.0, 78.0, .5, 4.0, 1.0):
        assert field.set_dielectric_constant_(epsilon) is field
        expected = next(r for r in row['references'] if r['variant'] == variant and r['epsilon'] == epsilon)
        energy, gradient = field.energy_gradient()
        assert field.dielectric_constant == epsilon
        assert field.positions == row['coords']
        assert energy == pytest.approx(expected['energy'], abs=1e-9, rel=1e-11)
        assert [v for p in gradient for v in p] == pytest.approx(expected['gradient'], abs=1e-9, rel=1e-11)
    assert field.energy_gradient() == before


@pytest.mark.parametrize('value', [0., -0., -1., float('nan'), float('inf'), -float('inf'), 5e-324])
def test_invalid_updates_preserve_state(value):
    row = DATA['cases'][0]
    mol = chematic.from_smiles(row['smiles'])
    field = mol.mmff_force_field(row['coords'], ignore_interfrag_interactions=False)
    field.set_dielectric_constant_(4.)
    before = field.energy_gradient()
    with pytest.raises(ValueError, match='MMFF dielectric'):
        field.set_dielectric_constant_(value)
    assert field.dielectric_constant == 4.
    assert field.energy_gradient() == before
    assert field.positions == row['coords']


def test_uff_refuses_mmff_only_screening():
    mol = chematic.from_smiles('CCO').add_hydrogens()
    row = DATA['cases'][2]
    field = mol.uff_force_field(row['coords'])
    before = field.energy_gradient()
    with pytest.raises(ValueError, match='UFF'):
        field.set_dielectric_constant_(4.)
    with pytest.raises(ValueError, match='UFF'):
        _ = field.dielectric_constant
    assert field.energy_gradient() == before
    assert field.positions == row['coords']


@pytest.mark.parametrize('variant', ['MMFF94', 'MMFF94s'])
def test_minimize_and_coordinate_update_use_current_dielectric(variant):
    row = DATA['cases'][4]
    mol = chematic.from_smiles(row['smiles']).add_hydrogens()
    field = mol.mmff_force_field(row['coords'], variant)
    _ = field.energy_gradient()
    field.set_dielectric_constant_(4.).set_dielectric_constant_(78.)
    moved = list(row['coords'][0]); moved[0] += .01
    field.set_position_(0, moved)
    fresh = mol.mmff_force_field(field.positions, variant)
    fresh.set_dielectric_constant_(78.)
    assert field.energy_gradient() == fresh.energy_gradient()
    assert field.minimize_(10) == fresh.minimize_(10)
    assert field.positions == fresh.positions
    assert field.dielectric_constant == 78.


def test_default_still_matches_existing_stateless_api():
    row = DATA['cases'][4]
    mol = chematic.from_smiles(row['smiles']).add_hydrogens()
    for variant in ('MMFF94', 'MMFF94s'):
        field = mol.mmff_force_field(row['coords'], variant)
        assert field.energy() == mol.rdkit_mmff_energy(row['coords'], variant)


def test_stub_declares_mutable_dielectric():
    text = (Path(__file__).resolve().parents[1] / 'python/chematic/__init__.pyi').read_text(encoding='utf-8')
    assert 'def dielectric_constant(self) -> float:' in text
    assert 'def set_dielectric_constant_(self, value: float) -> MolecularForceField:' in text
