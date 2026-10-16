"""Construction-time custom charges, independent buffered-Coulomb controls."""
import math
import itertools
import pytest
import chematic


@pytest.mark.parametrize('variant', ['MMFF94', 'MMFF94s'])
def test_supplied_charges_change_energy_and_gradient(variant):
    mol = chematic.from_smiles('[Na+].[Cl-]')
    xyz = [[0,0,0], [4,0,0]]
    zero = mol.mmff_force_field(xyz, variant, ignore_interfrag_interactions=False, partial_charges=[0,0])
    q = [0.3,-0.8]
    custom = mol.mmff_force_field(xyz, variant, ignore_interfrag_interactions=False, partial_charges=q)
    q[:] = [0,0]  # constructor owns its prepared coefficients
    for epsilon in [1,4,78,0.5,1]:
        custom.set_dielectric_constant_(epsilon)
        e,g = custom.energy_gradient()
        ze,zg = zero.energy_gradient()
        assert e-ze == pytest.approx(332.0716*(0.3*-0.8/epsilon)/4.05, abs=1e-12)
        dx = 332.0716*(0.3*-0.8/epsilon)/4.05**2
        assert [v-w for a,b in zip(g,zg) for v,w in zip(a,b)] == pytest.approx([dx,0,0,-dx,0,0], abs=1e-12)
    assert mol.formula == 'ClNa'


@pytest.mark.parametrize('charges', [[], [1], [1,2,3], [math.nan,0], [math.inf,0], [-math.inf,0], [1e300,1e300], [1e154,1e154], [1e-300,1e-300]])
def test_invalid_charges_raise(charges):
    mol = chematic.from_smiles('[Na+].[Cl-]')
    with pytest.raises(ValueError, match='partial charge'):
        mol.mmff_force_field([[0,0,0], [4,0,0]], ignore_interfrag_interactions=False, partial_charges=charges)


@pytest.mark.parametrize('charges', [['x',0], [None,0], 1])
def test_wrong_charge_types_raise(charges):
    with pytest.raises(TypeError):
        chematic.from_smiles('[Na+].[Cl-]').mmff_force_field([[0,0,0],[4,0,0]], partial_charges=charges)


def test_none_keeps_default_and_failure_does_not_mutate_other_fields():
    mol = chematic.from_smiles('[Na+].[Cl-]')
    xyz = [[0,0,0], [4,0,0]]
    default = mol.mmff_force_field(xyz, ignore_interfrag_interactions=False)
    before = default.energy_gradient()
    assert mol.mmff_force_field(xyz, ignore_interfrag_interactions=False, partial_charges=None).energy_gradient() == before
    with pytest.raises(ValueError):
        mol.mmff_force_field(xyz, partial_charges=[math.nan,0])
    assert default.energy_gradient() == before
    assert default.positions == xyz


def test_charge_iterables_are_bounded_and_preserve_source_errors():
    mol = chematic.from_smiles('[Na+].[Cl-]')
    xyz = [[0,0,0], [4,0,0]]
    consumed = []
    def endless():
        for index in itertools.count():
            consumed.append(index)
            yield 0.3
    with pytest.raises(ValueError, match='one entry per atom'):
        mol.mmff_force_field(xyz, partial_charges=endless())
    assert consumed == [0,1,2]
    def broken():
        yield 0.3
        raise RuntimeError('charge source failed')
    with pytest.raises(RuntimeError, match='charge source failed'):
        mol.mmff_force_field(xyz, partial_charges=broken())
    values = (q for q in [0.3,-0.8])
    assert mol.mmff_force_field(xyz, ignore_interfrag_interactions=False, partial_charges=values).energy_gradient() == mol.mmff_force_field(xyz, ignore_interfrag_interactions=False, partial_charges=[0.3,-0.8]).energy_gradient()


def test_explicit_hydrogen_count_and_initially_zero_pairs():
    mol = chematic.from_smiles('C.C').add_hydrogens()
    tetrahedron = [[.63,.63,.63], [-.63,-.63,.63], [-.63,.63,-.63], [.63,-.63,-.63]]
    xyz = [[0,0,0], [4,0,0]] + tetrahedron + [[x+4,y,z] for x,y,z in tetrahedron]
    with pytest.raises(ValueError, match='one entry per atom'):
        mol.mmff_force_field(xyz, partial_charges=[0.5,-0.5])
    zero = mol.mmff_force_field(xyz, ignore_interfrag_interactions=False, partial_charges=[0]*mol.atom_count)
    for value in [0.5,1e-12]:
        q = [value,-value]+[0]*(mol.atom_count-2)
        custom = mol.mmff_force_field(xyz, ignore_interfrag_interactions=False, partial_charges=q)
        # Subtracting vdW can hide the tiny term, so only use finite scaling
        # here; the Rust electrostatics-only test resolves the 1e-12 charge.
        if value == 0.5:
            assert custom.energy()-zero.energy() == pytest.approx(-332.0716*value*value/4.05, abs=1e-12)
        assert math.isfinite(custom.energy())


def test_pair_selection_is_fixed_when_coordinates_move():
    mol = chematic.from_smiles('[Na+].[Cl-]')
    custom = mol.mmff_force_field([[0,0,0],[4,0,0]], non_bonded_thresh=3, ignore_interfrag_interactions=False, partial_charges=[0.5,-0.5])
    assert custom.energy() == 0
    custom.set_positions_([[0,0,0],[2,0,0]])
    assert custom.energy() == 0


def test_custom_charge_field_minimization_uses_same_model():
    mol = chematic.from_smiles('[Na+].[Cl-]')
    q = [0.2,-0.2]
    custom = mol.mmff_force_field([[0,0,0],[4,0,0]], ignore_interfrag_interactions=False, partial_charges=q)
    custom.set_dielectric_constant_(4)
    status, energy = custom.minimize_(10)
    fresh = mol.mmff_force_field(custom.positions, ignore_interfrag_interactions=False, partial_charges=q).set_dielectric_constant_(4)
    assert status in (0,1)
    assert energy == pytest.approx(fresh.energy(), abs=1e-9)
    assert custom.energy_gradient() == fresh.energy_gradient()
