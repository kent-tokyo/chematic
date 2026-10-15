"""Native carbon enumeration budget and RDKit #9701 diagnostic contracts."""
import pytest
import chematic


@pytest.mark.parametrize('smiles,count',[
    ('NC12CC3CC(CC(C3)C1)C2',4),
    ('OC1C(O)C(O)C(O)C(O)C1O',9),
    ('CCCCCCCCCCCCCCCC',1),
    ('CCCCCCCCCCCCCCCC(O)F',2),
    ('FC(F)(Br)Cl',1),
    ('CC(C)(Br)Cl',1),
    ('[H]C([H])(F)Cl',1),
    ('[H]C([2H])(F)Cl',2),
    ('[2H]C(F)Cl',2),
    ('N[C@@H](C)C(=O)O',1),
    ('CC(O)C(O)C',3),
])
def test_expected_unique_assignments(smiles,count):
    mol=chematic.from_smiles(smiles)
    forms=mol.enumerate_stereoisomers()
    assert len(forms)==count
    keys={f.smiles for f in forms}
    assert len(keys)==count
    assert all(chematic.from_smiles(key).smiles==key for key in keys)


def test_explicit_h_aminoadamantane_retains_four_assignments():
    mol=chematic.from_smiles('NC12CC3CC(CC(C3)C1)C2').add_hydrogens()
    forms=mol.enumerate_stereoisomers()
    assert len(forms)==4
    assert all(f.atom_count==mol.atom_count for f in forms)


def test_name_and_typed_properties_are_independent_on_all_forms():
    mol=chematic.from_smiles('FC(Cl)(Br)I')
    mol.SetName('source')
    mol.SetIntProp('integer',-3)
    mol.SetBoolProp('bool',True)
    mol.SetDoubleProp('double',1.25)
    forms=mol.enumerate_stereoisomers()
    assert len(forms)==2
    for f in forms:
        assert f.GetName()=='source'
        assert f.GetIntProp('integer')==-3
        assert f.GetBoolProp('bool') is True
        assert f.GetDoubleProp('double')==1.25
    forms[0].SetName('changed')
    forms[0].SetIntProp('integer',9)
    assert mol.GetName()==forms[1].GetName()=='source'
    assert mol.GetIntProp('integer')==forms[1].GetIntProp('integer')==-3


def test_no_candidate_copy_retains_properties():
    mol=chematic.from_smiles('CCCCCCCCCCCC')
    mol.SetName('chain')
    forms=mol.enumerate_stereoisomers()
    assert len(forms)==1
    assert forms[0].GetName()=='chain'
    assert forms[0].smiles==mol.smiles


def test_actual_seven_candidate_limit_still_applies():
    mol=chematic.from_smiles('.'.join(['FC(Cl)(Br)I']*7))
    assert mol.enumerate_stereoisomers()==[]


def test_six_mapped_candidate_sites_have_64_assignments():
    smiles='.'.join(f'F[CH:{i}](Cl)Br' for i in range(1,7))
    assert len(chematic.from_smiles(smiles).enumerate_stereoisomers())==64


def test_existing_specified_configuration_is_preserved_with_new_center():
    mol=chematic.from_smiles('N[C@@H](C)C(O)F')
    forms=mol.enumerate_stereoisomers()
    assert len(forms)==2
    # The stereospecific query covers the fixed source center.
    query=chematic.compile_smarts('N[C@@H](C)C')
    assert all(query.matches(f,use_chirality=True) for f in forms)
