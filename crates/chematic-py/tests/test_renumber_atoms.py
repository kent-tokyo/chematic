"""Atom permutation contracts, including RDKit issue #9629."""
import itertools
import pytest
import chematic
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import rdMolDescriptors as D


@pytest.mark.parametrize('smiles', ['C1CCN2CCCC2C1','C1CCN2CCC2C1','C1CN2CCC1CC2','C1CC2CCC1N2'])
def test_bicyclic_counts_are_stable_with_warm_and_cold_caches(smiles):
    mol=Chem.MolFromSmiles(smiles)
    before=D.CalcNumAtomStereoCenters(mol)
    n=mol.GetNumAtoms()
    mol.GetRingInfo().AtomRings()
    mol.GetAtoms()
    mol.GetBonds()
    for order in [list(range(n)),list(reversed(range(n))),list(range(1,n))+[0]]:
        out=Chem.RenumberAtoms(mol,order)
        assert D.CalcNumAtomStereoCenters(out)==before
        assert Chem.MolToSmiles(out)==Chem.MolToSmiles(mol)
        assert out._atom_cache is None
        assert out._ring_info_cache is None
    if smiles=='C1CCN2CCCC2C1':
        assert before==1


def test_all_tetrahedral_permutations_keep_configuration_and_indexed_atoms():
    mol=chematic.from_smiles('F[C@](Cl)(Br)I')
    rows=mol.atom_table
    for order in itertools.permutations(range(5)):
        out=mol.renumber_atoms(order)
        assert out.smiles==mol.smiles
        assert out.num_stereocenters==1
        # Atom table's final degree/ring fields are graph-derived; each follows
        # the same source atom after a permutation.
        assert out.atom_table==[rows[i] for i in order]


@pytest.mark.parametrize('smiles', [
    'N[C@@H:7](C)C(=O)O','[H][C@](F)(Cl)Br','[2H][C@](F)(Cl)Br',
    'F/C=C/Cl','F/C=C\\Cl','CC=[C@AL1]=C(C)F',
    'Cl[Pt@SP1](Cl)(N)N','F[P@TB1](Cl)(Br)(I)N','F[Co@OH1](Cl)(Br)(I)(N)O',
])
def test_stereo_classes_and_isotopes_survive_reverse(smiles):
    mol=chematic.from_smiles(smiles)
    out=mol.renumber_atoms(list(reversed(range(mol.atom_count))))
    assert out.smiles==mol.smiles
    assert out.atom_count==mol.atom_count
    assert out.exact_mass==pytest.approx(mol.exact_mass)


def test_typed_properties_and_name_are_copied_independently():
    source=Chem.MolFromSmiles('N[C@@H](C)C(=O)O')
    source.SetName('source')
    source.SetIntProp('integer',-7)
    source.SetUnsignedProp('unsigned',9)
    source.SetDoubleProp('double',1.25)
    source.SetBoolProp('bool',True)
    source.SetProp('string','unchanged')
    out=Chem.RenumberAtoms(source,list(reversed(range(source.GetNumAtoms()))))
    assert out.GetName()=='source'
    assert out.GetIntProp('integer')==-7
    assert out.GetUnsignedProp('unsigned')==9
    assert out.GetDoubleProp('double')==1.25
    assert out.GetBoolProp('bool') is True
    assert out.GetProp('string')=='unchanged'
    out.SetName('copy')
    out.SetIntProp('integer',3)
    assert source.GetName()=='source'
    assert source.GetIntProp('integer')==-7


def test_get_num_atoms_accounts_for_explicit_h_in_permutations():
    explicit=Chem.MolFromSmiles('[H][C@](F)(Cl)Br')
    assert explicit.GetNumAtoms()==5
    assert explicit.GetNumAtoms(onlyHeavy=True)==4
    assert Chem.RenumberAtoms(explicit,list(range(5))).GetNumAtoms()==5
    mol=Chem.AddHs(Chem.MolFromSmiles('CO'))
    assert mol.GetNumAtoms()==6
    assert mol.GetNumAtoms(onlyHeavy=True)==2
    assert Chem.RenumberAtoms(mol,list(reversed(range(mol.GetNumAtoms())))).GetNumAtoms()==6


@pytest.mark.parametrize('order', [[],[0],[0,1],[0,1,1],[0,1,3],[-1,1,2],[0,1,2**32],[False,1,2],[0,1,1.5],[0,1,'2']])
def test_invalid_permutations_raise_value_error(order):
    with pytest.raises(ValueError):
        chematic.from_smiles('CCO').renumber_atoms(order)


def test_empty_graph_accepts_empty_permutation():
    mol=Chem.RWMol().GetMol()
    assert Chem.RenumberAtoms(mol,[]).GetNumAtoms()==0


def test_bond_indices_and_endpoint_direction_remain_fixed():
    mol=chematic.from_smiles('N->[Cu].F/C=C/Cl')
    order=list(reversed(range(mol.atom_count)))
    inverse={old:new for new,old in enumerate(order)}
    out=mol.renumber_atoms(order)
    for before,after in zip(mol.bond_table,out.bond_table):
        assert after[:2]==(inverse[before[0]],inverse[before[1]])
        assert after[2:]==before[2:]
