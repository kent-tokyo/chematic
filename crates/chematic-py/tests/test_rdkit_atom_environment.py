"""Original minimal atom environments; oracle/fixture hashes in validation."""
import ast
from pathlib import Path
import pytest
import chematic

@pytest.mark.parametrize("radius,enforce,bonds,distances",[
    (0,True,[],{0:0}),(1,True,[0],{0:0,1:1}),(2,True,[0,1],{0:0,1:1,2:2}),
    (3,True,[],{}),(3,False,[0,1],{0:0,1:1,2:2}),(2**32-1,False,[0,1],{0:0,1:1,2:2}),
])
def test_chain(radius,enforce,bonds,distances):
    assert chematic.from_smiles("CCO").rdkit_atom_environment(radius,0,enforce_size=enforce)=={"bond_indices":bonds,"atom_distances":distances}

@pytest.mark.parametrize("radius",[2,3])
def test_cyclic_pending_duplicates(radius):
    r=chematic.from_smiles("C1CCC1").rdkit_atom_environment(radius,0)
    assert r=={"bond_indices":[0,3,1,2],"atom_distances":{0:0,1:1,2:2,3:1}}

def test_h_root_filter_and_added_hs():
    m=chematic.from_smiles("[2H]C")
    assert m.rdkit_atom_environment(1,0)["bond_indices"]==[0]
    assert m.rdkit_atom_environment(1,1)["atom_distances"]=={}
    assert m.rdkit_atom_environment(1,1,use_hs=True)["bond_indices"]==[0]
    m=chematic.from_smiles("C").rdkit_add_hydrogens()
    assert m.rdkit_atom_environment(1,0)["bond_indices"]==[]
    assert m.rdkit_atom_environment(1,0,use_hs=True)["bond_indices"]==[0,1,2,3]

@pytest.mark.parametrize("kwargs",[{"max_visits":2},{"max_frontier":2},{"max_visits":0}])
def test_budgets(kwargs):
    with pytest.raises(ValueError,match="resource limit"):
        chematic.from_smiles("CC(C)C").rdkit_atom_environment(1,1,**kwargs)

def test_radius_zero_needs_no_frontier_budget():
    assert chematic.from_smiles("C").rdkit_atom_environment(0,0,max_visits=0,max_frontier=0)=={"bond_indices":[],"atom_distances":{0:0}}

def test_independent_results_and_source_properties():
    m=chematic.from_smiles("C1CCC1");before=m.smiles;m.SetIntProp("sentinel",9708)
    r=m.rdkit_atom_environment(2,0);r["bond_indices"].clear();r["atom_distances"][99]=88
    assert m.rdkit_atom_environment(2,0)["bond_indices"]==[0,3,1,2]
    assert 99 not in m.rdkit_atom_environment(0,0)["atom_distances"]
    assert m.smiles==before and m.GetIntProp("sentinel")==9708

@pytest.mark.parametrize("root",[1,99])
def test_indices_refer_to_h_suppressed_model(root):
    with pytest.raises(ValueError):chematic.from_smiles("[H]C([H])([H])[H]").rdkit_atom_environment(0,root)

def test_input_cap_empty_model_and_stubs():
    with pytest.raises(ValueError,match="1026 atoms"):
        chematic.from_smiles(".".join(["[H]C"]*513)).rdkit_atom_environment(0,0)
    with pytest.raises(ValueError):chematic.RWMol().GetMol().rdkit_atom_environment(0,0)
    p=Path(__file__).resolve().parents[1]/"python/chematic/__init__.pyi"
    cls=next(n for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n,ast.ClassDef) and n.name=="Mol")
    assert "rdkit_atom_environment" in {n.name for n in cls.body if isinstance(n,ast.FunctionDef)}
