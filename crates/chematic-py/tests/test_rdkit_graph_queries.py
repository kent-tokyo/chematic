"""Original minimal graph contracts; independent oracle script pins RDKit."""
import ast
from pathlib import Path
import pytest
import chematic

@pytest.mark.parametrize("smiles,weighted", [
    ("N->[Cu]",[[0,0],[1,0]]), ("[Cu]<-N",[[0,1],[0,0]]),
    ("C=C",[[0,2],[2,0]]), ("C#N",[[0,3],[3,0]]),
    ("[Mo]$[Mo]",[[0,4],[4,0]]),
])
def test_order_and_endpoint_valence(smiles,weighted):
    m=chematic.from_smiles(smiles)
    assert m.rdkit_adjacency_matrix(use_bo=True)==weighted
    assert m.rdkit_adjacency_matrix()==[[0,1],[1,0]]
    assert m.rdkit_shortest_path(0,1)==[0,1]
    assert m.rdkit_shortest_path(1,0)==[1,0]

def test_tie_order_direction_and_disconnected_result():
    m=chematic.from_smiles("C1CCC1.O")
    assert m.rdkit_shortest_path(1,3)==[1,0,3]
    assert m.rdkit_shortest_path(3,1)==[3,2,1]
    assert m.rdkit_shortest_path(0,4)==[]
    assert m.rdkit_adjacency_matrix()[0][4]==0

def test_h_removal_and_added_h_indices():
    m=chematic.from_smiles("[H]C([H])([H])[H]")
    assert m.rdkit_adjacency_matrix()==[[0]]
    with pytest.raises(ValueError): m.rdkit_shortest_path(0,1)
    h=m.rdkit_add_hydrogens()
    a=h.rdkit_adjacency_matrix()
    assert len(a)==5
    assert a[0]==[0,1,1,1,1]
    assert h.rdkit_shortest_path(1,4)==[1,0,4]

@pytest.mark.parametrize("start,end",[(0,0),(0,2),(3,0)])
def test_invalid_endpoints(start,end):
    with pytest.raises(ValueError):
        chematic.from_smiles("CC").rdkit_shortest_path(start,end)

@pytest.mark.parametrize("value",[-1,1,127,255,256])
def test_nonzero_fill_is_explicitly_unsupported(value):
    with pytest.raises(ValueError, match="nonzero adjacency empty_value"):
        chematic.from_smiles("CC.O").rdkit_adjacency_matrix(empty_value=value)

def test_empty_graph_return_shapes():
    m=chematic.RWMol().GetMol()
    assert m.rdkit_adjacency_matrix()==[]
    with pytest.raises(ValueError): m.rdkit_adjacency_matrix(empty_value=1)
    with pytest.raises(ValueError): m.rdkit_shortest_path(0,1)

def test_returned_arrays_and_source_properties_are_independent():
    m=chematic.from_smiles("CCO")
    before=m.smiles
    m.SetIntProp("sentinel",9707)
    a=m.rdkit_adjacency_matrix();a[0][1]=500
    path=m.rdkit_shortest_path(0,2);path.clear()
    assert m.rdkit_adjacency_matrix()[0][1]==1
    assert m.rdkit_shortest_path(0,2)==[0,1,2]
    assert m.smiles==before and m.GetIntProp("sentinel")==9707

def test_preparation_limit_applies_before_h_removal():
    m=chematic.from_smiles(".".join(["[H]C"]*513))
    for call in [lambda:m.rdkit_adjacency_matrix(),lambda:m.rdkit_shortest_path(0,1)]:
        with pytest.raises(ValueError,match="1026 atoms"): call()

def test_stub_methods_are_declared_on_mol():
    p=Path(__file__).resolve().parents[1]/"python/chematic/__init__.pyi"
    cls=next(n for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n,ast.ClassDef) and n.name=="Mol")
    assert {"rdkit_adjacency_matrix","rdkit_shortest_path"}<={n.name for n in cls.body if isinstance(n,ast.FunctionDef)}
