"""Original graph/title reuse controls inspired by Open Babel #3032."""
import ast
from pathlib import Path

import chematic
from chematic import rdkit_compat as compat
import pytest


SMILES = ['CCO', 'c1ccccc1', 'F[C@](Cl)(Br)I', '[13CH3:7][NH3+]']


def editable(source, binding):
    return chematic.RWMol(source) if binding == 'native' else compat.RWMol(compat.Mol(source))


def native_snapshot(rw, binding):
    snapshot = rw.GetMol()
    return snapshot if binding == 'native' else snapshot._mol


@pytest.mark.parametrize('binding', ['native', 'compat'])
@pytest.mark.parametrize('smiles', SMILES)
def test_complete_reset_preserves_source_and_snapshots(smiles, binding):
    source = chematic.from_smiles(smiles)
    if '@' in smiles:
        source = source.with_stereo_groups([('and', 17, [1], [])])
    source.SetName('old title 日本語')
    source.SetProp('text', 'retain on source')
    source.SetIntProp('integer', 17)
    source.SetDoubleProp('double', 0.75)
    source.SetBoolProp('flag', True)
    before = (source.smiles, source.formula, source.mw, source.tpsa,
              source.stereo_groups, source.GetPropsAsDict())
    rw = editable(source, binding)
    snapshot = native_snapshot(rw, binding)
    assert rw.Clear() is None
    empty = native_snapshot(rw, binding)
    assert rw.GetNumAtoms() == rw.GetNumBonds() == 0
    assert rw.GetName() == empty.GetName() == ''
    assert empty.GetPropNames() == [] and not empty.HasProp('_Name')
    assert empty.stereo_groups == []
    fresh = chematic.RWMol().GetMol()
    assert empty.formula == fresh.formula
    assert empty.mw == fresh.mw == 0 and empty.tpsa == fresh.tpsa == 0
    for unchanged in [source, snapshot]:
        assert (unchanged.smiles, unchanged.formula, unchanged.mw, unchanged.tpsa,
                unchanged.stereo_groups, unchanged.GetPropsAsDict()) == before
    # Reuse the same editor, including index zero. Old isotope/map/stereo and
    # properties must not survive into a new two-atom graph with no title.
    assert rw.AddAtom(6) == 0
    assert rw.AddAtom(8) == 1
    assert rw.AddBond(0, 1, 'SINGLE') == 1
    rebuilt = native_snapshot(rw, binding)
    fresh_rebuilt = chematic.from_smiles('CO')
    assert rebuilt.smiles == fresh_rebuilt.smiles
    assert rebuilt.to_mol_block() == fresh_rebuilt.to_mol_block()
    assert rebuilt.GetName() == '' and rebuilt.GetPropNames() == []
    assert rebuilt.stereo_groups == []
    assert empty.atom_count == 0
    assert snapshot.GetName() == 'old title 日本語'


@pytest.mark.parametrize('binding', ['native', 'compat'])
def test_empty_clear_is_idempotent_and_repeated_reuse_has_no_title_leak(binding):
    rw = editable(chematic.from_smiles('C'), binding)
    for i in range(100):
        rw.SetName(f'cycle {i}')
        assert rw.GetName() == f'cycle {i}'
        rw.Clear()
        rw.Clear()
        assert rw.GetNumAtoms() == rw.GetNumBonds() == 0
        assert rw.GetName() == ''
        assert native_snapshot(rw, binding).GetPropNames() == []
        assert rw.AddAtom(7 if i % 2 else 8) == 0
        assert rw.GetName() == ''


def test_native_clear_stub_and_compat_annotation():
    root = Path(__file__).resolve().parents[1] / 'python/chematic'
    tree = ast.parse((root / '__init__.pyi').read_text(encoding='utf-8'))
    owner = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'RWMol')
    method = next(n for n in owner.body if isinstance(n, ast.FunctionDef) and n.name == 'Clear')
    assert [arg.arg for arg in method.args.args] == ['self']
    assert ast.unparse(method.returns) == 'None'
    assert compat.RWMol.Clear.__annotations__['return'] in (None, 'None')
