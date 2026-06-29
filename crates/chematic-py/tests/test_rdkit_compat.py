"""Tests for chematic.rdkit_compat — RDKit API compatibility layer."""
import pytest
import chematic
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import (
    Descriptors, rdMolDescriptors, DataStructs, ExplicitBitVect,
    Atom, Bond, BondType, RingInfo,
)


# ---------------------------------------------------------------------------
# Import style
# ---------------------------------------------------------------------------

def test_import_style():
    assert hasattr(Chem, "MolFromSmiles")
    assert hasattr(Chem, "SDMolSupplier")
    assert hasattr(Chem, "SDWriter")
    assert Descriptors is not None
    assert rdMolDescriptors is not None
    assert DataStructs is not None


# ---------------------------------------------------------------------------
# Mol property methods
# ---------------------------------------------------------------------------

def test_setprop_getprop():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("ID", "ethanol")
    assert mol.GetProp("ID") == "ethanol"


def test_hasprop():
    mol = Chem.MolFromSmiles("CCO")
    assert not mol.HasProp("x")
    mol.SetProp("x", "1")
    assert mol.HasProp("x")


def test_getpropsasdict():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("A", "1")
    mol.SetProp("B", "2")
    d = mol.GetPropsAsDict()
    assert d["A"] == "1" and d["B"] == "2"


def test_getpropnames():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("foo", "bar")
    assert "foo" in mol.GetPropNames()


def test_clearprop():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("x", "1")
    mol.ClearProp("x")
    assert not mol.HasProp("x")


def test_clearprop_missing_noop():
    mol = Chem.MolFromSmiles("CCO")
    mol.ClearProp("nonexistent")  # must not raise


def test_setintprop():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetIntProp("rank", 42)
    assert mol.GetProp("rank") == "42"


def test_setdoubleprop():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetDoubleProp("score", 3.14)
    assert mol.GetProp("score") == "3.14"


def test_setboolprop():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetBoolProp("active", True)
    assert mol.GetProp("active") == "1"
    mol.SetBoolProp("active", False)
    assert mol.GetProp("active") == "0"


def test_getprop_missing_raises():
    mol = Chem.MolFromSmiles("CCO")
    with pytest.raises(KeyError):
        mol.GetProp("missing")


# ---------------------------------------------------------------------------
# SDWriter + SDMolSupplier roundtrip
# ---------------------------------------------------------------------------

def test_sdf_roundtrip_props(tmp_path):
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol.SetProp("ID", "benzene")
    mol.SetProp("MW", "78.11")

    out = tmp_path / "out.sdf"
    with Chem.SDWriter(str(out)) as w:
        w.write(mol)

    mols = list(Chem.SDMolSupplier(str(out)))
    assert len(mols) == 1
    m = mols[0]
    assert m is not None
    assert m.GetProp("ID") == "benzene"
    assert m.GetProp("MW") == "78.11"


def test_sdf_name_in_header_not_sd_field(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "ethanol")
    mol.SetProp("Activity", "7.2")

    out = tmp_path / "out.sdf"
    with Chem.SDWriter(str(out)) as w:
        w.write(mol)

    raw = out.read_text()
    assert "ethanol" in raw              # name in MOL header
    assert "> <_Name>" not in raw        # must not appear as SD field
    assert "> <Activity>" in raw


def test_setprops_filters_output(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("ID", "ethanol")
    mol.SetProp("MW", "46.07")
    mol.SetProp("Source", "internal")

    out = tmp_path / "out.sdf"
    w = Chem.SDWriter(str(out))
    w.SetProps(["ID", "MW"])
    w.write(mol)
    w.close()

    raw = out.read_text()
    assert "> <ID>" in raw
    assert "> <MW>" in raw
    assert "> <Source>" not in raw


def test_setprops_empty_writes_no_sd_fields(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("ID", "ethanol")

    out = tmp_path / "out.sdf"
    w = Chem.SDWriter(str(out))
    w.SetProps([])
    w.write(mol)
    w.close()

    raw = out.read_text()
    assert "> <ID>" not in raw
    assert "$$$$" in raw


def test_typed_props_roundtrip(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    mol.SetIntProp("rank", 7)
    mol.SetDoubleProp("score", 1.5)
    mol.SetBoolProp("active", True)

    out = tmp_path / "out.sdf"
    with Chem.SDWriter(str(out)) as w:
        w.write(mol)

    mols = list(Chem.SDMolSupplier(str(out)))
    m = mols[0]
    assert m.GetProp("rank") == "7"
    assert m.GetProp("score") == "1.5"
    assert m.GetProp("active") == "1"


# ---------------------------------------------------------------------------
# flush / close idempotency
# ---------------------------------------------------------------------------

def test_flush_and_close_idempotent(tmp_path):
    out = tmp_path / "out.sdf"
    w = Chem.SDWriter(str(out))
    w.write(Chem.MolFromSmiles("CCO"))
    w.flush()
    w.flush()   # second flush must not crash
    w.close()
    w.close()   # second close must not crash


# ---------------------------------------------------------------------------
# SetKekulize / SetForceV3000 are no-ops (must not crash)
# ---------------------------------------------------------------------------

def test_noop_writer_flags(tmp_path):
    out = tmp_path / "out.sdf"
    w = Chem.SDWriter(str(out))
    w.SetKekulize(False)
    w.SetForceV3000(True)
    w.write(Chem.MolFromSmiles("c1ccccc1"))
    w.close()
    assert out.exists()


# ---------------------------------------------------------------------------
# Context manager
# ---------------------------------------------------------------------------

def test_context_manager(tmp_path):
    out = tmp_path / "out.sdf"
    with Chem.SDWriter(str(out)) as w:
        w.write(Chem.MolFromSmiles("c1ccccc1"))
    assert out.exists()
    mols = list(Chem.SDMolSupplier(str(out)))
    assert len(mols) == 1


# ---------------------------------------------------------------------------
# SDMolSupplier __len__
# ---------------------------------------------------------------------------

def test_sdmolsupplier_len(tmp_path):
    out = tmp_path / "out.sdf"
    with Chem.SDWriter(str(out)) as w:
        for smi in ["CCO", "c1ccccc1", "CC(=O)O"]:
            w.write(Chem.MolFromSmiles(smi))
    assert len(Chem.SDMolSupplier(str(out))) == 3


# ---------------------------------------------------------------------------
# Descriptors / rdMolDescriptors smoke tests
# ---------------------------------------------------------------------------

def test_descriptors_mw():
    mol = Chem.MolFromSmiles("CCO")
    assert abs(Descriptors.MolWt(mol) - 46.07) < 0.1


def test_rdmoldescriptors_tpsa():
    mol = Chem.MolFromSmiles("CCO")
    assert rdMolDescriptors.CalcTPSA(mol) >= 0


def test_tanimoto_self():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    assert DataStructs.TanimotoSimilarity(fp, fp) == 1.0


# ---------------------------------------------------------------------------
# ExplicitBitVect
# ---------------------------------------------------------------------------

def test_explicit_bit_vect_basic():
    bv = ExplicitBitVect(16)
    assert bv.GetNumBits() == 16
    assert bv.GetBit(0) is False

    bv.SetBit(0)
    bv.SetBit(7)
    bv.SetBit(15)
    assert bv.GetBit(0) is True
    assert bv.GetBit(7) is True
    assert bv.GetBit(15) is True
    assert bv.GetBit(1) is False

    assert bv.GetOnBits() == [0, 7, 15]
    bs = bv.ToBitString()
    assert len(bs) == 16
    assert bs[0] == "1" and bs[7] == "1" and bs[15] == "1" and bs[1] == "0"


def test_explicit_bit_vect_index_error():
    bv = ExplicitBitVect(8)
    with pytest.raises(IndexError):
        bv.GetBit(8)
    with pytest.raises(IndexError):
        bv.SetBit(-1)


# ---------------------------------------------------------------------------
# GetMorganFingerprintAsBitVect returns ExplicitBitVect
# ---------------------------------------------------------------------------

def test_morgan_returns_explicit_bit_vect():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    assert isinstance(fp, ExplicitBitVect)
    assert fp.GetNumBits() == 2048
    bs = fp.ToBitString()
    assert len(bs) == 2048
    assert "1" in bs  # benzene must set some bits


def test_morgan_radius3():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 3)
    assert isinstance(fp, ExplicitBitVect)


# ---------------------------------------------------------------------------
# Tanimoto with ExplicitBitVect
# ---------------------------------------------------------------------------

def test_tanimoto_identical_bitvect():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    assert DataStructs.TanimotoSimilarity(fp, fp) == 1.0


def test_tanimoto_disjoint():
    bv1 = ExplicitBitVect(8)
    bv2 = ExplicitBitVect(8)
    bv1.SetBit(0)
    bv2.SetBit(7)
    assert DataStructs.TanimotoSimilarity(bv1, bv2) == 0.0


def test_bulk_tanimoto():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    results = DataStructs.BulkTanimotoSimilarity(fp, [fp, fp])
    assert results == [1.0, 1.0]


def test_bulk_tanimoto_different_mols():
    benzene = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles("c1ccccc1"), 2
    )
    ethanol = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles("CCO"), 2
    )
    results = DataStructs.BulkTanimotoSimilarity(benzene, [benzene, ethanol])
    assert results[0] == 1.0
    assert 0.0 <= results[1] < 1.0


# ---------------------------------------------------------------------------
# Unsupported options fail loudly
# ---------------------------------------------------------------------------

def test_morgan_use_features_raises():
    mol = Chem.MolFromSmiles("c1ccccc1")
    with pytest.raises(NotImplementedError):
        rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True)


def test_morgan_unknown_kwarg_raises():
    mol = Chem.MolFromSmiles("c1ccccc1")
    with pytest.raises(TypeError):
        rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, unknownParam=True)


# ---------------------------------------------------------------------------
# Mol / Atom / Bond read-only surface
# ---------------------------------------------------------------------------

def test_atom_bond_counts():
    mol = Chem.MolFromSmiles("CCO")  # ethanol: C, C, O = 3 heavy atoms, 2 bonds
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2


def test_atom_symbols():
    mol = Chem.MolFromSmiles("CCO")
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    assert "C" in symbols
    assert "O" in symbols


def test_atomic_num():
    mol = Chem.MolFromSmiles("CCO")
    nums = {a.GetSymbol(): a.GetAtomicNum() for a in mol.GetAtoms()}
    assert nums["C"] == 6
    assert nums["O"] == 8


def test_formal_charge_neutral():
    mol = Chem.MolFromSmiles("CCO")
    assert all(a.GetFormalCharge() == 0 for a in mol.GetAtoms())


def test_aromatic_atoms_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    assert all(a.GetIsAromatic() for a in mol.GetAtoms())


def test_aromatic_atoms_ethanol():
    mol = Chem.MolFromSmiles("CCO")
    assert not any(a.GetIsAromatic() for a in mol.GetAtoms())


def test_aromatic_bonds_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    assert all(b.GetIsAromatic() for b in mol.GetBonds())
    assert all(b.GetBondType() == BondType.AROMATIC for b in mol.GetBonds())


def test_bond_type_single():
    mol = Chem.MolFromSmiles("CC")
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == BondType.SINGLE
    assert bond.GetBondTypeAsDouble() == 1.0


def test_bond_type_double():
    mol = Chem.MolFromSmiles("C=C")
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == BondType.DOUBLE
    assert bond.GetBondTypeAsDouble() == 2.0


def test_bond_type_aromatic_double():
    mol = Chem.MolFromSmiles("c1ccccc1")
    assert mol.GetBondWithIdx(0).GetBondTypeAsDouble() == 1.5


def test_degree():
    mol = Chem.MolFromSmiles("CO")  # methanol: C has degree 1 (one O neighbor)
    o_atom = next(a for a in mol.GetAtoms() if a.GetSymbol() == "O")
    assert o_atom.GetDegree() == 1


def test_total_num_hs():
    mol = Chem.MolFromSmiles("CCO")  # CH3-CH2-OH: first C has 3 implicit Hs
    atoms = list(mol.GetAtoms())
    # The methyl carbon (degree 1 in ethanol) has 3 implicit Hs
    methyl_c = next(a for a in atoms if a.GetSymbol() == "C" and a.GetDegree() == 1)
    assert methyl_c.GetTotalNumHs() == 3


def test_is_in_ring():
    benzene = Chem.MolFromSmiles("c1ccccc1")
    assert all(a.IsInRing() for a in benzene.GetAtoms())
    ethanol = Chem.MolFromSmiles("CCO")
    assert not any(a.IsInRing() for a in ethanol.GetAtoms())


def test_get_atom_with_idx():
    mol = Chem.MolFromSmiles("CCO")
    atom = mol.GetAtomWithIdx(0)
    assert isinstance(atom, Atom)
    assert atom.GetIdx() == 0


def test_get_bond_with_idx():
    mol = Chem.MolFromSmiles("CCO")
    bond = mol.GetBondWithIdx(0)
    assert isinstance(bond, Bond)
    assert bond.GetIdx() == 0


def test_atom_index_error():
    mol = Chem.MolFromSmiles("CCO")
    with pytest.raises(IndexError):
        mol.GetAtomWithIdx(-1)
    with pytest.raises(IndexError):
        mol.GetAtomWithIdx(99)


def test_bond_index_error():
    mol = Chem.MolFromSmiles("CCO")
    with pytest.raises(IndexError):
        mol.GetBondWithIdx(999)


def test_get_begin_end_atom():
    mol = Chem.MolFromSmiles("CO")
    bond = mol.GetBondWithIdx(0)
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    assert isinstance(a1, Atom)
    assert isinstance(a2, Atom)
    assert {a1.GetIdx(), a2.GetIdx()} == {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}


def test_get_other_atom_idx():
    mol = Chem.MolFromSmiles("CO")
    bond = mol.GetBondWithIdx(0)
    a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    assert bond.GetOtherAtomIdx(a1) == a2
    assert bond.GetOtherAtomIdx(a2) == a1


def test_get_other_atom_idx_invalid():
    mol = Chem.MolFromSmiles("CO")
    bond = mol.GetBondWithIdx(0)
    with pytest.raises(ValueError):
        bond.GetOtherAtomIdx(99)


def test_pyridine_n_atomic_num():
    mol = Chem.MolFromSmiles("c1ccncc1")
    atomic_nums = [a.GetAtomicNum() for a in mol.GetAtoms()]
    assert 7 in atomic_nums  # nitrogen


# ---------------------------------------------------------------------------
# RingInfo
# ---------------------------------------------------------------------------

def test_ring_info_acyclic():
    mol = Chem.MolFromSmiles("CCO")
    ri = mol.GetRingInfo()
    assert isinstance(ri, RingInfo)
    assert ri.NumRings() == 0
    assert ri.AtomRings() == ()
    assert ri.BondRings() == ()


def test_ring_info_cyclohexane():
    mol = Chem.MolFromSmiles("C1CCCCC1")
    ri = mol.GetRingInfo()
    assert ri.NumRings() == 1
    assert len(ri.AtomRings()) == 1
    assert len(ri.AtomRings()[0]) == 6
    assert len(ri.BondRings()[0]) == 6


def test_ring_info_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    ri = mol.GetRingInfo()
    assert ri.NumRings() == 1
    assert len(ri.BondRings()[0]) == 6


def test_ring_info_naphthalene():
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    ri = mol.GetRingInfo()
    assert ri.NumRings() == 2
    # The shared bond must appear in 2 bond rings
    flat_bonds = [b for ring in ri.BondRings() for b in ring]
    # At least one bond appears twice (shared bond)
    assert max(flat_bonds.count(b) for b in flat_bonds) == 2


def test_atom_rings_type():
    mol = Chem.MolFromSmiles("c1ccccc1")
    ar = mol.GetRingInfo().AtomRings()
    assert isinstance(ar, tuple)
    assert isinstance(ar[0], tuple)


def test_bond_rings_type():
    mol = Chem.MolFromSmiles("c1ccccc1")
    br = mol.GetRingInfo().BondRings()
    assert isinstance(br, tuple)
    assert isinstance(br[0], tuple)


def test_num_atom_rings_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    ri = mol.GetRingInfo()
    assert all(ri.NumAtomRings(i) == 1 for i in range(mol.GetNumAtoms()))


def test_num_bond_rings_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    ri = mol.GetRingInfo()
    assert all(ri.NumBondRings(i) == 1 for i in range(mol.GetNumBonds()))


def test_ring_info_invalid_idx():
    mol = Chem.MolFromSmiles("c1ccccc1")
    ri = mol.GetRingInfo()
    with pytest.raises(IndexError):
        ri.NumAtomRings(-1)
    with pytest.raises(IndexError):
        ri.NumBondRings(999)


def test_bond_is_in_ring_cyclohexane():
    mol = Chem.MolFromSmiles("C1CCCCC1")
    assert all(b.IsInRing() for b in mol.GetBonds())


def test_bond_not_in_ring_ethanol():
    mol = Chem.MolFromSmiles("CCO")
    assert not any(b.IsInRing() for b in mol.GetBonds())


def test_atom_is_in_ring_size_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    atom = mol.GetAtomWithIdx(0)
    assert atom.IsInRingSize(6) is True
    assert atom.IsInRingSize(5) is False


# ---------------------------------------------------------------------------
# Morgan FP nBits folding
# ---------------------------------------------------------------------------

def test_morgan_nbits_1024():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    assert fp.GetNumBits() == 1024
    assert len(fp.ToBitString()) == 1024


def test_morgan_nbits_4096():
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096)
    assert fp.GetNumBits() == 4096


def test_morgan_nbits_tanimoto():
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    assert DataStructs.TanimotoSimilarity(fp, fp) == 1.0


def test_morgan_nbits_zero_raises():
    mol = Chem.MolFromSmiles("c1ccccc1")
    with pytest.raises(ValueError):
        rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=0)


# ---------------------------------------------------------------------------
# ConvertToNumpyArray
# ---------------------------------------------------------------------------

def test_convert_to_numpy():
    import numpy as np
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    arr = DataStructs.ConvertToNumpyArray(fp)
    assert arr.shape == (2048,)
    on = fp.GetOnBits()
    assert all(arr[i] == 1 for i in on)
    assert int(arr.sum()) == len(on)


def test_convert_to_numpy_dest():
    import numpy as np
    mol = Chem.MolFromSmiles("c1ccccc1")
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    dest = np.zeros(2048, dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, dest)
    assert int(dest.sum()) == len(fp.GetOnBits())


# ---------------------------------------------------------------------------
# GetSubstructMatch / GetSubstructMatches RDKit shape
# ---------------------------------------------------------------------------

def test_get_substruct_match():
    mol = Chem.MolFromSmiles("c1ccccc1")
    m = mol.GetSubstructMatch("c1ccccc1")
    assert isinstance(m, tuple)
    assert len(m) == 6


def test_get_substruct_match_none():
    mol = Chem.MolFromSmiles("CCO")
    assert mol.GetSubstructMatch("c1ccccc1") == ()


def test_get_substruct_matches_tuple():
    mol = Chem.MolFromSmiles("c1ccccc1O")
    matches = mol.GetSubstructMatches("c")
    assert isinstance(matches, tuple)
    assert all(isinstance(m, tuple) for m in matches)
    assert len(matches) == 6  # six aromatic carbons


def test_get_substruct_matches_uniquify():
    mol = Chem.MolFromSmiles("c1ccccc1")
    # Without uniquify a symmetric pattern yields more matches than with it
    uniq = mol.GetSubstructMatches("c1ccccc1", uniquify=True)
    full = mol.GetSubstructMatches("c1ccccc1", uniquify=False)
    assert len(uniq) <= len(full)
    assert len(uniq) == 1  # single ring, one atom set


# ---------------------------------------------------------------------------
# Morgan bitInfo
# ---------------------------------------------------------------------------

def test_morgan_bitinfo_populated():
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    bitInfo = {}
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=bitInfo)
    assert len(bitInfo) > 0
    assert len(bitInfo) == len(fp.GetOnBits())


def test_morgan_bitinfo_shape():
    mol = Chem.MolFromSmiles("c1ccccc1")
    bitInfo = {}
    rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=bitInfo)
    for bit, envs in bitInfo.items():
        assert isinstance(envs, tuple)
        for env in envs:
            assert isinstance(env, tuple)
            assert len(env) == 2  # (atom_idx, radius)


def test_morgan_bitinfo_keys_are_onbits():
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    bitInfo = {}
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=bitInfo)
    for bit in bitInfo:
        assert fp.GetBit(bit) is True


def test_morgan_bitinfo_radius_range():
    mol = Chem.MolFromSmiles("c1ccncc1")
    bitInfo = {}
    rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=bitInfo)
    for envs in bitInfo.values():
        for atom_idx, radius in envs:
            assert 0 <= radius <= 2


def test_morgan_bitinfo_atom_range():
    mol = Chem.MolFromSmiles("c1ccncc1")
    n = mol.GetNumAtoms()
    bitInfo = {}
    rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=bitInfo)
    for envs in bitInfo.values():
        for atom_idx, radius in envs:
            assert 0 <= atom_idx < n


def test_morgan_bitinfo_none_noop():
    mol = Chem.MolFromSmiles("c1ccccc1")
    # bitInfo=None must work without a dict (default path)
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    assert fp.GetNumBits() == 2048


def test_morgan_bitinfo_folded():
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    bitInfo = {}
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, 2, nBits=1024, bitInfo=bitInfo
    )
    assert fp.GetNumBits() == 1024
    assert all(bit < 1024 for bit in bitInfo)
    for bit in bitInfo:
        assert fp.GetBit(bit) is True


def test_morgan_bitinfo_chirality_raises():
    mol = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
    with pytest.raises(NotImplementedError):
        rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, 2, useChirality=True, bitInfo={}
        )
