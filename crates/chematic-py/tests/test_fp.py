"""Tests for fingerprint methods and similarity functions."""
import pytest
import chematic


@pytest.fixture(scope="module")
def aspirin():
    return chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")


@pytest.fixture(scope="module")
def benzene():
    return chematic.from_smiles("c1ccccc1")


@pytest.fixture(scope="module")
def toluene():
    return chematic.from_smiles("Cc1ccccc1")


# ---------------------------------------------------------------------------
# ECFP4 / ECFP6
# ---------------------------------------------------------------------------

def test_ecfp4_length(aspirin):
    fp = aspirin.ecfp4()
    assert isinstance(fp, bytes)
    assert len(fp) == 256


def test_ecfp6_length(aspirin):
    fp = aspirin.ecfp6()
    assert len(fp) == 256


def test_ecfp4_different_mols(aspirin, benzene):
    assert aspirin.ecfp4() != benzene.ecfp4()


def test_ecfp4_same_mol_deterministic(aspirin):
    assert aspirin.ecfp4() == aspirin.ecfp4()


def test_ecfp4_numpy(aspirin):
    import numpy as np
    fp = aspirin.ecfp4_numpy()
    assert fp.shape == (2048,)
    assert fp.dtype == np.uint8


# ---------------------------------------------------------------------------
# FCFP4
# ---------------------------------------------------------------------------

def test_fcfp4_length(aspirin):
    fp = aspirin.fcfp4()
    assert len(fp) == 256


def test_fcfp4_differs_from_ecfp4(aspirin):
    # FCFP4 and ECFP4 should generally differ
    assert aspirin.fcfp4() != aspirin.ecfp4()


# ---------------------------------------------------------------------------
# Atom-pair / Torsion
# ---------------------------------------------------------------------------

def test_atom_pair_fp_length(aspirin):
    assert len(aspirin.atom_pair_fp()) == 256


def test_torsion_fp_length(aspirin):
    assert len(aspirin.torsion_fp()) == 256


def test_atom_pair_fp_different(aspirin, benzene):
    assert aspirin.atom_pair_fp() != benzene.atom_pair_fp()


def test_rdkit_torsion_fp_length(aspirin):
    assert len(aspirin.rdkit_torsion_fp()) == 256


def test_rdkit_torsion_fp_different(aspirin, benzene):
    assert aspirin.rdkit_torsion_fp() != benzene.rdkit_torsion_fp()


def test_rdkit_torsion_fp_deterministic(aspirin):
    assert aspirin.rdkit_torsion_fp() == aspirin.rdkit_torsion_fp()


def test_rdkit_torsion_fp_independent_of_native_torsion_fp(aspirin):
    # Separate opt-in function -- must not be the same bytes as the native
    # (non-RDKit) scheme just because both happen to be 256-byte outputs.
    assert aspirin.rdkit_torsion_fp() != aspirin.torsion_fp()


def test_rdkit_torsion_fp_matches_rdkit_oracle_on_simple_molecules():
    # Regression pin for the 3 real bugs found and fixed this round (missing
    # -2 topological torsion correction, double-counted torsion paths,
    # missing 3-membered-ring closure entries) -- these specific molecules
    # were used to isolate each bug and must stay bit-exact.
    import chematic

    cases = [
        "CCCC",  # linear butane, no ring
        "C1CC1",  # bare cyclopropane -- only a ring-closure torsion exists
        "CC1CC1",  # methylcyclopropane -- symmetric ring substitution
        "COC1CC1",  # cyclopropyl ether
    ]
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return  # rdkit not installed in this environment -- skip
    for smi in cases:
        m_chem = chematic.from_smiles(smi)
        m_rd = Chem.MolFromSmiles(smi)
        rd_bits = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(
            m_rd, nBits=2048
        ).ToBitString()
        chem_bytes = m_chem.rdkit_torsion_fp()
        chem_bits = "".join(
            format(byte, "08b")[::-1] for byte in chem_bytes
        )[:2048]
        assert chem_bits == rd_bits, f"{smi}: parity mismatch vs RDKit oracle"


# ---------------------------------------------------------------------------
# MACCS
# ---------------------------------------------------------------------------

def test_maccs_length(aspirin):
    # MACCS 166-bit → 21 bytes (ceil(166/8))
    assert len(aspirin.maccs()) == 21


def test_maccs_numpy_shape(aspirin):
    import numpy as np
    fp = aspirin.maccs_numpy()
    assert fp.shape == (166,)
    assert fp.dtype == np.uint8


def test_maccs_different_mols(aspirin, benzene):
    assert aspirin.maccs() != benzene.maccs()


# ---------------------------------------------------------------------------
# ECFP4 chiral
# ---------------------------------------------------------------------------

def test_ecfp4_chiral_length(aspirin):
    assert len(aspirin.ecfp4_chiral()) == 256


# ---------------------------------------------------------------------------
# Tanimoto similarity
# ---------------------------------------------------------------------------

def test_tanimoto_identical():
    fp = chematic.from_smiles("c1ccccc1").ecfp4()
    assert chematic.tanimoto(fp, fp) == pytest.approx(1.0)


def test_tanimoto_different(aspirin, benzene):
    sim = chematic.tanimoto(aspirin.ecfp4(), benzene.ecfp4())
    assert 0.0 <= sim < 1.0


def test_tanimoto_similar_molecules(benzene, toluene):
    sim = chematic.tanimoto(benzene.ecfp4(), toluene.ecfp4())
    # RDKit's Morgan(r=2, 2048 bits) gives 0.27 for this pair too — a folded
    # 2048-bit ECFP4 has real bit collisions, so 0.3 overstated the analogy.
    assert sim > 0.2  # toluene is a close analog of benzene


def test_tanimoto_length_mismatch():
    with pytest.raises(ValueError):
        chematic.tanimoto(b"\x00" * 256, b"\x00" * 21)


# ---------------------------------------------------------------------------
# B5: Layered fingerprint (SMARTS/structural layers)
# ---------------------------------------------------------------------------

def test_layered_fp_vs_ecfp4(aspirin):
    # Layered FP should differ from ECFP4 (different encoding)
    assert aspirin.layered_fp() != aspirin.ecfp4()


def test_layered_fp_tanimoto(benzene, toluene):
    sim = chematic.tanimoto(benzene.layered_fp(), toluene.layered_fp())
    assert 0.0 <= sim <= 1.0


# ---------------------------------------------------------------------------
# Avalon fingerprint
# ---------------------------------------------------------------------------

def test_avalon_fp_length(aspirin):
    fp = aspirin.avalon_fp()
    assert isinstance(fp, bytes)
    assert len(fp) == 256


def test_avalon_fp_deterministic(aspirin):
    assert aspirin.avalon_fp() == aspirin.avalon_fp()


def test_avalon_fp_different_mols(aspirin, benzene):
    assert aspirin.avalon_fp() != benzene.avalon_fp()


def test_avalon_fp_tanimoto_similarity_ordering(benzene, toluene, aspirin):
    sim_toluene = chematic.tanimoto(benzene.avalon_fp(), toluene.avalon_fp())
    sim_aspirin = chematic.tanimoto(benzene.avalon_fp(), aspirin.avalon_fp())
    assert sim_toluene > sim_aspirin
