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


@pytest.mark.parametrize(
    ("packed_method", "numpy_method", "n_bits"),
    [
        ("ecfp6", "ecfp6_numpy", 2048),
        ("ecfp4_chiral", "ecfp4_chiral_numpy", 2048),
        ("fcfp4", "fcfp4_numpy", 2048),
        ("fcfp6", "fcfp6_numpy", 2048),
        ("atom_pair_fp", "atom_pair_fp_numpy", 2048),
        ("torsion_fp", "torsion_fp_numpy", 2048),
        ("rdkit_torsion_fp", "rdkit_torsion_fp_numpy", 2048),
        ("rdkit_atom_pair_fp", "rdkit_atom_pair_fp_numpy", 2048),
        ("rdkit_layered_fp", "rdkit_layered_fp_numpy", 2048),
        ("rdkit_pattern_fp", "rdkit_pattern_fp_numpy", 2048),
        ("rdkit_rdk_fp", "rdkit_rdk_fp_numpy", 2048),
        ("rdkit_ecfp4", "rdkit_ecfp4_numpy", 2048),
        ("path_fp", "path_fp_numpy", 2048),
        ("topo_path_fp", "topo_path_fp_numpy", 2048),
        ("avalon_fp", "avalon_fp_numpy", 2048),
        ("rdkit_avalon_fp", "rdkit_avalon_fp_numpy", 2048),
        ("pattern_fp", "pattern_fp_numpy", 2048),
        ("pharmacophore_fp", "pharmacophore_fp_numpy", 2048),
        ("erg", "erg_numpy", 2048),
    ],
)
def test_binary_fingerprint_numpy_matches_packed_bits(
    aspirin, packed_method, numpy_method, n_bits
):
    import numpy as np

    packed = getattr(aspirin, packed_method)()
    expected = np.unpackbits(
        np.frombuffer(packed, dtype=np.uint8), bitorder="little"
    )[:n_bits]
    actual = getattr(aspirin, numpy_method)()

    assert actual.shape == (n_bits,)
    assert actual.dtype == np.uint8
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("count_simulation", [False, True])
def test_configurable_morgan_numpy_matches_packed_bits(aspirin, count_simulation):
    import numpy as np

    packed = aspirin.rdkit_ecfp_config(
        radius=2,
        nbits=512,
        include_chirality=True,
        count_simulation=count_simulation,
    )
    actual = aspirin.rdkit_ecfp_config_numpy(
        radius=2,
        nbits=512,
        include_chirality=True,
        count_simulation=count_simulation,
    )
    expected = np.unpackbits(
        np.frombuffer(packed, dtype=np.uint8), bitorder="little"
    )

    assert actual.shape == (512,)
    assert actual.dtype == np.uint8
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("n_bits", [0, 7, 2049])
def test_rdkit_avalon_numpy_rejects_non_byte_aligned_widths(aspirin, n_bits):
    with pytest.raises(ValueError, match="positive multiple of 8"):
        aspirin.rdkit_avalon_fp_numpy(n_bits)


def test_3d_pharmacophore_numpy_matches_packed_bits(aspirin):
    import numpy as np

    coords = [[float(i), 0.0, 0.0] for i in range(aspirin.rdkit_num_atoms())]
    packed = aspirin.pharmacophore_fp_3d(coords)
    actual = aspirin.pharmacophore_fp_3d_numpy(coords)
    expected = np.unpackbits(
        np.frombuffer(packed, dtype=np.uint8), bitorder="little"
    )

    assert actual.shape == (2048,)
    assert actual.dtype == np.uint8
    np.testing.assert_array_equal(actual, expected)


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


def test_rdkit_atom_pair_fp_length(aspirin):
    assert len(aspirin.rdkit_atom_pair_fp()) == 256


def test_rdkit_atom_pair_fp_different(aspirin, benzene):
    assert aspirin.rdkit_atom_pair_fp() != benzene.rdkit_atom_pair_fp()


def test_rdkit_atom_pair_fp_deterministic(aspirin):
    assert aspirin.rdkit_atom_pair_fp() == aspirin.rdkit_atom_pair_fp()


def test_rdkit_atom_pair_fp_independent_of_native_atom_pair_fp(aspirin):
    # Separate opt-in function -- must not be the same bytes as the native
    # (non-RDKit) scheme just because both happen to be 256-byte outputs.
    assert aspirin.rdkit_atom_pair_fp() != aspirin.atom_pair_fp()


def test_rdkit_atom_code_and_torsion_path_explanations_match_oracle():
    mol = chematic.from_smiles("C=CC(=O)O")
    expected = [
        ("C", 1, 1),
        ("C", 2, 1),
        ("C", 3, 1),
        ("O", 1, 1),
        ("O", 1, 0),
    ]
    for atom_index, fields in enumerate(expected):
        code = mol.rdkit_atom_code(atom_index)
        assert chematic.explain_rdkit_atom_code(code) == fields

    path = [0, 1, 2, 4]
    score = mol.rdkit_topological_torsion_path_code(path)
    assert chematic.explain_rdkit_topological_torsion_code(score) == [
        ("C", 1, 1),
        ("C", 2, 1),
        ("C", 3, 1),
        ("O", 1, 0),
    ]
    assert score == mol.rdkit_topological_torsion_path_code(path[::-1])

    try:
        from rdkit import Chem
        from rdkit.Chem.AtomPairs import Torsions, Utils
    except ImportError:
        return
    rd_mol = Chem.MolFromSmiles("C=CC(=O)O")
    assert score == Torsions.pyScorePath(rd_mol, path, 4)
    assert tuple(chematic.explain_rdkit_topological_torsion_code(score)) == (
        Torsions.ExplainPathScore(score, 4)
    )
    for atom_index in range(rd_mol.GetNumAtoms()):
        assert mol.rdkit_atom_code(atom_index) == Utils.GetAtomCode(
            rd_mol.GetAtomWithIdx(atom_index)
        )


def test_rdkit_atom_code_diagnostics_reject_invalid_inputs():
    mol = chematic.from_smiles("CC")
    with pytest.raises(ValueError, match="out of range"):
        mol.rdkit_atom_code(2)
    with pytest.raises(ValueError, match="must not be empty"):
        mol.rdkit_topological_torsion_path_code([])
    with pytest.raises(ValueError, match="at most 7"):
        chematic.explain_rdkit_topological_torsion_code(0, 8)


def test_rdkit_atom_pair_fp_matches_rdkit_oracle_on_simple_molecules():
    # Regression pin: these were the molecules used to verify the atom-pair
    # port (shares its atom-invariant scheme with rdkit_torsion_fp, so these
    # overlap that fingerprint's own repro set) -- must stay bit-exact.
    import chematic

    cases = [
        "CCO",
        "CCCC",
        "C1CC1",
        "CC1CC1",
        "COC1CC1",
    ]
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
    except ImportError:
        return  # rdkit not installed in this environment -- skip
    for smi in cases:
        m_chem = chematic.from_smiles(smi)
        m_rd = Chem.MolFromSmiles(smi)
        rd_bits = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            m_rd, nBits=2048
        ).ToBitString()
        chem_bytes = m_chem.rdkit_atom_pair_fp()
        chem_bits = "".join(
            format(byte, "08b")[::-1] for byte in chem_bytes
        )[:2048]
        assert chem_bits == rd_bits, f"{smi}: parity mismatch vs RDKit oracle"


def _bits(fp_bytes):
    return set(i for i in range(2048) if (fp_bytes[i // 8] >> (i % 8)) & 1)


def test_rdkit_pattern_fp_substructure_bits_are_a_subset():
    # Pattern fingerprints exist for substructure screening: a substructure's
    # own bits must be a subset of its parent molecule's bits, or screening
    # would produce false negatives. Checked on whole *concrete* molecules
    # (this port doesn't support fingerprinting a SMARTS query molecule
    # itself -- see the module's own doc comment) where one is a literal
    # substructure of the other.
    import chematic

    pairs = [
        ("c1ccccc1", "c1ccc2ccccc2c1"),
        ("CCO", "CCOC"),
        ("c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"),
    ]
    for sub, full in pairs:
        b_sub = _bits(chematic.from_smiles(sub).rdkit_pattern_fp())
        b_full = _bits(chematic.from_smiles(full).rdkit_pattern_fp())
        assert b_sub <= b_full, f"{sub} bits not a subset of {full} bits"


def test_rdkit_pattern_fp_length(aspirin):
    assert len(aspirin.rdkit_pattern_fp()) == 256


def test_rdkit_pattern_fp_different(aspirin, benzene):
    assert aspirin.rdkit_pattern_fp() != benzene.rdkit_pattern_fp()


def test_rdkit_pattern_fp_deterministic(aspirin):
    assert aspirin.rdkit_pattern_fp() == aspirin.rdkit_pattern_fp()


def test_rdkit_pattern_fp_independent_of_native_pattern_fp(aspirin):
    # Separate opt-in function -- must not be the same bytes as the native
    # (non-RDKit) scheme just because both happen to be 256-byte outputs.
    assert aspirin.rdkit_pattern_fp() != aspirin.pattern_fp()


def test_rdkit_pattern_fp_matches_rdkit_oracle_on_simple_molecules():
    # Regression pin, including the Kekule-notation heteroaromatic that
    # exposed the aromaticity-perception bug this port originally had (a raw
    # bond.order == Aromatic check silently missed every ring bond of a
    # Kekule-written aromatic ring) -- must stay bit-exact.
    import chematic

    cases = [
        "CCO",
        "CCCC",
        "C1CC1",
        "c1ccccc1",
        "c1ccc2ccccc2c1",
        "S1C2=CC3=CC=CC=C3C=C2N=C1C4=CC=CC=C4",
    ]
    try:
        from rdkit import Chem
    except ImportError:
        return  # rdkit not installed in this environment -- skip
    for smi in cases:
        m_chem = chematic.from_smiles(smi)
        m_rd = Chem.MolFromSmiles(smi)
        rd_bits = Chem.PatternFingerprint(m_rd, fpSize=2048).ToBitString()
        chem_bytes = m_chem.rdkit_pattern_fp()
        chem_bits = "".join(
            format(byte, "08b")[::-1] for byte in chem_bytes
        )[:2048]
        assert chem_bits == rd_bits, f"{smi}: parity mismatch vs RDKit oracle"


def test_rdkit_rdk_fp_length(aspirin):
    assert len(aspirin.rdkit_rdk_fp()) == 256


def test_rdkit_rdk_fp_different(aspirin, benzene):
    assert aspirin.rdkit_rdk_fp() != benzene.rdkit_rdk_fp()


def test_rdkit_rdk_fp_deterministic(aspirin):
    assert aspirin.rdkit_rdk_fp() == aspirin.rdkit_rdk_fp()


def _on_bits(packed: bytes) -> list[int]:
    return [
        byte_index * 8 + bit
        for byte_index, byte in enumerate(packed)
        for bit in range(8)
        if byte & (1 << bit)
    ]


def test_rdkit_rdk_fp_ignore_atoms_matches_rdkit_2026_09_1():
    mol = chematic.from_smiles("CCO")
    assert _on_bits(mol.rdkit_rdk_fp(ignore_atoms=[2])) == [1308, 1772]
    assert _on_bits(mol.rdkit_rdk_fp(ignore_atoms=[1])) == []
    assert mol.rdkit_rdk_fp(ignore_atoms=[]) == mol.rdkit_rdk_fp()


def test_rdkit_rdk_fp_ignore_atoms_rejects_invalid_index():
    with pytest.raises(ValueError, match="out of range"):
        chematic.from_smiles("CCO").rdkit_rdk_fp(ignore_atoms=[3])


def test_rdkit_rdk_fp_independent_of_native_path_fp(aspirin):
    # Separate opt-in function -- must not be the same bytes as chematic's own
    # pre-existing, non-bit-exact linear-path approximation just because both
    # happen to be 256-byte outputs.
    assert aspirin.rdkit_rdk_fp() != aspirin.path_fp()


def test_rdkit_rdk_fp_matches_rdkit_oracle_on_simple_molecules():
    # Regression pin for the branched-subgraph enumeration (root-at-min-bond-index
    # backtracking), path-local (not molecular) atom degree in the bond hash, and
    # RDKit's own weakened Mersenne Twister used for the second bit per feature --
    # these specific molecules were used to validate each mechanism and must stay
    # bit-exact.
    import chematic

    cases = [
        "CC",
        "C=C",
        "C#C",
        "CCC",
        "CC(C)C",
        "c1ccccc1",
        "C1CC1",
        "CCO",
        "CC(=O)O",
        "c1ccc2ccccc2c1",
        "CC(=O)Oc1ccccc1C(=O)O",
    ]
    try:
        from rdkit import Chem
    except ImportError:
        return  # rdkit not installed in this environment -- skip
    for smi in cases:
        m_chem = chematic.from_smiles(smi)
        m_rd = Chem.MolFromSmiles(smi)
        rd_bits = Chem.RDKFingerprint(
            m_rd, minPath=1, maxPath=7, fpSize=2048
        ).ToBitString()
        chem_bytes = m_chem.rdkit_rdk_fp()
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


def test_rdkit_layered_fp_length(aspirin):
    assert len(aspirin.rdkit_layered_fp()) == 256


def test_rdkit_layered_fp_different(aspirin, benzene):
    assert aspirin.rdkit_layered_fp() != benzene.rdkit_layered_fp()


def test_rdkit_layered_fp_deterministic(aspirin):
    assert aspirin.rdkit_layered_fp() == aspirin.rdkit_layered_fp()


def test_rdkit_layered_fp_independent_of_native_layered_fp(aspirin):
    # Separate opt-in function -- must not be the same bytes as chematic's own
    # pre-existing, non-bit-exact scheme just because both happen to be
    # 256-byte outputs.
    assert aspirin.rdkit_layered_fp() != aspirin.layered_fp()


def test_rdkit_layered_fp_matches_rdkit_oracle_on_simple_molecules():
    # Regression pin, including the isotope-labeled-hydrogen case that exposed
    # a real bug this port originally had (LayeredFingerprintMol enumerates
    # with useHs=false, unlike rdkit_rdk_fp's useHs=true -- a naive useHs=true
    # port silently included bonds touching an explicit deuterium atom).
    import chematic

    cases = [
        "CC",
        "C=C",
        "C#C",
        "CCC",
        "CC(C)C",
        "c1ccccc1",
        "C1CC1",
        "CCO",
        "CC(=O)O",
        "c1ccc2ccccc2c1",
        "CC(=O)Oc1ccccc1C(=O)O",
        "[2H]C([2H])([2H])NC=O",
    ]
    try:
        from rdkit import Chem
    except ImportError:
        return  # rdkit not installed in this environment -- skip
    for smi in cases:
        m_chem = chematic.from_smiles(smi)
        m_rd = Chem.MolFromSmiles(smi)
        rd_bits = Chem.LayeredFingerprint(
            m_rd, minPath=1, maxPath=7, fpSize=2048
        ).ToBitString()
        chem_bytes = m_chem.rdkit_layered_fp()
        chem_bits = "".join(
            format(byte, "08b")[::-1] for byte in chem_bytes
        )[:2048]
        assert chem_bits == rd_bits, f"{smi}: parity mismatch vs RDKit oracle"


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


def test_rdkit_avalon_fp_matches_rdkit(aspirin):
    # RDKit 2026.03.1 GetAvalonFP(aspirin, nBits=512) on-bits
    raw = bytes(aspirin.rdkit_avalon_fp(512))
    assert len(raw) == 64
    bits = [i * 8 + j for i, b in enumerate(raw) for j in range(8) if b >> j & 1]
    assert bits[:5] == [13, 15, 41, 67, 85] and len(bits) == 54
