"""Basic smoke tests for the chematic Python bindings."""

import chematic
import pytest


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def test_from_smiles_benzene():
    mol = chematic.from_smiles("c1ccccc1")
    assert mol.formula == "C6H6"
    assert abs(mol.mw - 78.11) < 0.05


def test_from_smiles_ethanol():
    mol = chematic.from_smiles("CCO")
    assert mol.formula == "C2H6O"


def test_from_smiles_invalid():
    with pytest.raises(ValueError):
        chematic.from_smiles("[C")  # unclosed bracket


def test_is_valid_smiles():
    assert chematic.is_valid_smiles("CCO") is True
    assert chematic.is_valid_smiles("[C") is False  # unclosed bracket


def test_repr():
    mol = chematic.from_smiles("CCO")
    # canonical SMILES of ethanol is OCC
    assert "OCC" in repr(mol)


# ---------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------

def test_aspirin_descriptors():
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    assert abs(mol.mw - 180.16) < 0.05
    assert mol.hbd == 1
    assert mol.hba == 3  # Lipinski O+N count (3 oxygens that qualify)
    assert mol.rotatable_bonds == 3
    assert mol.heavy_atoms == 13
    assert mol.ring_count == 1
    assert mol.aromatic_ring_count == 1


def test_lipinski():
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert mol.lipinski_passes is True


def test_qed_range():
    mol = chematic.from_smiles("c1ccccc1")
    assert 0.0 <= mol.qed <= 1.0


def test_descriptors_dict():
    mol = chematic.from_smiles("CCO")
    d = mol.descriptors()
    assert isinstance(d, dict)
    assert "mw" in d
    assert "logp" in d
    assert "tpsa" in d
    assert "qed" in d
    assert "lipinski_passes" in d
    assert "bbb_score" in d


# ---------------------------------------------------------------------------
# pKa / ADMET
# ---------------------------------------------------------------------------

def test_pka_aspirin():
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    pka = mol.pka()
    assert "most_acidic" in pka
    assert "most_basic" in pka
    # aspirin is acidic
    assert pka["most_acidic"] is not None


def test_admet_aspirin():
    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    admet = mol.admet()
    assert "bbb" in admet
    assert "caco2" in admet
    assert "herg_risk" in admet
    assert "cyp3a4_risk" in admet
    assert isinstance(admet["bbb"], bool)
    assert isinstance(admet["caco2"], float)


# ---------------------------------------------------------------------------
# Fingerprints & similarity
# ---------------------------------------------------------------------------

def test_ecfp4_length():
    mol = chematic.from_smiles("c1ccccc1")
    fp = mol.ecfp4()
    assert isinstance(fp, bytes)
    assert len(fp) == 256


def test_maccs_length():
    mol = chematic.from_smiles("c1ccccc1")
    fp = mol.maccs()
    assert isinstance(fp, bytes)
    assert len(fp) == 21


def test_tanimoto_self():
    mol = chematic.from_smiles("c1ccccc1")
    fp = mol.ecfp4()
    assert abs(chematic.tanimoto(fp, fp) - 1.0) < 1e-9


def test_tanimoto_different():
    mol1 = chematic.from_smiles("c1ccccc1")      # benzene
    mol2 = chematic.from_smiles("c1cccnc1")      # pyridine
    sim = chematic.tanimoto(mol1.ecfp4(), mol2.ecfp4())
    assert 0.0 <= sim < 1.0


def test_tanimoto_length_mismatch():
    with pytest.raises(ValueError):
        chematic.tanimoto(b"\x00" * 10, b"\x00" * 20)


# ---------------------------------------------------------------------------
# SMARTS
# ---------------------------------------------------------------------------

def test_smarts_match_hydroxyl():
    mol = chematic.from_smiles("CCO")
    assert chematic.smarts_match("[OH]", mol) is True
    assert chematic.smarts_match("[NH2]", mol) is False


def test_smarts_match_aromatic():
    mol = chematic.from_smiles("c1ccccc1")
    assert chematic.smarts_match("c1ccccc1", mol) is True


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def test_svg_returns_string():
    mol = chematic.from_smiles("c1ccccc1")
    svg = mol.svg()
    assert isinstance(svg, str)
    assert "<svg" in svg


# ---------------------------------------------------------------------------
# Transformations
# ---------------------------------------------------------------------------

def test_scaffold():
    mol = chematic.from_smiles("CC1CCN(CC1)c1ncnc2[nH]ccc12")
    scaf = mol.scaffold()
    assert isinstance(scaf, chematic.Mol)


def test_inchi():
    mol = chematic.from_smiles("c1ccccc1")
    assert mol.inchi.startswith("InChI=")
    assert len(mol.inchikey) == 27


def test_version():
    assert hasattr(chematic, "__version__")
    assert chematic.__version__ == "0.3.2"


# ---------------------------------------------------------------------------
# numpy fingerprints (Mol methods)
# ---------------------------------------------------------------------------

def test_ecfp4_numpy_shape():
    import numpy as np
    mol = chematic.from_smiles("c1ccccc1")
    fp = mol.ecfp4_numpy()
    assert fp.shape == (2048,)
    assert fp.dtype == np.uint8
    assert set(fp.tolist()).issubset({0, 1})


def test_maccs_numpy_shape():
    import numpy as np
    mol = chematic.from_smiles("c1ccccc1")
    fp = mol.maccs_numpy()
    assert fp.shape == (166,)
    assert fp.dtype == np.uint8


# ---------------------------------------------------------------------------
# chematic.bulk — batch API
# ---------------------------------------------------------------------------

def test_bulk_parse():
    mols = chematic.bulk.parse(["c1ccccc1", "CCO", "[C"])
    assert len(mols) == 3
    assert isinstance(mols[0], chematic.Mol)
    assert isinstance(mols[1], chematic.Mol)
    assert mols[2] is None  # invalid SMILES


def test_bulk_ecfp4_shape():
    import numpy as np
    smiles = ["c1ccccc1", "CCO", "CC(=O)O"]
    fps = chematic.bulk.ecfp4(smiles)
    assert fps.shape == (3, 2048)
    assert fps.dtype == np.uint8
    assert set(fps.flatten().tolist()).issubset({0, 1})


def test_bulk_ecfp4_skips_invalid():
    import numpy as np
    smiles = ["c1ccccc1", "[C", "CCO"]  # middle one invalid
    fps = chematic.bulk.ecfp4(smiles)
    assert fps.shape == (2, 2048)  # 2 valid molecules


def test_bulk_maccs_shape():
    import numpy as np
    smiles = ["c1ccccc1", "CCO"]
    fps = chematic.bulk.maccs(smiles)
    assert fps.shape == (2, 166)
    assert fps.dtype == np.uint8


def test_bulk_descriptors():
    smiles = ["c1ccccc1", "CCO"]
    descs = chematic.bulk.descriptors(smiles)
    assert len(descs) == 2
    assert "mw" in descs[0]
    assert "logp" in descs[0]
    assert "qed" in descs[0]
    assert abs(descs[0]["mw"] - 78.11) < 0.1  # benzene


def test_bulk_tanimoto_shape():
    import numpy as np
    smiles_a = ["c1ccccc1", "CCO"]
    smiles_b = ["c1ccccc1", "CCO", "CC(=O)O"]
    sim = chematic.bulk.tanimoto(smiles_a, smiles_b)
    assert sim.shape == (2, 3)
    assert sim.dtype == np.float32
    # diagonal self-similarity
    assert abs(sim[0, 0] - 1.0) < 1e-5  # benzene vs benzene


def test_bulk_tanimoto_self_similarity():
    import numpy as np
    smiles = ["c1ccccc1", "CCO", "CC(=O)O"]
    sim = chematic.bulk.tanimoto(smiles, smiles)
    assert sim.shape == (3, 3)
    for i in range(3):
        assert abs(sim[i, i] - 1.0) < 1e-5, f"self-sim[{i}] = {sim[i, i]}"


def test_bulk_tanimoto_search():
    import numpy as np
    scores = chematic.bulk.tanimoto_search("c1ccccc1", ["c1ccccc1", "CCO", "c1cccnc1"])
    assert len(scores) == 3
    assert scores.dtype == np.float32
    assert abs(scores[0] - 1.0) < 1e-5  # self
    assert 0.0 <= scores[1] <= 1.0
    assert 0.0 <= scores[2] <= 1.0


# ---------------------------------------------------------------------------
# ESOL + standardize (Mol methods)
# ---------------------------------------------------------------------------

def test_esol():
    mol = chematic.from_smiles("CCO")  # ethanol
    sol = mol.esol
    assert isinstance(sol, float)
    # ethanol is quite soluble, ESOL should be > -3
    assert sol > -3.0


def test_standardize():
    mol = chematic.from_smiles("[NH4+].[Cl-]")  # ammonium chloride salt
    std = mol.standardize()
    assert isinstance(std, chematic.Mol)


# ---------------------------------------------------------------------------
# SDF streaming
# ---------------------------------------------------------------------------

SDF_TEXT = """aspirin
  chematic

 13 13  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -3.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971   -2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  8 10  1  0  0  0  0
 10 11  2  0  0  0  0
  2 12  1  0  0  0  0
 12 13  2  0  0  0  0
M  END
> <Activity>
5.3

> <Source>
test

$$$$
"""


def test_iter_sdf_str():
    records = list(chematic.iter_sdf_str(SDF_TEXT))
    assert len(records) >= 1
    rec = records[0]
    assert isinstance(rec, chematic.SdfRecord)
    assert isinstance(rec.mol, chematic.Mol)
    assert rec.name == "aspirin"
    props = rec.properties()
    assert "Activity" in props
    assert props["Activity"] == "5.3"
    assert rec.get("Activity") == "5.3"
    assert rec.get("NonExistent") is None


def test_iter_sdf_str_smiles():
    records = list(chematic.iter_sdf_str(SDF_TEXT))
    rec = records[0]
    # should be a valid SMILES string
    assert len(rec.smiles) > 0
    mol2 = chematic.from_smiles(rec.smiles)
    assert mol2.heavy_atoms > 0


# ---------------------------------------------------------------------------
# SimilarityIndex (LSH)
# ---------------------------------------------------------------------------

def test_similarity_index_basic():
    idx = chematic.SimilarityIndex()
    assert len(idx) == 0
    i = idx.add("c1ccccc1")
    assert i == 0
    j = idx.add("CCO")
    assert j == 1
    assert len(idx) == 2


def test_similarity_index_from_smiles():
    smiles = ["c1ccccc1", "CCO", "CC(=O)O", "c1cccnc1"]
    idx = chematic.SimilarityIndex.from_smiles(smiles)
    assert len(idx) == 4


def test_similarity_index_search():
    smiles_db = ["c1ccccc1", "CCO", "c1cccnc1", "c1ccoc1", "CC(=O)O"]
    idx = chematic.SimilarityIndex.from_smiles(smiles_db)
    # benzene should find itself (index 0, similarity ~1.0)
    hits = idx.search("c1ccccc1", threshold=0.5)
    assert len(hits) >= 1
    indices, scores = zip(*hits)
    assert 0 in indices  # benzene itself


def test_similarity_index_get_smiles():
    idx = chematic.SimilarityIndex()
    idx.add("c1ccccc1")
    idx.add("CCO")
    assert idx.get_smiles(0) == "c1ccccc1"
    assert idx.get_smiles(1) == "CCO"


def test_similarity_index_invalid_smiles():
    idx = chematic.SimilarityIndex()
    with pytest.raises(ValueError):
        idx.add("[C")  # invalid SMILES
