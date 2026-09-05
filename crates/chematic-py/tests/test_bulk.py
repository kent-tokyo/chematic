"""Tests for bulk parallel processing (chematic.bulk module)."""
import pytest
import chematic


SMILES_LIST = [
    "c1ccccc1",            # benzene
    "CCO",                  # ethanol
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine
]


def test_bulk_parse():
    mols = chematic.bulk.parse(SMILES_LIST)
    assert len(mols) == 4
    for m in mols:
        assert m is not None
        assert isinstance(m, chematic.Mol)


def test_bulk_parse_with_invalid():
    smiles_with_bad = SMILES_LIST + ["NOT_VALID???"]
    mols = chematic.bulk.parse(smiles_with_bad)
    assert len(mols) == 5
    # Invalid SMILES returns None
    assert mols[-1] is None


def test_bulk_ecfp4_matrix():
    import numpy as np
    matrix = chematic.bulk.ecfp4(SMILES_LIST)
    assert matrix.shape == (4, 2048)
    assert matrix.dtype == np.uint8


def test_bulk_descriptors_rows():
    rows = chematic.bulk.descriptors(SMILES_LIST)
    assert len(rows) == 4
    for row in rows:
        assert isinstance(row, dict)
        assert "mw" in row
        assert "logp" in row
        assert "tpsa" in row


def test_bulk_descriptors_array_selected_columns():
    import numpy as np

    result = chematic.bulk.descriptors_array(SMILES_LIST, ["mw", "tpsa", "hbd"])
    assert set(result) == {"mw", "tpsa", "hbd"}
    assert result["mw"].shape == (4,)
    assert result["mw"].dtype == np.float64
    assert result["hbd"].shape == (4,)

    with pytest.raises(ValueError, match="unknown column"):
        chematic.bulk.descriptors_array(SMILES_LIST, ["not_a_descriptor"])


def test_bulk_tanimoto_matrix():
    import numpy as np
    sim_matrix = chematic.bulk.tanimoto_matrix(SMILES_LIST)
    n = len(SMILES_LIST)
    assert sim_matrix.shape == (n, n)
    assert sim_matrix.dtype in (np.float32, np.float64)
    # Diagonal should be 1.0 (self-similarity)
    for i in range(n):
        assert sim_matrix[i][i] == pytest.approx(1.0, abs=0.01)


def test_bulk_tanimoto_matrix_symmetry():
    import numpy as np
    sim = chematic.bulk.tanimoto_matrix(SMILES_LIST[:3])
    assert sim.shape == (3, 3)
    # Symmetry: sim[i][j] == sim[j][i]
    for i in range(3):
        for j in range(3):
            assert abs(sim[i][j] - sim[j][i]) < 1e-6
