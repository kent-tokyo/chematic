import pytest

import chematic


def _mols(*smiles: str) -> list[chematic.Mol]:
    return [chematic.from_smiles(value) for value in smiles]


def test_labeled_rgroup_rows_align_symmetric_columns_like_rdkit() -> None:
    rows = chematic.rgroup_decompose_labeled(
        "c1cc([*:1])ccc1[*:2]",
        _mols("Cc1ccc(CC)cc1", "CCc1ccc(N)cc1"),
    )

    assert rows[0] is not None
    assert rows[1] is not None
    assert "CC" in rows[0]["R1"]
    assert "C" in rows[0]["R2"]
    assert "CC" in rows[1]["R1"]
    assert "N" in rows[1]["R2"]
    assert "*:1" in rows[0]["R1"]
    assert "*:2" in rows[0]["R2"]


def test_labeled_rgroup_columns_preserve_unmatched_rows() -> None:
    columns = chematic.rgroup_decompose_columns(
        "c1cc([*:1])ccc1[*:2]",
        _mols("Cc1ccc(CC)cc1", "CCCC"),
    )

    assert set(columns) == {"Core", "R1", "R2"}
    assert all(len(values) == 2 for values in columns.values())
    assert columns["Core"][1] is None
    assert columns["R1"][1] is None
    assert columns["R2"][1] is None


def test_labeled_rgroup_refuses_ambiguous_label_contracts() -> None:
    with pytest.raises(ValueError, match="appears more than once"):
        chematic.rgroup_decompose_labeled(
            "[c:1]1cccc[c:1]1", _mols("Cc1ccccc1")
        )

    with pytest.raises(ValueError, match="must be terminal"):
        chematic.rgroup_decompose_labeled(
            "c1cc([*:2]C)ccc1", _mols("Cc1ccccc1")
        )
