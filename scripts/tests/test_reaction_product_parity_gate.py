"""The reaction evidence gate must not turn invalid products into parity wins."""

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem

from scripts.reaction_product_parity_gate import canonical_product_set, rdkit_products


def test_product_tuple_keeps_multiplicity() -> None:
    water = Chem.MolFromSmiles("O")
    assert canonical_product_set([(water, water)], "two_waters") == [["O", "O"]]


def test_invalid_rdkit_product_is_not_comparable() -> None:
    # Neutral pyridine N cannot acquire methyl while remaining aromatic;
    # RDKit can still emit a product Mol from this template before sanitizing.
    with pytest.raises(ValueError, match="unsanitizable product"):
        rdkit_products("[#7:1]>>[#7:1]C", ["c1ccncc1"], "invalid_pyridinium")
