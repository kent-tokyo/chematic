import pytest

import chematic


def test_activity_cliffs_rejects_mismatched_input_lengths():
    mols = [chematic.from_smiles("CC"), chematic.from_smiles("CCC")]
    with pytest.raises(ValueError, match="must have the same length"):
        chematic.activity_cliffs(mols, [1.0])
