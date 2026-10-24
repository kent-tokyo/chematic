"""MCS query output retains the alternatives omitted by concrete Mol projection."""

import pytest
import chematic as ch


@pytest.mark.parametrize("left,right,options", [
    ("NO", "CO", {"atom_compare": "any"}),
    ("C=N", "C#N", {"bond_compare": "any"}),
    ("C1CC1", "N1CC1", {"atom_compare": "any"}),
])
def test_mcs_smarts_keeps_alternatives(left, right, options):
    smarts, timed_out = ch.find_mcs_smarts_checked(
        [ch.from_smiles(left), ch.from_smiles(right)], **options
    )
    assert timed_out is False
    assert smarts is not None and "," in smarts


def test_mcs_smarts_reports_resource_limit():
    with pytest.raises(ch.ChematicInputError) as error:
        ch.find_mcs_smarts_checked(
            [ch.from_smiles("CCC"), ch.from_smiles("CCC")], max_search_nodes=0
        )
    assert error.value.code == "mcs_resource_limit"
