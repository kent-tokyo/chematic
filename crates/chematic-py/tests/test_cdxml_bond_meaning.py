"""CDXML graph class/direction plus explicit refusal for unsupported orders."""
import json
from pathlib import Path

import chematic
import pytest

CASES = json.loads((Path(__file__).parents[3] / "validation/cdxml_bond_meaning/cases.json").read_text())["cases"]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_public_cdxml_bond_contract(case):
    if case["group"] == "refusal":
        for infer in (False, True):
            with pytest.raises(ValueError, match="unsupported CDXML bond order"):
                chematic.from_cdxml(case["input"], infer_nondirectional_stereo=infer)
        for target in ("smiles", "cdxml"):
            with pytest.raises(ValueError, match="unsupported CDXML bond order"):
                chematic.convert_format(case["input"], "cdxml", target)
        return
    oracle = case["expected"]["smiles"]
    expected = chematic.from_smiles(oracle).smiles
    for infer in (False, True):
        assert chematic.from_cdxml(case["input"], infer_nondirectional_stereo=infer).smiles == expected
    assert chematic.convert_format(case["input"], "cdxml", "smiles") == expected
    output = chematic.convert_format(case["input"], "cdxml", "cdxml")
    assert chematic.from_cdxml(output).smiles == expected
    # Isolate output from the input CDXML reader by starting from authored SMILES.
    output = chematic.convert_format(oracle, "smiles", "cdxml")
    assert chematic.from_cdxml(output).smiles == expected
