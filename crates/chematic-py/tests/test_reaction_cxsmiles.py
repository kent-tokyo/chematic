"""RDKit #9732: enhanced groups follow the emitted atom, not old offsets."""
import json
from pathlib import Path
import re

import pytest
import chematic

PACKET = json.loads((Path(__file__).parents[3] / "validation/rdkit-9732-2026-10-24/cases.json").read_text())


def membership(text):
    base, _, extension = text.partition("|")
    atoms = []
    for role, side in enumerate(re.split(r"(?<!-)>", base.strip())):
        for component in filter(None, side.split(".")):
            mol = chematic.from_smiles(component)
            cip = {entry["atom_idx"]: entry["descriptor"] for entry in mol.cip_stereo()}
            atoms.extend((role, mol.smiles, cip.get(i)) for i in range(mol.atom_count))
    result = []
    for group in re.finditer(r"(a|&[0-9]+|o[0-9]+):([0-9]+(?:,[0-9]+)*)", extension):
        for index in map(int, group[2].split(",")):
            assert index < len(atoms)
            role, component, cip = atoms[index]
            assert cip is not None, (text, index)
            result.append((role, component, group[1], cip))
    return sorted(result)


@pytest.mark.parametrize("case", PACKET["positive"], ids=lambda c: c["id"])
@pytest.mark.parametrize("canonical", [False, True])
def test_group_membership_and_role_survive_normalization(case, canonical):
    output = chematic.normalize_reaction_cxsmiles(case["input"], canonical=canonical)
    assert membership(output) == membership(case["input"])
    assert chematic.normalize_reaction_cxsmiles(output, canonical=canonical) == output


@pytest.mark.parametrize("case", PACKET["negative"], ids=lambda c: c["id"])
def test_unsupported_or_ambiguous_input_is_refused(case):
    with pytest.raises(ValueError):
        chematic.normalize_reaction_cxsmiles(case["input"])


def test_original_example_and_bytes_have_same_correct_member():
    source = PACKET["positive"][0]["input"]
    output = chematic.normalize_reaction_cxsmiles(source)
    assert output.endswith(" |&1:1|")
    assert output == chematic.normalize_reaction_cxsmiles(source.encode())
    with pytest.raises(ValueError):
        chematic.normalize_reaction_cxsmiles(b"\xff")
