"""RDKit #9359: expose existing conservative CDXML inference to Python."""
from pathlib import Path
import re

import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
TEXT = (ROOT / "validation/rdkit_issues/cdxml/nondirectional-stereo.cdxml").read_text(encoding="utf-8")
AMBIGUOUS = (ROOT / "validation/rdkit_issues/cdxml/ambiguous-nondirectional.cdxml").read_text(encoding="utf-8")


def codes(mol):
    return [(v["atom_idx"], v["descriptor"]) for v in mol.cip_stereo()]


@pytest.mark.parametrize("display", ["Bold", "Hash", "Dash"])
def test_default_and_explicit_false_preserve_unassigned_stereo(display):
    text = TEXT.replace("Bold", display)
    default = chematic.from_cdxml(text)
    explicit = chematic.from_cdxml(cdxml_str=text, infer_nondirectional_stereo=False)
    assert codes(default) == codes(explicit) == []
    assert default.smiles == explicit.smiles
    assert default.atom_count == 6


@pytest.mark.parametrize("display,expected", [("Bold", "S"), ("Hash", "R"), ("Dash", "R")])
def test_opt_in_has_the_directional_reference_configuration(display, expected):
    mol = chematic.from_cdxml(TEXT.replace("Bold", display), infer_nondirectional_stereo=True)
    assert codes(mol) == [(0, expected)]
    assert mol.atom_count == 6


@pytest.mark.parametrize("display", ["Bold", "Hash", "Dash"])
def test_nondirectional_endpoint_reversal_does_not_invert_stereo(display):
    text = TEXT.replace("Bold", display)
    reversed_text = text.replace('B="1" E="2" Display=', 'B="2" E="1" Display=')
    a = chematic.from_cdxml(text, infer_nondirectional_stereo=True)
    b = chematic.from_cdxml(reversed_text, infer_nondirectional_stereo=True)
    assert codes(a) == codes(b)
    assert a.smiles == b.smiles


@pytest.mark.parametrize("display", ["WedgeBegin", ""])
def test_directional_and_plain_displays_do_not_depend_on_the_flag(display):
    text = TEXT.replace(' Display="Bold"', f' Display="{display}"' if display else "")
    a = chematic.from_cdxml(text)
    b = chematic.from_cdxml(text, infer_nondirectional_stereo=True)
    assert codes(a) == codes(b)
    assert a.smiles == b.smiles


@pytest.mark.parametrize("text", [AMBIGUOUS, re.sub(r' p="[^"]*"', "", TEXT)])
def test_ambiguous_or_missing_geometry_abstains(text):
    mol = chematic.from_cdxml(text, infer_nondirectional_stereo=True)
    assert codes(mol) == []


def test_first_fragment_contract_remains_unchanged():
    text = TEXT.replace("</page>", '<fragment><n id="99" Element="8"/></fragment></page>')
    assert chematic.from_cdxml(text, infer_nondirectional_stereo=True).atom_count == 6


def test_option_is_keyword_only():
    with pytest.raises(TypeError):
        chematic.from_cdxml(TEXT, True)


def test_invalid_document_still_raises_value_error():
    with pytest.raises(ValueError):
        chematic.from_cdxml('<CDXML><fragment><n id="1" Element="999"/></fragment></CDXML>', infer_nondirectional_stereo=True)
