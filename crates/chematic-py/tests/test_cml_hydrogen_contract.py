"""Open Babel #2655 and total-vs-explicit CML hydrogen accounting."""
from pathlib import Path
import xml.etree.ElementTree as ET

import chematic
import pytest


CH2 = (Path(__file__).resolve().parents[3] /
       "validation/openbabel_issues/cml/explicit-ch2.cml").read_text(encoding="utf-8")


def with_count(value):
    return CH2.replace('id="carbon" elementType="C"',
                       f'id="carbon" elementType="C" hydrogenCount="{value}"')


@pytest.mark.parametrize("reader", [chematic.from_cml, chematic.from_cml_strict])
@pytest.mark.parametrize("count,formula", [(None, "CH2"), (2, "CH2"), (3, "CH3"), (4, "CH4")])
def test_ch2_and_total_count_preserve_composition(reader, count, formula):
    mol = reader(CH2 if count is None else with_count(count))
    assert mol.formula == formula
    assert reader(mol.to_cml()).formula == formula


@pytest.mark.parametrize("smiles", [
    "[H]C", "[H][C]", "[CH2]", "[13CH2][2H]", "[NH3+][H]",
    "[H]OC", "[H]c1ccccc1", "[H][H]",
])
def test_writer_preserves_mixed_hydrogens_isotopes_and_charge(smiles):
    mol = chematic.from_smiles(smiles)
    output = mol.to_cml()
    atoms = ET.fromstring(output).findall(".//{*}atom")
    assert atoms and all("hydrogenCount" in atom.attrib for atom in atoms)
    again = chematic.from_cml_strict(output)
    assert again.formula == mol.formula
    assert again.exact_mass == pytest.approx(mol.exact_mass, abs=1e-12)
    assert again.smiles == mol.smiles
    assert again.to_cml() == output


@pytest.mark.parametrize("reader", [chematic.from_cml, chematic.from_cml_strict])
@pytest.mark.parametrize("count", ["-1", "NaN", "3.5", "1", "258", "184467440737095516160"])
def test_bad_or_contradictory_count_is_an_explicit_error(reader, count):
    with pytest.raises(ValueError, match="hydrogenCount"):
        reader(with_count(count))


def test_methane_with_one_explicit_h_emits_total_four():
    xml = chematic.from_smiles("[H]C").to_cml()
    carbon = next(a for a in ET.fromstring(xml).findall(".//{*}atom")
                  if a.attrib["elementType"] == "C")
    assert carbon.attrib["hydrogenCount"] == "4"
