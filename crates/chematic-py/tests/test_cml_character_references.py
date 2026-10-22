"""Original CML lexical/typed-field packet; refusal is separate from coverage."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] /
                     "validation/cml_character_references/cases.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("reader", [chematic.from_cml, chematic.from_cml_strict])
@pytest.mark.parametrize("case", PACKET["cases"], ids=lambda c: c["id"])
def test_cml_lexical_variants_and_explicit_rejections(reader, case):
    if case["expected_refusal"]:
        with pytest.raises(ValueError):
            reader(case["source"])
        return
    mol = reader(case["source"])
    root = ET.fromstring(mol.to_cml())
    atoms = root.findall(".//{*}atom")
    assert len(atoms) == 2
    for atom, expected in zip(atoms, PACKET["expected"]["atoms"]):
        a = atom.attrib
        assert a["elementType"] == expected["element"]
        assert int(a.get("formalCharge", 0)) == expected["charge"]
        assert int(a.get("isotopeNumber", 0)) == expected["isotope"]
        assert int(a["hydrogenCount"]) == expected["hydrogens"]
    # Python graph-only CML APIs intentionally do not retain coordinates.
