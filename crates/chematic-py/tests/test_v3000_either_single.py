"""Open Babel #2429: an explicit either carrier must not acquire E/Z."""
from pathlib import Path

import pytest
import chematic

SOURCE = (Path(__file__).resolve().parents[3] / "validation/openbabel-2429-2026-10-24/either-single.mol").read_text()
COORDS = [[-1, 1], [0, 0], [1, 0], [2, -1]]


def variant(side, reverse, mirror, three_d):
    text = SOURCE.replace(" CFG=2", "")
    a, b = ("10", "20") if side == "left" else ("30", "40")
    atom1, atom2 = (b, a) if reverse else (a, b)
    bond_id = "50" if side == "left" else "70"
    text = text.replace(f"{bond_id} 1 {a} {b}", f"{bond_id} 1 {atom1} {atom2} CFG=2")
    if mirror:
        text = text.replace("40 C 2 -1", "40 C 2 1")
    if three_d:
        text = text.replace(" 2D", " 3D").replace("10 C -1 1 0", "10 C -1 1 0.5")
    return text


@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("mirror", [False, True])
@pytest.mark.parametrize("three_d", [False, True])
def test_either_carriers_keep_alkene_stereo_unknown(side, reverse, mirror, three_d):
    source = variant(side, reverse, mirror, three_d)
    mol = chematic.from_mol_v3000(source)
    assert mol.atom_count == 4
    assert mol.formula == "C4H8"
    assert "/" not in mol.smiles and "\\" not in mol.smiles
    before = mol.to_moljson()
    output = mol.to_mol_v3000(COORDS, strict=True)
    reread = chematic.from_mol_v3000(output)
    assert reread.smiles == mol.smiles
    assert mol.to_moljson() == before
    # Exact source retention is a separate pre-existing API; Mol alone does
    # not own the V3000 source metadata containing the wavy bond's location.
    assert chematic.roundtrip_mol_v3000_source(source) == source


@pytest.mark.parametrize("cfg", ["", " CFG=2"])
def test_unmarked_and_direct_double_controls(cfg):
    source = SOURCE.replace(" CFG=2", "").replace("60 2 20 30", f"60 2 20 30{cfg}")
    smiles = chematic.from_mol_v3000(source).smiles
    assert (("/" in smiles) or ("\\" in smiles)) == (cfg == "")
