"""Open Babel #2333: format-local XYZ species normalization."""
from pathlib import Path

import pytest
import chematic


@pytest.mark.parametrize("symbol,canonical", [("CL","Cl"),("BR","Br"),("CA","Ca"),("FE","Fe"),("NA","Na"),("OG","Og")])
def test_aliases_in_plain_and_extended_bindings(symbol, canonical):
    coords = [[1., -2.5, 1.25]]
    for reader in ["plain", "extended"]:
        text = f"1\ncase\n{symbol} 1 -2.5 1.25\n"
        if reader == "plain":
            mol, actual = chematic.from_xyz(text)
        else:
            record = chematic.from_extxyz(text)
            mol, actual = record["mol"], record["coords"]
        assert actual == coords
        assert mol.to_xyz(actual).splitlines()[2].split()[0] == canonical
        assert mol.atom_count == 1


@pytest.mark.parametrize("token", ["cl","cL","ZZ","CL1","CLO","17","D","T","*","ＣＬ"])
def test_invalid_tokens_remain_errors(token):
    text = f"1\ninvalid\n{token} 0 0 0\n"
    for reader in [chematic.from_xyz, chematic.from_extxyz]:
        with pytest.raises(ValueError, match="unknown element"):
            reader(text)


def test_only_species_changes_not_extra_string_fields():
    record = chematic.from_extxyz('1\nProperties=species:S:1:pos:R:3:name:S:1 label=CL\nCL 1 2 3 CA\n')
    assert record["properties"]["name"] == [["CA"]]
    assert record["info"]["label"] == "CL"
    assert record["coords"] == [[1.,2.,3.]]
    assert record["mol"].to_xyz(record["coords"]).splitlines()[2].split()[0] == "Cl"


def test_original_explicit_h_fixture_keeps_all_atoms():
    path = Path(__file__).resolve().parents[3] / "validation/openbabel_issues/xyz-elements/chloromethane-uppercase.xyz"
    mol, coords = chematic.from_xyz(path.read_text(encoding="utf-8"))
    assert mol.atom_count == 5
    assert len(coords) == 5
    assert [r.split()[0] for r in mol.to_xyz(coords).splitlines()[2:]] == ["C","Cl","H","H","H"]


def test_trajectory_binding_reads_aliases_in_every_frame():
    records = chematic.from_extxyz_all("1\nfirst\nCL 0 0 0\n1\nsecond\nBR 1 2 3\n")
    assert [r["mol"].to_xyz(r["coords"]).splitlines()[2].split()[0] for r in records] == ["Cl","Br"]


def test_smiles_keeps_its_own_case_rules():
    with pytest.raises(ValueError):
        chematic.from_smiles("CL")
    assert chematic.from_smiles("CCl").atom_count == 2
