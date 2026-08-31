"""Tests for the common format conversion bridge."""

import pytest

import chematic


def test_convert_graph_formats_round_trip():
    mol2 = chematic.convert_format("CCO", "smiles", "mol2")
    assert "@<TRIPOS>MOLECULE" in mol2
    assert chematic.are_identical(
        chematic.from_smiles(chematic.convert_format(mol2, ".mol2", "smi")),
        chematic.from_smiles("CCO"),
    )

    cjson = chematic.convert_format(mol2, "mol2", "cjson")
    assert chematic.convert_format(cjson, "cjson", "moljson")
    assert chematic.are_identical(
        chematic.from_smiles(chematic.convert_format(cjson, "cjson", "smiles")),
        chematic.from_smiles("CCO"),
    )


def test_convert_coordinate_formats_use_input_coordinates():
    xyz = """3
ethanol
C 0.0 0.0 0.0
C 1.5 0.0 0.0
O 2.8 0.0 0.0
"""
    pdb = chematic.convert_format(xyz, "xyz", "pdb")
    assert "ATOM" in pdb or "HETATM" in pdb
    pdbqt = chematic.convert_format(xyz, "xyz", "pdbqt", name="ETH")
    assert "REMARK  NAME = ETH" in pdbqt

    xyz2 = chematic.convert_format(pdb, "pdb", "xyz", comment="round-trip")
    assert "round-trip" in xyz2


def test_convert_requires_coordinates_for_3d_output():
    with pytest.raises(ValueError, match="PDB output requires coords"):
        chematic.convert_format("CCO", "smiles", "pdb")

    with pytest.raises(ValueError, match="charges must have the same length"):
        chematic.convert_format(
            "CCO", "smiles", "pdbqt", coords=[[0, 0, 0]] * 3, charges=[0.0]
        )


def test_convert_rejects_unknown_formats():
    with pytest.raises(ValueError, match="unsupported molecular format"):
        chematic.convert_format("CCO", "smiles", "foo")
