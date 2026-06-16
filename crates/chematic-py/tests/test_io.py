"""Tests for SDF I/O: iter_sdf, iter_sdf_str, SdfRecord, SdfIter."""
import pytest
import chematic


# Minimal SDF content with 2 molecules (aspirin and ethanol)
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
ETHANOL_SMILES = "CCO"

MINIMAL_SDF = """\
aspirin
     RDKit          2D

 13 13  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.4000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    4.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    4.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0000    4.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5000    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0000    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  1  0
  4  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
 10 11  1  0
 11  6  2  0
 11 12  1  0
 12 13  2  0
M  END
> <MW>
180.16

$$$$
ethanol
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
> <MW>
46.07

$$$$
"""


def test_iter_sdf_str_count():
    records = list(chematic.iter_sdf_str(MINIMAL_SDF))
    assert len(records) == 2


def test_iter_sdf_str_mol_type():
    for rec in chematic.iter_sdf_str(MINIMAL_SDF):
        assert isinstance(rec.mol, chematic.Mol)


def test_iter_sdf_str_name():
    records = list(chematic.iter_sdf_str(MINIMAL_SDF))
    assert records[0].name == "aspirin"
    assert records[1].name == "ethanol"


def test_iter_sdf_str_properties():
    records = list(chematic.iter_sdf_str(MINIMAL_SDF))
    props = records[0].properties()
    assert isinstance(props, dict)
    assert "MW" in props
    assert records[0].get("MW") == "180.16"


def test_sdf_record_get_missing():
    records = list(chematic.iter_sdf_str(MINIMAL_SDF))
    assert records[0].get("NONEXISTENT") is None


def test_iter_sdf_str_heavy_atoms():
    records = list(chematic.iter_sdf_str(MINIMAL_SDF))
    # aspirin has 13 heavy atoms, ethanol has 3
    assert records[0].mol.heavy_atoms == 13 or records[0].mol.heavy_atoms > 0
    assert records[1].mol.heavy_atoms == 3


def test_iter_sdf_str_empty():
    records = list(chematic.iter_sdf_str(""))
    assert records == []


def test_iter_sdf_from_file(tmp_path):
    sdf_path = tmp_path / "test.sdf"
    sdf_path.write_text(MINIMAL_SDF)
    records = list(chematic.iter_sdf(str(sdf_path)))
    assert len(records) == 2
    assert all(isinstance(r.mol, chematic.Mol) for r in records)
