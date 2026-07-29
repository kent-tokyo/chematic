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


# ── parse_sdf_with_coords: issue #171 (blank MOL name line) ─────────────────

MOL_A = """mol_a
  chematic

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"""

# Live RDKit 2026.03.3 Chem.MolToMolBlock() output for an unnamed molecule
# (AllChem.Compute2DCoords(Chem.MolFromSmiles("CC"))) -- blank name line.
RDKIT_BLANK_NAME_MOL = (
    "\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
    "   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\nM  END\n"
)

MALFORMED_RECORD = "broken\n  prog\n\n  NOTNUM  0  0 V2000\nM  END\n"


def test_parse_sdf_with_coords_blank_name_middle_record():
    sdf = f"{MOL_A}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{MOL_A}$$$$\n"
    records = chematic.parse_sdf_with_coords(sdf)
    assert len(records) == 3
    _, name0, _ = records[0]
    _, name1, coords1 = records[1]
    _, name2, _ = records[2]
    assert name0 == "mol_a"
    assert name1 == ""
    assert len(coords1) == 2
    assert name2 == "mol_a"


def test_parse_sdf_with_coords_skips_malformed_records():
    # Malformed records are silently skipped, matching iter_sdf's behaviour --
    # not a regression introduced by delegating to SdfRecordReader.
    sdf = f"{MOL_A}$$$$\n{MALFORMED_RECORD}$$$$\n{MOL_A}$$$$\n"
    records = chematic.parse_sdf_with_coords(sdf)
    assert len(records) == 2
