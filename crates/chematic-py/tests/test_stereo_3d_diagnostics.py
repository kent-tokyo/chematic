import json
import math
from pathlib import Path

import chematic
import pytest


def test_selenium_stereo_reports_satisfied_violated_and_unevaluable():
    path = (Path(__file__).resolve().parents[2] / "chematic-3d" / "src" /
            "stereo_constraints_se_neutral_oracle.json")
    row = json.loads(path.read_text(encoding="utf-8"))["rows"][0]
    mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
    coords = row["coords"]
    for geometry, expected in ((coords, "satisfied"),
                               ([[-x, y, z] for x, y, z in coords], "violated"),
                               ([[0.0, 0.0, 0.0]] * len(coords), "unevaluable")):
        report = mol.stereo_3d_diagnostics(geometry)
        assert report["reported"] == 1
        assert report[expected] == 1
        assert report["atom_stereo"][0]["status"] == expected
        reason = report["atom_stereo"][0]["reason"]
        assert (reason is not None) == (expected == "unevaluable")
        assert report["double_bond_stereo"] == []


def test_double_bond_and_square_planar_are_reported_separately():
    mol = chematic.from_smiles("F/C=C/F")
    trans = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0]]
    cis = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]]
    assert mol.stereo_3d_diagnostics(trans)["double_bond_stereo"][0]["status"] == "satisfied"
    assert mol.stereo_3d_diagnostics(cis)["double_bond_stereo"][0]["status"] == "violated"
    pt = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    row = pt.stereo_3d_diagnostics(coords)["atom_stereo"][0]
    assert row == {"atom": 0, "status": "unevaluable", "reason": "UnsupportedCoordination"}


def test_stereo_diagnostics_empty_and_invalid_coordinates():
    mol = chematic.from_smiles("CCO")
    report = mol.stereo_3d_diagnostics([[0, 0, 0]] * 3)
    assert (report["reported"], report["violated"], report["unevaluable"]) == (0, 0, 0)
    with pytest.raises(ValueError, match="finite and match"):
        mol.stereo_3d_diagnostics([[0, 0, 0]])
    with pytest.raises(ValueError, match="finite and match"):
        mol.stereo_3d_diagnostics([[math.nan, 0, 0]] * 3)


def test_repair_stereo_3d_recovers_development_ammonium():
    smiles = "C[N@@H+](CCNc1ccc2c3c(cccc13)C(=O)N(CCOCCO)C2=O)Cc1ccccc1B(O)O"
    mol = chematic.from_smiles(smiles)
    before = mol.generate_3d_etkdg()
    assert mol.stereo_3d_diagnostics(before)["violated"] == 1
    result = mol.repair_stereo_3d(before)
    assert result["ok"] and result["fully_verified"]
    assert result["failures"] == [] and result["unevaluable"] == []
    assert result["repaired"][0]["kind"] == "atom"
    assert result["repaired"][0]["index"] == 1
    assert result["repaired"][0]["atoms_moved"] > 0
    assert result["repaired"][0]["max_displacement"] > 0
    assert mol.stereo_3d_diagnostics(result["coordinates"])["satisfied"] == 1


def test_repair_stereo_3d_reports_ring_failure_and_unknown_stereo():
    mol = chematic.from_smiles("C1[C@@H]2CCC1CC2")
    result = mol.repair_stereo_3d(mol.generate_3d_etkdg())
    assert result["ok"] is False and result["fully_verified"] is False
    assert result["coordinates"] is None
    assert result["failures"] == [{"kind": "atom", "index": 1,
                                  "reason": "NoBridgeEligibleSubstituent"}]
    pt = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    geometry = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    unknown = pt.repair_stereo_3d(geometry)
    assert unknown["ok"] and not unknown["fully_verified"]
    assert unknown["unevaluable"] == [{"kind": "atom", "index": 0}]
    with pytest.raises(ValueError, match="finite and match"):
        pt.repair_stereo_3d([[float("nan"), 0, 0]] * 5)


def test_repaired_ammonium_agrees_with_rdkit_stereo_from_coordinates():
    Chem = pytest.importorskip("rdkit.Chem")
    smiles = "C[N@@H+](CCNc1ccc2c3c(cccc13)C(=O)N(CCOCCO)C2=O)Cc1ccccc1B(O)O"
    mol = chematic.from_smiles(smiles)
    geometry = mol.generate_3d_etkdg()
    repaired = mol.repair_stereo_3d(geometry)["coordinates"]
    reference = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
    expected = reference.GetAtomWithIdx(1).GetProp("_CIPCode")
    codes = []
    for coords in (geometry, repaired):
        molecule = Chem.Mol(reference)
        conf = Chem.Conformer(molecule.GetNumAtoms())
        for idx, point in enumerate(coords):
            conf.SetAtomPosition(idx, point)
        molecule.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(molecule, replaceExistingTags=True)
        Chem.AssignStereochemistry(molecule, force=True, cleanIt=True)
        codes.append(molecule.GetAtomWithIdx(1).GetProp("_CIPCode"))
    assert codes[0] != expected and codes[1] == expected


def test_repair_stereo_rejects_new_clash_pairs_even_if_count_drops():
    # Exposed COD development record 532: 2 old clashes become 1 new clash.
    smiles = "C[C@@H]1OC(=O)N(C(=O)[C@H](Cc2ccccc2)NC(=O)OC(C)(C)C)[C@H]1C(=O)N[C@@H](CC(=O)OCc1ccccc1)c1ccccc1"
    mol = chematic.from_smiles(smiles)
    coords = mol.generate_3d_etkdg()
    legacy = mol.repair_stereo_3d(coords)
    assert legacy["ok"] and legacy["clash_check"] is None
    result = mol.repair_stereo_3d(coords, clash_cutoff=0.8)
    assert not result["ok"] and result["coordinates"] is None
    assert result["failures"] == [{"kind": "geometry", "reason": "NewNonbondedClashes"}]
    report = result["clash_check"]
    assert report["status"] == "rejected"
    assert (report["before_count"], report["after_count"]) == (2, 1)
    assert report["new"][0]["atoms"] == [3, 26]
    assert report["new"][0]["distance"] < 0.8
    assert len(report["resolved"]) == 2


def test_repair_stereo_clash_screen_success_failure_and_invalid_cutoff():
    mol = chematic.from_smiles("C[N@H+](CC)CCC")
    coords = mol.generate_3d_etkdg()
    passed = mol.repair_stereo_3d(coords, clash_cutoff=0.8)
    assert passed["ok"] and passed["clash_check"]["status"] == "passed"
    assert passed["clash_check"]["new"] == []
    for cutoff in (0, -1, math.nan, math.inf):
        with pytest.raises(ValueError, match="positive"):
            mol.repair_stereo_3d(coords, clash_cutoff=cutoff)
    ring = chematic.from_smiles("C1[C@@H]2CCC1CC2")
    failed = ring.repair_stereo_3d(ring.generate_3d_etkdg(), clash_cutoff=0.8)
    assert not failed["ok"] and failed["coordinates"] is None
    assert failed["clash_check"] == {"cutoff": 0.8, "status": "not_evaluated_repair_failed"}
