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


def test_alternative_stereo_repair_recovers_checked_development_case():
    smiles = "CC[C@H](NS(=O)(=O)c1ccc(C)cc1)c1n[nH]c(=S)o1"
    mol = chematic.from_smiles(smiles)
    coords = mol.generate_3d_etkdg()
    baseline = mol.repair_stereo_3d(coords, clash_cutoff=0.8)
    assert not baseline["ok"] and baseline["clash_check"]["status"] == "rejected"
    result = mol.repair_stereo_3d(coords, clash_cutoff=0.8, max_alternative_attempts=3)
    assert result["ok"] and result["fully_verified"]
    assert result["retry"]["attempted"] == 2
    assert result["retry"]["selected_variant"] == 2
    assert result["retry"]["clash_rejections"] == 1
    assert result["clash_check"]["new"] == []
    assert mol.stereo_3d_diagnostics(result["coordinates"])["violated"] == 0
    repeat = mol.repair_stereo_3d(coords, clash_cutoff=0.8, max_alternative_attempts=3)
    assert repeat == result
    for a, b, *_ in mol.bond_table:
        assert abs(math.dist(coords[a], coords[b]) -
                   math.dist(result["coordinates"][a], result["coordinates"][b])) < 1e-10


def test_alternative_stereo_repair_rejects_distorted_ligand_geometry():
    smiles = "O=C1OC(c2ccccc2)=N[C@]1(Cc1ccccc1)SCc1ccccc1"
    mol = chematic.from_smiles(smiles)
    coords = mol.generate_3d_etkdg()
    result = mol.repair_stereo_3d(coords, clash_cutoff=0.8, max_alternative_attempts=3)
    assert not result["ok"] and result["coordinates"] is None
    assert result["retry"]["tetrahedral_geometry_rejections"] == 3
    with pytest.raises(ValueError, match="clash_cutoff is required"):
        mol.repair_stereo_3d(coords, max_alternative_attempts=1)
    with pytest.raises(ValueError, match="at most 16"):
        mol.repair_stereo_3d(coords, clash_cutoff=0.8, max_alternative_attempts=17)


def test_alternative_stereo_repair_agrees_with_rdkit_on_recovered_case():
    Chem = pytest.importorskip("rdkit.Chem")
    smiles = "CC[C@H](NS(=O)(=O)c1ccc(C)cc1)c1n[nH]c(=S)o1"
    mol = chematic.from_smiles(smiles)
    result = mol.repair_stereo_3d(mol.generate_3d_etkdg(), clash_cutoff=0.8,
                                  max_alternative_attempts=3)
    assert result["ok"]
    reference = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
    expected = reference.GetAtomWithIdx(2).GetProp("_CIPCode")
    conf = Chem.Conformer(reference.GetNumAtoms())
    for index, xyz in enumerate(result["coordinates"]):
        conf.SetAtomPosition(index, xyz)
    reference.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(reference, replaceExistingTags=True)
    Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
    assert reference.GetAtomWithIdx(2).GetProp("_CIPCode") == expected


def test_distorted_prototype_illustrates_verifier_oracle_boundary():
    Chem = pytest.importorskip("rdkit.Chem")
    fixture = json.loads(Path(__file__).with_name("stereo_distorted_prototype_719.json").read_text())
    mol = chematic.from_smiles(fixture["smiles"])
    assert mol.stereo_3d_diagnostics(fixture["coords"])["atom_stereo"][0]["status"] == "satisfied"
    for geometry, raw in ((fixture["coords"], "satisfied"),
                          ([[-x, y, z] for x, y, z in fixture["coords"]], "violated")):
        assert mol.stereo_3d_diagnostics(geometry)["atom_stereo"][0]["status"] == raw
        screened = mol.stereo_3d_diagnostics(geometry, require_regular_geometry=True)
        assert screened["atom_stereo"][0] == {"atom": fixture["center"],
            "status": "unevaluable", "reason": "DistortedTetrahedralGeometry"}
        assert (screened["satisfied"], screened["violated"], screened["unevaluable"]) == (0, 0, 1)
    reference = Chem.MolFromSmiles(fixture["smiles"])
    Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
    index = fixture["center"]
    expected = reference.GetAtomWithIdx(index).GetProp("_CIPCode")
    conf = Chem.Conformer(reference.GetNumAtoms())
    for i, point in enumerate(fixture["coords"]):
        conf.SetAtomPosition(i, point)
    reference.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(reference, replaceExistingTags=True)
    Chem.AssignStereochemistry(reference, force=True, cleanIt=True)
    assert reference.GetAtomWithIdx(index).GetProp("_CIPCode") != expected


def test_tetrahedral_geometry_screen_preserves_regular_and_double_bond_results():
    smiles = "CC[C@H](NS(=O)(=O)c1ccc(C)cc1)c1n[nH]c(=S)o1"
    mol = chematic.from_smiles(smiles)
    repair = mol.repair_stereo_3d(mol.generate_3d_etkdg(), clash_cutoff=0.8,
                                  max_alternative_attempts=3)
    assert repair["ok"]
    coords = repair["coordinates"]
    assert mol.stereo_3d_diagnostics(coords, require_regular_geometry=True) == mol.stereo_3d_diagnostics(coords)
    mirror = [[-x, y, z] for x, y, z in coords]
    assert mol.stereo_3d_diagnostics(mirror, require_regular_geometry=True) == mol.stereo_3d_diagnostics(mirror)
    alkene = chematic.from_smiles("F/C=C/F")
    geometry = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0]]
    assert alkene.stereo_3d_diagnostics(geometry, require_regular_geometry=True) == alkene.stereo_3d_diagnostics(geometry)
    with pytest.raises(ValueError, match="finite and match"):
        mol.stereo_3d_diagnostics([[math.nan, 0, 0]] * len(coords), require_regular_geometry=True)


def test_near_planar_implicit_h_carbon_is_not_confidently_assigned():
    mol = chematic.from_smiles("F[C@H](Cl)Br")
    coords = [[-0.8770239853135249, 0.4924418600749846, 0.9596809511921262],
              [0.0, 0.0, 0.0],
              [-0.2096233151177489, 0.7302430680925626, -0.6898994627180793],
              [0.6835160275334422, -0.6881890778227236, -0.2249292230318938]]
    assert mol.stereo_3d_diagnostics(coords)["atom_stereo"][0]["status"] == "violated"
    strict = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
    assert strict["atom_stereo"][0] == {"atom": 1, "status": "unevaluable",
                                         "reason": "DistortedTetrahedralGeometry"}
    for scale in (0.5, 3.0):
        transformed = [[scale * -x + 7.0, scale * y - 5.0, scale * z + 1.0]
                       for x, y, z in coords]
        assert mol.stereo_3d_diagnostics(transformed, require_regular_geometry=True)["atom_stereo"][0]["status"] == "unevaluable"
    regular = [[1., 1., 1.], [0., 0., 0.], [-1., -1., 1.], [-1., 1., -1.]]
    assert mol.stereo_3d_diagnostics(regular, require_regular_geometry=True) == mol.stereo_3d_diagnostics(regular)
    Chem = pytest.importorskip("rdkit.Chem")
    reference = Chem.MolFromSmiles("F[C@H](Cl)Br")
    conf = Chem.Conformer(reference.GetNumAtoms())
    for index, point in enumerate(coords):
        conf.SetAtomPosition(index, point)
    reference.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(reference, replaceExistingTags=True)
    assert reference.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED


def test_near_planar_carbon_screen_agrees_for_explicit_and_implicit_hydrogen():
    implicit = chematic.from_smiles("F[C@H](Cl)Br")
    explicit = chematic.from_smiles("F[C@]([H])(Cl)Br")
    coords = [[-0.8770239853135249, 0.4924418600749846, 0.9596809511921262],
              [0., 0., 0.],
              [-0.2096233151177489, 0.7302430680925626, -0.6898994627180793],
              [0.6835160275334422, -0.6881890778227236, -0.2249292230318938]]
    vectors = [coords[i] for i in (0, 2, 3)]
    unit = [[v / math.sqrt(sum(x * x for x in vector)) for v in vector] for vector in vectors]
    opposite = [-sum(v[i] for v in unit) for i in range(3)]
    norm = math.sqrt(sum(x * x for x in opposite))
    hydrogen = [x / norm for x in opposite]
    explicit_coords = [coords[0], coords[1], hydrogen, coords[2], coords[3]]
    for mol, geometry in ((implicit, coords), (explicit, explicit_coords)):
        assert mol.stereo_3d_diagnostics(geometry)["atom_stereo"][0]["status"] == "violated"
        screened = mol.stereo_3d_diagnostics(geometry, require_regular_geometry=True)
        assert screened["atom_stereo"][0]["status"] == "unevaluable"
        assert screened["atom_stereo"][0]["reason"] == "DistortedTetrahedralGeometry"

    # The extra condition is limited to H directions within the documented
    # cosine tolerance of the reconstructed direction.
    for degrees, expected in ((0.001, "unevaluable"), (0.1, "violated")):
        angle = math.radians(degrees)
        hx, hy, hz = hydrogen
        rotated = [hx * math.cos(angle) - hy * math.sin(angle),
                   hx * math.sin(angle) + hy * math.cos(angle), hz]
        geometry = [coords[0], coords[1], rotated, coords[2], coords[3]]
        screened = explicit.stereo_3d_diagnostics(geometry, require_regular_geometry=True)
        assert screened["atom_stereo"][0]["status"] == expected

    Chem = pytest.importorskip("rdkit.Chem")
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    for smiles, geometry in (("F[C@H](Cl)Br", coords),
                             ("F[C@]([H])(Cl)Br", explicit_coords)):
        ref = Chem.MolFromSmiles(smiles, parser)
        conf = Chem.Conformer(ref.GetNumAtoms())
        for index, point in enumerate(geometry):
            conf.SetAtomPosition(index, point)
        ref.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
        assert ref.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED


def test_independently_placed_hydrogen_is_not_treated_as_reconstructed():
    mol = chematic.from_smiles("F[C@]([H])(Cl)Br")
    coords = [[-0.31716336694990643, -1.1529430890295795, 0.8131938328381386],
              [0., 0., 0.],
              [0.4044183324438519, 0.48023863138805034, 1.248358878123955],
              [0.6088488519188352, 0.646845493145074, 0.15207343144295907],
              [-0.6224489427854606, -0.38682559639918296, -0.5246302514739701]]
    screened = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
    assert screened["atom_stereo"][0]["status"] == "violated"
    Chem = pytest.importorskip("rdkit.Chem")
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    ref = Chem.MolFromSmiles("F[C@]([H])(Cl)Br", parser)
    expected = ref.GetAtomWithIdx(1).GetChiralTag()
    conf = Chem.Conformer(ref.GetNumAtoms())
    for index, point in enumerate(coords):
        conf.SetAtomPosition(index, point)
    ref.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
    assert ref.GetAtomWithIdx(1).GetChiralTag() != expected
    assert ref.GetAtomWithIdx(1).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED


def test_repair_can_reject_distorted_success_without_returning_coordinates():
    mol = chematic.from_smiles("O=C1OC(c2ccccc2)=N[C@]1(Cc1ccccc1)SCc1ccccc1")
    coords = mol.generate_3d_etkdg()
    default = mol.repair_stereo_3d(coords)
    assert default["ok"] and default["regular_geometry_check"] is None
    assert mol.stereo_3d_diagnostics(default["coordinates"], require_regular_geometry=True)["unevaluable"] == 1
    result = mol.repair_stereo_3d(coords, require_regular_geometry=True)
    assert not result["ok"] and not result["fully_verified"] and result["coordinates"] is None
    assert result["regular_geometry_check"] == {"status": "rejected", "distorted_atoms": [11]}
    assert result["failures"] == [{"kind": "atom", "index": 11, "reason": "DistortedTetrahedralGeometry"}]
    assert result["retry"]["selected_variant"] is None


def test_repair_regular_geometry_screen_success_and_failure():
    mol = chematic.from_smiles("C[N@H+](CC)CCC")
    result = mol.repair_stereo_3d(mol.generate_3d_etkdg(), require_regular_geometry=True)
    assert result["ok"] and result["fully_verified"]
    assert result["regular_geometry_check"] == {"status": "passed", "distorted_atoms": []}
    ring = chematic.from_smiles("C1[C@@H]2CCC1CC2")
    result = ring.repair_stereo_3d(ring.generate_3d_etkdg(), require_regular_geometry=True)
    assert not result["ok"] and result["coordinates"] is None
    assert result["regular_geometry_check"] == {"status": "not_evaluated_repair_failed"}
    pt = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    result = pt.repair_stereo_3d(coords, require_regular_geometry=True)
    assert result["ok"] and not result["fully_verified"]
    assert result["regular_geometry_check"]["status"] == "passed"
    assert result["unevaluable"] == [{"kind": "atom", "index": 0}]


def test_repair_retries_when_only_regular_geometry_check_rejects_baseline():
    # Exposed COD development record 2110: baseline has no new clash, but
    # variant 2 fixes its distorted tetrahedral geometry.
    smiles = "CC(=O)Oc1cc2c(s1)CCN([C@H](C(=O)C1CC1)c1ccccc1F)C2"
    mol = chematic.from_smiles(smiles)
    coords = mol.generate_3d_etkdg()
    baseline = mol.repair_stereo_3d(coords, clash_cutoff=0.8)
    assert baseline["ok"] and baseline["clash_check"]["status"] == "passed"
    rejected = mol.repair_stereo_3d(coords, clash_cutoff=0.8, require_regular_geometry=True)
    assert not rejected["ok"] and rejected["regular_geometry_check"]["status"] == "rejected"
    result = mol.repair_stereo_3d(coords, clash_cutoff=0.8, max_alternative_attempts=3,
                                  require_regular_geometry=True)
    assert result["ok"] and result["fully_verified"]
    assert result["retry"]["attempted"] == 2 and result["retry"]["selected_variant"] == 2
    assert result["regular_geometry_check"] == {"status": "passed", "distorted_atoms": []}
    assert result["clash_check"]["status"] == "passed"
    assert mol.stereo_3d_diagnostics(result["coordinates"], require_regular_geometry=True)["unevaluable"] == 0


def test_repair_reports_clash_and_tetrahedral_distortion_together():
    mol = chematic.from_smiles("O=C1OC(c2ccccc2)=N[C@]1(Cc1ccccc1)SCc1ccccc1")
    result = mol.repair_stereo_3d(mol.generate_3d_etkdg(), clash_cutoff=0.8,
                                  max_alternative_attempts=3, require_regular_geometry=True)
    assert not result["ok"] and result["coordinates"] is None
    assert result["clash_check"]["status"] == "rejected"
    assert result["regular_geometry_check"] == {"status": "rejected", "distorted_atoms": [11]}
    assert result["failures"] == [
        {"kind": "geometry", "reason": "NewNonbondedClashes"},
        {"kind": "atom", "index": 11, "reason": "DistortedTetrahedralGeometry"}]
    assert result["retry"]["attempted"] == 3
    assert result["retry"]["tetrahedral_geometry_rejections"] == 3
