import math

import pytest

import chematic


def assert_same_indexed_graph(mol, reference):
    assert [row[1] for row in mol.atom_table] == [atom.GetAtomicNum() for atom in reference.GetAtoms()]
    bonds = {tuple(sorted((row[0], row[1]))) for row in mol.bond_table}
    reference_bonds = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
                       for bond in reference.GetBonds()}
    assert bonds == reference_bonds


def test_square_planar_geometry_reports_ideal_and_degenerate_coordinates():
    mol = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    report = mol.square_planar_geometry_diagnostics(coords)[0]
    assert report["center"] == 0
    assert report["permutation"] == "SP1"
    assert report["neighbors"] == [1, 2, 3, 4]
    assert report["unavailable"] is None
    assert {(a, b) for a, b, _ in report["trans"]} == {(1, 3), (2, 4)}
    assert all(angle == pytest.approx(180.0) for _, _, angle in report["trans"])
    assert all(angle == pytest.approx(90.0) for _, _, angle in report["cis"])
    report = mol.square_planar_geometry_diagnostics([[0, 0, 0]] * 5)[0]
    assert report["unavailable"] == "degenerate_ligand_vector"
    assert report["trans"] == report["cis"] == []
    with pytest.raises(ValueError, match="finite and match"):
        mol.square_planar_geometry_diagnostics(coords[:-1])
    with pytest.raises(ValueError, match="finite and match"):
        mol.square_planar_geometry_diagnostics([[math.nan, 0, 0]] + coords[1:])


@pytest.mark.parametrize("pattern", ["[Pt@SP{}](F)(Cl)(Br)I", "F[Pt@SP{}](Cl)Br"])
@pytest.mark.parametrize("permutation", [1, 2, 3])
def test_square_planar_embedding_enforces_final_angle_policy(pattern, permutation):
    mol = chematic.from_smiles(pattern.format(permutation)).add_hydrogens()
    coords = None
    failures = []
    for seed in range(1, 11):
        try:
            coords = mol.rdkit_embed_square_planar(random_seed=seed)
            break
        except RuntimeError as error:
            failures.append(str(error))
    assert coords is not None, failures
    report = mol.square_planar_geometry_diagnostics(coords)[0]
    assert report["unavailable"] is None
    assert all(angle >= 150.0 for _, _, angle in report["trans"])
    assert all(60.0 <= angle <= 120.0 for _, _, angle in report["cis"])


def test_square_planar_embedding_rejects_invalid_or_unsupported_input():
    with pytest.raises(ValueError, match="no square-planar declaration"):
        chematic.from_smiles("CCO").rdkit_embed_square_planar()
    with pytest.raises(ValueError, match="unnumbered_permutation"):
        chematic.from_smiles("[Pt@SP](F)(Cl)(Br)I").rdkit_embed_square_planar()
    with pytest.raises(ValueError, match="unsupported_neighbor_order"):
        chematic.from_smiles("F[Pt@SP3]Cl").rdkit_embed_square_planar()
    mol = chematic.from_smiles("F[Pt@SP3](Cl)Br")
    for trans, cis in [(120.0, 30.0), (float("nan"), 30.0), (150.0, -1.0)]:
        with pytest.raises(ValueError, match="invalid or overlapping"):
            mol.rdkit_embed_square_planar(trans_min_degrees=trans, cis_tolerance_degrees=cis)


@pytest.mark.parametrize("pattern", ["[Pt@SP{}](F)(Cl)(Br)I", "F[Pt@SP{}](Cl)(Br)I"])
@pytest.mark.parametrize("permutation", [1, 2, 3])
def test_square_planar_angles_match_rdkit_coordinates(pattern, permutation):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem, rdMolTransforms

    smiles = pattern.format(permutation)
    reference = Chem.MolFromSmiles(smiles)
    assert AllChem.EmbedMolecule(reference, randomSeed=42) == 0
    center = next(atom.GetIdx() for atom in reference.GetAtoms() if atom.GetSymbol() == "Pt")
    conformer = reference.GetConformer()
    coords = [list(conformer.GetAtomPosition(i)) for i in range(reference.GetNumAtoms())]
    mol = chematic.from_smiles(smiles)
    assert_same_indexed_graph(mol, reference)
    report = mol.square_planar_geometry_diagnostics(coords)[0]
    assert report["center"] == center
    assert report["permutation"] == f"SP{permutation}"
    assert report["unavailable"] is None
    assert len(report["trans"]) == 2 and len(report["cis"]) == 4
    for a, b, angle in report["trans"] + report["cis"]:
        assert angle == pytest.approx(rdMolTransforms.GetAngleDeg(conformer, a, center, b), abs=1e-10)
    # This is an oracle check for these six fixed structures, not a general geometry cutoff.
    assert min(angle for _, _, angle in report["trans"]) > 150
    assert max(angle for _, _, angle in report["cis"]) < 120


@pytest.mark.parametrize("smiles", [
    "F[Pt@SP3](Cl)Br",
    "O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F",
    "Brc1ccc2c3c(cccc13)[Se@@](c1ccccc1)=[Se@SP3]2Br",
])
def test_square_planar_vacancy_measures_only_real_ligands(smiles):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem, rdMolTransforms

    reference = Chem.MolFromSmiles(smiles)
    assert AllChem.EmbedMolecule(reference, randomSeed=42) == 0
    conf = reference.GetConformer()
    coords = [list(conf.GetAtomPosition(i)) for i in range(reference.GetNumAtoms())]
    mol = chematic.from_smiles(smiles)
    assert_same_indexed_graph(mol, reference)
    report = mol.square_planar_geometry_diagnostics(coords)[0]
    assert report["unavailable"] is None
    assert report["neighbors"].count(None) == 1
    assert len(report["trans"]) == 1 and len(report["cis"]) == 2
    for a, b, angle in report["trans"] + report["cis"]:
        assert angle == pytest.approx(rdMolTransforms.GetAngleDeg(conf, a, report["center"], b), abs=1e-10)


def test_square_planar_retry_preserves_first_success():
    mol = chematic.from_smiles("[Pt@SP3](F)(Cl)(Br)I").add_hydrogens()
    for seed in range(1, 11):
        try:
            first = mol.rdkit_embed_square_planar(random_seed=seed)
            break
        except RuntimeError:
            continue
    else:
        pytest.fail("no initial embedding met the angle policy")
    result = mol.rdkit_embed_square_planar_with_retry(random_seed=seed)
    assert result["attempts"] == 1 and result["used_tightened_bounds"] is False
    assert result["initial_rejection"] is None
    assert result["coordinates"] == first


def test_square_planar_retry_recovers_rejected_embedding():
    mol = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I").add_hydrogens()
    for seed in range(1, 16):
        try:
            mol.rdkit_embed_square_planar(random_seed=seed)
        except RuntimeError:
            try:
                result = mol.rdkit_embed_square_planar_with_retry(random_seed=seed)
                break
            except RuntimeError:
                continue
    else:
        pytest.fail("no retry recovered a rejected embedding")
    assert result["attempts"] == 2 and result["used_tightened_bounds"] is True
    assert "outside the selected limits" in result["initial_rejection"]
    report = mol.square_planar_geometry_diagnostics(result["coordinates"])[0]
    assert all(angle >= 150.0 for _, _, angle in report["trans"])
    assert all(60.0 <= angle <= 120.0 for _, _, angle in report["cis"])
    with pytest.raises(ValueError, match="no square-planar declaration"):
        chematic.from_smiles("CCO").rdkit_embed_square_planar_with_retry()


def test_square_planar_nonplanarity_can_exceed_the_angle_screen():
    mol = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    angle = math.radians(20)
    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0],
              [-math.cos(angle), 0, math.sin(angle)], [0, -math.cos(angle), math.sin(angle)]]
    report = mol.square_planar_geometry_diagnostics(coords)[0]
    assert all(value >= 150 for _, _, value in report["trans"])
    assert all(60 <= value <= 120 for _, _, value in report["cis"])
    assert report["max_unit_triple_product"] == pytest.approx(math.sin(angle), abs=1e-12)


def test_square_planar_checked_embedding_requires_explicit_limit_and_keeps_coordinates():
    mol = chematic.from_smiles("[Pt@SP3](F)(Cl)(Br)I").add_hydrogens()
    with pytest.raises(TypeError):
        mol.rdkit_embed_square_planar_checked()
    for limit in (-0.1, 1.1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="between 0 and 1"):
            mol.rdkit_embed_square_planar_checked(max_unit_triple_product=limit)
    for seed in range(1, 11):
        try:
            base = mol.rdkit_embed_square_planar_with_retry(random_seed=seed)
            report = mol.square_planar_geometry_diagnostics(base["coordinates"])[0]
            if report["max_unit_triple_product"] > 1e-9:
                break
        except RuntimeError:
            continue
    else:
        pytest.fail("no evaluable nonplanar embedding")
    checked = mol.rdkit_embed_square_planar_checked(max_unit_triple_product=1.0, random_seed=seed)
    assert checked["coordinates"] == base["coordinates"]
    assert checked["attempts"] == base["attempts"]
    assert checked["other_stereo_checked"] == 0
    assert checked["max_unit_triple_product"] == pytest.approx(report["max_unit_triple_product"])
    with pytest.raises(RuntimeError, match="unit triple product"):
        mol.rdkit_embed_square_planar_checked(
            max_unit_triple_product=report["max_unit_triple_product"] / 2, random_seed=seed
        )


@pytest.mark.parametrize("smiles", [
    "[Pt@SP3](F)(Cl)(Br)I.C[C@H](O)F",
    "F[Pt@SP1](Cl)Br.F/C=C/F",
])
def test_square_planar_checked_embedding_checks_other_declared_stereo(smiles):
    mol = chematic.from_smiles(smiles).add_hydrogens()
    errors = []
    for seed in range(1, 11):
        try:
            result = mol.rdkit_embed_square_planar_checked(
                random_seed=seed, max_unit_triple_product=1.0
            )
            break
        except RuntimeError as error:
            errors.append(str(error))
    else:
        pytest.fail(f"no checked embedding: {errors}")
    assert result["other_stereo_checked"] == 1
    assert result["coordinates"] == mol.rdkit_embed_square_planar_with_retry(random_seed=seed)["coordinates"]
