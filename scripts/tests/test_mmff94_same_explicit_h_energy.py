import pytest

from scripts.mmff94_same_explicit_h_energy import (
    evaluate,
    gradient_diagnostic,
    percentile_nearest_rank,
    rdkit_term_energy_breakdown,
    summarize,
    terminal,
)


class QuadraticEnergy:
    def mmff94_bounded_analytic_gradient(self, coords):
        return [[2.0 * value for value in point] for point in coords]

    def mmff94_energy_breakdown(self, coords):
        return {"total": sum(value * value for point in coords for value in point)}


def test_nearest_rank_percentile_is_deterministic():
    assert percentile_nearest_rank([], 9, 10) is None
    assert percentile_nearest_rank([3.0, 1.0, 2.0, 4.0], 9, 10) == 4.0


def test_summary_preserves_terminal_row_accounting():
    rows = [
        terminal(
            {"primary_category": "small"},
            0,
            "ok",
            abs_delta_kcal_mol=0.5,
        ),
        terminal({}, 1, "embed_failure"),
    ]
    result = summarize(rows)
    assert result["row_accounting"] == {
        "input_count": 2,
        "terminal_count": 2,
        "status_counts": {"ok": 1, "embed_failure": 1},
    }
    assert result["comparable_rows"] == 1
    assert result["within_1_kcal_mol"] == 1
    assert result["within_5_kcal_mol"] == 1
    assert result["gradient_diagnostic"] == {
        "rows": 0,
        "input_indices": [],
        "max_abs_error_kcal_mol_angstrom": None,
        "max_scaled_error": None,
    }


def test_gradient_diagnostic_matches_quadratic_energy():
    result = gradient_diagnostic(QuadraticEnergy(), [[1.0, -2.0, 0.5]], 1e-5)
    assert result["components"] == 3
    assert result["max_abs_error_kcal_mol_angstrom"] < 1e-9
    assert result["max_scaled_error"] < 1e-9


def test_summary_attributes_residuals_to_terms():
    terms = {key: 0.0 for key in ("bond", "angle", "stretch_bend", "oop", "torsion", "vdw", "electrostatic")}
    rdkit = dict(terms, angle=1.0, vdw=2.0)
    big = dict(terms, angle=3.5, oop=-0.05)
    rows = [
        terminal(
            {"id": "big"},
            0,
            "ok",
            abs_delta_kcal_mol=3.45,
            rdkit_energy_kcal_mol=3.0,
            rdkit_term_energies_kcal_mol=rdkit,
            term_delta_kcal_mol=big,
            dominant_term="angle",
        ),
        terminal(
            {"id": "small"},
            1,
            "ok",
            abs_delta_kcal_mol=0.01,
            rdkit_energy_kcal_mol=3.0,
            rdkit_term_energies_kcal_mol=rdkit,
            term_delta_kcal_mol=dict(terms, bond=0.01),
            dominant_term="bond",
        ),
    ]
    result = summarize(rows)
    assert result["per_term"]["angle"]["max_abs_delta_kcal_mol"] == 3.5
    assert result["per_term"]["angle"]["rows_above_1_kcal_mol"] == 1
    assert result["per_term"]["bond"]["rows_above_0_1_kcal_mol"] == 0
    assert result["rdkit_term_sum_vs_total_max_abs_kcal_mol"] == 0.0
    assert result["dominant_term_counts_above_1_kcal_mol"] == {"angle": 1}
    assert [row["id"] for row in result["rows_above_1_kcal_mol"]] == ["big"]


def test_rdkit_isolated_terms_sum_to_total_without_changing_coordinates():
    import pytest
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Includes an O-H donor, an aromatic ring, and a carbonyl; exercising
    # all seven terms rather than testing the setter loop against a mock.
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
    assert AllChem.EmbedMolecule(mol, randomSeed=43) == 0
    properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    force_field = AllChem.MMFFGetMoleculeForceField(mol, properties, confId=0)
    before = mol.GetConformer().GetPositions().copy()
    total = force_field.CalcEnergy()
    terms = rdkit_term_energy_breakdown(mol)
    assert set(terms) == {
        "bond", "angle", "stretch_bend", "torsion", "oop", "vdw", "electrostatic"
    }
    assert sum(terms.values()) == pytest.approx(total, abs=1e-9)
    assert force_field.CalcEnergy() == total
    assert (mol.GetConformer().GetPositions() == before).all()


def test_explicit_oxygen_hydrogen_vdw_and_charge_transfer_match_rdkit():
    import pytest
    import chematic
    from rdkit import Chem
    from rdkit.Chem import AllChem

    for smiles in ["NCC(=O)O", "ON(C)C", "Oc1ccccc1", "CC(=O)Nc1ccc(O)cc1"]:
        rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=43) == 0
        properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
        mol = chematic.from_smiles(smiles).add_hydrogens()
        # Check correspondence before sharing coordinates; matching element
        # lists alone would miss permutations of two atoms of the same element.
        assert [a["element"] for a in mol.depict_data()["atoms"]] == [
            a.GetSymbol() for a in rdkit_mol.GetAtoms()
        ]
        assert {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table} == {
            tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
            for b in rdkit_mol.GetBonds()
        }
        types = mol.mmff94_numeric_atom_types()
        charges = mol.mmff94_charges_numeric()
        for atom in rdkit_mol.GetAtoms():
            if atom.GetAtomicNum() == 1 and any(
                nb.GetAtomicNum() == 8 for nb in atom.GetNeighbors()
            ):
                i = atom.GetIdx()
                assert types[i] == properties.GetMMFFAtomType(i), smiles
                assert charges[i] == pytest.approx(properties.GetMMFFPartialCharge(i)), smiles
        coords = rdkit_mol.GetConformer().GetPositions().tolist()
        assert mol.mmff94_energy_breakdown(coords)["vdw"] == pytest.approx(
            rdkit_term_energy_breakdown(rdkit_mol)["vdw"], abs=1e-9
        ), smiles


def test_same_element_sequence_cannot_hide_an_indexed_graph_mismatch(monkeypatch):
    import chematic

    # Both isomers have C5H12 and the same AddHs element sequence. The
    # branching still makes sharing atom-indexed coordinates invalid.
    branched = chematic.from_smiles("CC(C)CC")
    monkeypatch.setattr(chematic, "from_smiles", lambda _: branched)
    result = evaluate({"smiles": "CCCCC"}, 0, 43, set(), 1e-5)
    assert result["status"] == "atom_order_mismatch"
    assert result["reason"] == "indexed_bond_graph_mismatch"


@pytest.mark.parametrize("smiles", [
    "c1ccccc1-c1ccccc1",
    "C=CC=C",
    "c1ccccc1C[n+]1ccc(N)c2ccccc21",
    "CC(=O)Nc1ccc(O)cc1",
    "N#CC",
])
def test_angle_and_stretch_bend_terms_match_rdkit_on_residual_motifs(smiles):
    import pytest
    import chematic
    from rdkit import Chem
    from rdkit.Chem import AllChem

    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=43) == 0
    # Exercise the linear-center cosine potential away from equilibrium;
    # a nearly straight embedded nitrile would conceal the wrong formula.
    if smiles == "N#CC":
        rdkit_mol.GetConformer().SetAtomPosition(0, (0.0, 1.1, 0.0))
        rdkit_mol.GetConformer().SetAtomPosition(1, (0.0, 0.0, 0.0))
        rdkit_mol.GetConformer().SetAtomPosition(2, (1.3, -0.4, 0.0))
    mol = chematic.from_smiles(smiles).add_hydrogens()
    assert [a["element"] for a in mol.depict_data()["atoms"]] == [
        a.GetSymbol() for a in rdkit_mol.GetAtoms()
    ]
    assert {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table} == {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        for b in rdkit_mol.GetBonds()
    }
    coords = rdkit_mol.GetConformer().GetPositions().tolist()
    actual = mol.mmff94_energy_breakdown(coords)
    expected = rdkit_term_energy_breakdown(rdkit_mol)
    for term in ["angle", "stretch_bend"]:
        assert actual[term] == pytest.approx(expected[term], abs=1e-8), (smiles, term)


@pytest.mark.parametrize("smiles,twist", [
    ("C1CC1", None), ("C1CO1", None), ("C1CN1", None), ("CC1CC1", None),
    ("C1CCC1", None), ("CCCC", None),
    ("C=[N+](C)c1ccccc1", (0, 1, 3, 4)),
    ("CC=[N+](C)c1ccccc1", (1, 2, 4, 5)),
])
def test_torsion_terms_match_rdkit_for_small_rings_and_iminium(smiles, twist):
    import chematic
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=43) == 0
    if twist is not None:
        # Embedded iminium groups are nearly planar: move away from the
        # minimum so a missing conjugated barrier cannot hide at zero energy.
        rdMolTransforms.SetDihedralDeg(rdkit_mol.GetConformer(), *twist, 70.0)
    mol = chematic.from_smiles(smiles).add_hydrogens()
    properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol)
    assert list(mol.mmff94_numeric_atom_types()) == [
        properties.GetMMFFAtomType(i)
        for i in range(rdkit_mol.GetNumAtoms())
    ]
    assert {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table} == {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        for b in rdkit_mol.GetBonds()
    }
    coords = rdkit_mol.GetConformer().GetPositions().tolist()
    actual = mol.mmff94_energy_breakdown(coords)
    expected = rdkit_term_energy_breakdown(rdkit_mol)
    assert actual["torsion"] == pytest.approx(expected["torsion"], abs=1e-8)
    gradient = gradient_diagnostic(mol, coords, 1e-5)
    assert gradient["max_scaled_error"] < 1e-5


@pytest.mark.parametrize("smiles,displaced", [
    ("C=C(C)C", 1), ("O=C(N)c1ccccc1", 1),
    ("C=CC=[N+](C)C", 2), ("CN(C)C=O", 1),
    ("[O-][n+]1ccccc1", 1),
])
def test_oop_energy_matches_rdkit_away_from_planarity(smiles, displaced):
    import chematic
    from rdkit import Chem
    from rdkit.Chem import AllChem

    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=43) == 0
    properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol)
    position = rdkit_mol.GetConformer().GetAtomPosition(displaced)
    rdkit_mol.GetConformer().SetAtomPosition(displaced, (position.x, position.y, position.z + 0.4))
    mol = chematic.from_smiles(smiles).add_hydrogens()
    assert list(mol.mmff94_numeric_atom_types()) == [
        properties.GetMMFFAtomType(i) for i in range(rdkit_mol.GetNumAtoms())
    ]
    assert {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table} == {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in rdkit_mol.GetBonds()
    }
    coords = rdkit_mol.GetConformer().GetPositions().tolist()
    assert mol.mmff94_energy_breakdown(coords)["oop"] == pytest.approx(
        rdkit_term_energy_breakdown(rdkit_mol)["oop"], abs=1e-9
    )
    assert gradient_diagnostic(mol, coords, 1e-5)["max_scaled_error"] < 1e-5
    # Directly compare the complete force-field gradient as well as the
    # internal central difference, on these deliberately nonplanar motifs.
    force_field = AllChem.MMFFGetMoleculeForceField(rdkit_mol, properties, confId=0)
    expected = force_field.CalcGrad()
    actual = [value for point in mol.mmff94_bounded_analytic_gradient(coords) for value in point]
    assert max(abs(a - b) / (1.0 + abs(b)) for a, b in zip(actual, expected, strict=True)) < 1e-7


def test_vdw_energy_retains_pairs_beyond_ten_angstroms():
    import chematic
    from rdkit import Chem
    from rdkit.Chem import AllChem

    smiles = "CCCCCCCCCCCC"
    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=43) == 0
    mol = chematic.from_smiles(smiles).add_hydrogens()
    assert {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table} == {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in rdkit_mol.GetBonds()
    }
    # Enlarge the existing embedding uniformly to ensure nonexcluded distant
    # pairs exist without creating disconnected fragments or atom reordering.
    coords = (rdkit_mol.GetConformer().GetPositions() * 2.0).tolist()
    for index, position in enumerate(coords):
        rdkit_mol.GetConformer().SetAtomPosition(index, position)
    assert mol.mmff94_energy_breakdown(coords)["vdw"] == pytest.approx(
        rdkit_term_energy_breakdown(rdkit_mol)["vdw"], abs=1e-9
    )


@pytest.mark.parametrize("smiles", [
    "O", "NC(=[NH2+])N", "[NH2+]=C(N)N", "CNC(=[NH2+])N",
    "CN(C)C(=[N+](C)C)N(C)C", "NNC(=[NH2+])N", "N=C(N)N",
    "CC(=[NH2+])N", "CNC=[NH2+]", "NC(=[NH2+])NC(=[NH2+])N",
    "N=C(N)NC(=[NH2+])N", "CO", "COC", "[OH-]", "[OH3+]",
    "C=[OH+]", "C=[NH2+]", "C[NH3+]", "NC(=[NH2+])NCCCCC",
    "NC(=[NH2+])[N+](=O)[O-]",
])
def test_cation_water_holdout_types_charges_energy_and_gradient_match_rdkit(smiles):
    import chematic
    from rdkit import Chem
    from rdkit.Chem import AllChem

    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=43) == 0
    properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
    mol = chematic.from_smiles(smiles).add_hydrogens()
    assert [a["element"] for a in mol.depict_data()["atoms"]] == [
        a.GetSymbol() for a in rdkit_mol.GetAtoms()
    ]
    assert {(min(a, b), max(a, b), order) for a, b, order, *_ in mol.bond_table} == {
        (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
         max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType()))
        for b in rdkit_mol.GetBonds()
    }
    assert list(mol.mmff94_numeric_atom_types()) == [
        properties.GetMMFFAtomType(i) for i in range(rdkit_mol.GetNumAtoms())
    ]
    assert list(mol.mmff94_charges_numeric()) == pytest.approx([
        properties.GetMMFFPartialCharge(i) for i in range(rdkit_mol.GetNumAtoms())
    ], abs=1e-12)
    assert sum(mol.mmff94_charges_numeric()) == pytest.approx(
        sum(atom.GetFormalCharge() for atom in rdkit_mol.GetAtoms()), abs=1e-12
    )
    # Perturb every coordinate deterministically; ideal planar cations and
    # equilibrium water angles must not hide a missing force-field term.
    coords = rdkit_mol.GetConformer().GetPositions().tolist()
    for index, point in enumerate(coords):
        point[2] += 0.03 * ((index % 5) - 2)
        rdkit_mol.GetConformer().SetAtomPosition(index, point)
    actual = mol.mmff94_energy_breakdown(coords)
    expected = rdkit_term_energy_breakdown(rdkit_mol)
    for term, energy in expected.items():
        assert actual[term] == pytest.approx(energy, abs=1e-8), (smiles, term)
    force_field = AllChem.MMFFGetMoleculeForceField(rdkit_mol, properties, confId=0)
    assert actual["total"] == pytest.approx(force_field.CalcEnergy(), abs=1e-8)
    gradient = [value for point in mol.mmff94_bounded_analytic_gradient(coords) for value in point]
    reference = force_field.CalcGrad()
    assert max(abs(a - b) / (1.0 + abs(b)) for a, b in zip(gradient, reference, strict=True)) < 1e-7
    assert gradient_diagnostic(mol, coords, 1e-5)["max_scaled_error"] < 1e-5
