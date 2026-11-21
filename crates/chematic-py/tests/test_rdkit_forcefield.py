"""Mol.rdkit_mmff_* / rdkit_uff_* against RDKit 2026.03.1 force fields.

Reference values: ``Chem.AddHs(MolFromSmiles(SMILES))`` embedded with
``EmbedMolecule(randomSeed=42)``; MMFF via ``MMFFGetMoleculeForceField(m,
MMFFGetMoleculeProperties(m, mmffVariant=...), nonBondedThresh=100,
ignoreInterfragInteractions=True)``, UFF via ``UFFGetMoleculeForceField(m,
vdwThresh=10, ignoreInterfragInteractions=True)``, then ``Initialize();
Minimize(maxIts=200); CalcEnergy()``.

The recorded Linux comparison is bit-identical. Portable assertions allow
the same host-libm/minimizer tolerances as the Rust tests.
"""

import json
import math
from pathlib import Path

import pytest

import chematic

SMILES = "CC(=O)Nc1ccc(O)cc1"
COORDS = [
    [3.762175028044453, 0.30538563337732244, -0.09745815343294277],
    [2.3329369300899403, 0.11610000936619912, -0.5176254421062628],
    [2.1306032462731035, -0.09593752049575305, -1.7310611678300563],
    [1.3484147775946786, 0.19222362621777822, 0.4794175681721041],
    [-0.04651062589597437, 0.03729534973500948, 0.2549743150386899],
    [-0.6601760831775039, -1.210120317536445, 0.3277557282622417],
    [-2.027875532733807, -1.3550238887272026, 0.10667566446336685],
    [-2.822456517029374, -0.2670983319336878, -0.1922446077538659],
    [-4.187714927031642, -0.3854347185261268, -0.416753109003994],
    [-2.215264059602086, 0.990686079372787, -0.2676816459777285],
    [-0.8506489653837782, 1.1179323355827404, -0.04458045101770329],
    [3.8141231841990857, 1.3434050655543417, 0.3024234778452463],
    [4.477057254518884, 0.10425801972676844, -0.9145679992818092],
    [3.909133114355366, -0.42656371299342444, 0.7351586531314641],
    [1.6266288380623064, 0.37726586194308837, 1.4846206751064692],
    [-0.03711779601633895, -2.0820010342352497, 0.5648832746587801],
    [-2.5130630664729114, -2.314409067782831, 0.16030971994515217],
    [-4.834672252546139, -0.3986600801567403, 0.37067475360645646],
    [-2.8406901294426152, 1.8542155191372227, -0.5040194618516627],
    [-0.3648824178056287, 2.096481172374245, -0.10090179197395394],
]


@pytest.fixture
def mol():
    return chematic.from_smiles(SMILES).add_hydrogens()


@pytest.mark.parametrize(
    "variant, energy, opt_energy, x0",
    [
        ("MMFF94", 22.565180175565263, -12.775974493568183, 3.6705868960020265),
        ("MMFF94s", 22.959370426002092, -11.453744586799015, 3.66682072353845),
    ],
)
def test_mmff_matches_rdkit_with_platform_tolerance(mol, variant, energy, opt_energy, x0):
    assert mol.rdkit_mmff_energy(COORDS, variant) == pytest.approx(energy, abs=1e-12)
    status, e, coords = mol.rdkit_mmff_optimize(COORDS, 200, variant)
    assert status == 0
    assert e == pytest.approx(opt_energy, abs=1e-9)
    assert coords[0][0] == pytest.approx(x0, abs=1e-6)
    assert len(coords) == len(COORDS)


def test_uff_matches_rdkit_with_platform_tolerance(mol):
    assert mol.rdkit_uff_has_all_params()
    assert mol.rdkit_uff_energy(COORDS) == pytest.approx(40.395092089572906, abs=1e-12)
    status, e, coords = mol.rdkit_uff_optimize(COORDS, 200)
    assert status == 0
    assert e == pytest.approx(19.492409722746174, abs=1e-9)
    assert coords[0][0] == pytest.approx(3.764669178348734, abs=1e-6)
    labels = mol.rdkit_uff_atom_labels()
    assert labels[:5] == ["C_3", "C_R", "O_R", "N_R", "C_R"]


@pytest.mark.parametrize("variant", ["MMFF94", "MMFF94s"])
def test_persistent_mmff_force_field_reuses_prepared_terms(mol, variant):
    field = mol.mmff_force_field(COORDS, variant)
    assert isinstance(field, chematic.MolecularForceField)
    assert field.kind == variant
    assert field.num_atoms == len(COORDS)
    assert field.positions == COORDS
    assert field.energy() == pytest.approx(mol.rdkit_mmff_energy(COORDS, variant), abs=1e-12)

    energy, gradient = field.energy_gradient()
    assert energy == pytest.approx(field.energy(), abs=1e-12)
    assert len(gradient) == len(COORDS)
    assert all(len(point) == 3 and all(math.isfinite(x) for x in point) for point in gradient)

    moved = list(COORDS[0])
    moved[0] += 0.01
    assert field.set_position_(0, moved) is field
    assert field.position(0) == moved
    changed_energy = field.energy()
    assert changed_energy != energy
    assert field.set_positions_(COORDS) is field
    assert field.energy() == pytest.approx(energy, abs=1e-12)

    status, optimized_energy = field.minimize_()
    assert status == 0
    assert optimized_energy < energy
    assert len(field.positions) == len(COORDS)


def test_consistent_torsion_gradient_available_for_mmff(mol):
    mmff = mol.mmff_force_field(COORDS)
    compatible = mmff.gradient()
    consistent = mmff.gradient_with_consistent_torsions()
    consistent_planes = mmff.gradient_with_consistent_planes()
    assert len(consistent) == len(compatible)
    for a, b in zip(compatible, consistent):
        assert a == pytest.approx(b, abs=1e-8)
    for a, b in zip(compatible, consistent_planes):
        assert a == pytest.approx(b, abs=1e-8)
    uff = mol.uff_force_field(COORDS)
    with pytest.raises(ValueError, match="require an MMFF"):
        uff.gradient_with_consistent_torsions()
    with pytest.raises(ValueError, match="require an MMFF"):
        uff.gradient_with_consistent_planes()


def test_consistent_plane_gradient_near_linear_amine():
    fixture_path = Path(__file__).resolve().parents[2] / "chematic-ff/tests/mmff94_near_linear_amine.json"
    fixture = json.loads(fixture_path.read_text())
    mol = chematic.from_smiles(fixture["smiles"]).add_hydrogens()
    field = mol.mmff_force_field(fixture["coords"])
    compatible = field.gradient()[5][1]
    consistent = field.gradient_with_consistent_planes()[5][1]
    coords = [point[:] for point in fixture["coords"]]
    step = 1e-8
    coords[5][1] += step
    plus = field.set_positions_(coords).energy()
    coords[5][1] -= 2 * step
    minus = field.set_positions_(coords).energy()
    numeric = (plus - minus) / (2 * step)
    assert abs(consistent - numeric) < 0.01 * abs(numeric)
    assert abs(compatible - numeric) > 1000


def test_persistent_uff_force_field_and_input_validation(mol):
    field = mol.uff_force_field(COORDS)
    assert field.kind == "UFF"
    assert field.energy() == pytest.approx(mol.rdkit_uff_energy(COORDS), abs=1e-12)
    assert len(field.gradient()) == len(COORDS)

    with pytest.raises(IndexError):
        field.position(len(COORDS))
    with pytest.raises(ValueError, match="finite"):
        field.set_position_(0, [math.nan, 0.0, 0.0])
    with pytest.raises(ValueError, match="expected"):
        field.set_positions_(COORDS[:-1])
    with pytest.raises(ValueError, match="non-negative"):
        mol.uff_force_field(COORDS, vdw_thresh=-1.0)
    with pytest.raises(ValueError, match="MMFF94"):
        mol.mmff_force_field(COORDS, variant="MMFF95")


def test_uff_missing_params_and_bad_input(mol):
    # RDKit 2026.03.1 UFFHasAllMoleculeParams after AddHs.
    for smi, expected in [("C=[Se]", False), ("[Cu]C", False), ("C[Se]C", True), ("C[Mg]Cl", True)]:
        assert chematic.from_smiles(smi).add_hydrogens().rdkit_uff_has_all_params() is expected
    with pytest.raises(ValueError):
        mol.rdkit_uff_energy(COORDS[:-1])
    with pytest.raises(ValueError):
        mol.rdkit_mmff_optimize(COORDS, 200, "MMFF95")


def test_heavy_conformer_can_be_ranked_after_explicit_hydrogen_refinement():
    mol = chematic.from_smiles("CCO")
    coords = mol.generate_3d()
    refined, all_atom_coords, energy, converged, iterations, residual = (
        mol.refine_mmff94_with_explicit_hydrogens(coords)
    )
    assert len(refined) == len(coords) == 3
    assert len(all_atom_coords) == len(mol.add_hydrogens().atom_table)
    assert math.isfinite(energy)
    assert math.isfinite(residual) and residual >= 0.0
    assert isinstance(converged, bool)
    assert 0 <= iterations <= 200

    with pytest.raises(ValueError, match="expected 3 atom coordinates"):
        mol.refine_mmff94_with_explicit_hydrogens(coords[:-1])
