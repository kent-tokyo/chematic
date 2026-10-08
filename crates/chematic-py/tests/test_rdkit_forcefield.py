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


def test_uff_missing_params_and_bad_input(mol):
    # RDKit 2026.03.1 UFFHasAllMoleculeParams after AddHs.
    for smi, expected in [("C=[Se]", False), ("[Cu]C", False), ("C[Se]C", True), ("C[Mg]Cl", True)]:
        assert chematic.from_smiles(smi).add_hydrogens().rdkit_uff_has_all_params() is expected
    with pytest.raises(ValueError):
        mol.rdkit_uff_energy(COORDS[:-1])
    with pytest.raises(ValueError):
        mol.rdkit_mmff_optimize(COORDS, 200, "MMFF95")
