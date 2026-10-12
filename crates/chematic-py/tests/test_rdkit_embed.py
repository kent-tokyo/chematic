"""Seeded ETKDG embedding and bounds checks pinned to RDKit 2026.09.1.

The complete versioned A6 lane performs corpus-scale coordinate, energy, and
bounds comparisons. These unit tests protect deterministic embedding, atom
counts, and representative RDKit 2026.09.1 bounds without making all-platform
bit-identity a unit-test contract.
"""
import pytest

import chematic

EMBED_42 = {
    'CC(=O)Nc1ccc(O)cc1': [
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
    ],
    'C[C@H](N)C(=O)O': [
        [-0.8592407634599415, 1.1514661472637542, -0.16522593775029595],
        [-0.1853754145711437, -0.15987068859904846, -0.5099590208150512],
        [-0.802105522995388, -1.1734690014909304, 0.28318342058733675],
        [1.2541597265949056, -0.13057261963056138, -0.1652195934721842],
        [1.972011877866091, 0.7167166775677084, -0.7762230953089931],
        [1.821236640219406, -0.9882817051797496, 0.7878718150743558],
        [-0.2750021524916736, 1.5837702964536673, 0.6895632730001544],
        [-0.9005780761432235, 1.8285273117575287, -1.0365760821155205],
        [-1.8862461561959238, 0.996256849997651, 0.2550251212858728],
        [-0.33699872662889885, -0.34379839275665236, -1.583884215673667],
        [-1.3606582832409595, -0.8755071898300736, 1.0831481019512659],
        [-1.0289647171341654, -2.0490815471739134, -0.18946831775330764],
        [2.587761568180817, -0.5561561383794941, 1.3277645309899992],
    ],
    'CC(=O)[O-].[Na+]': [
        [-0.649681113578404, -0.05323018228431236, -0.12806740711050485],
        [0.8166224987437282, 0.022958662150776684, 0.15995277160566965],
        [1.4411856388140436, -0.9629086179421968, 0.6350507380783897],
        [1.509380961298548, 1.2011989416423827, -0.09927670436629632],
        [0.0, 0.0, 0.0],
        [-0.9464615638849146, -1.0531510530947263, -0.5225972149389126],
        [-0.9259870457122426, 0.7539736605594457, -0.8416937017286762],
        [-1.2450593756807549, 0.09115858896863063, 0.7966315184603254],
    ],
}

BOUNDS_CCO = [
    [0.0, 1.524, 2.4156792843222545, 1.119400794877744, 1.119400794877744, 1.119400794877744, 2.1950665941551053, 2.1950665941551053, 3.2472182203294455],
    [1.504, 0.0, 1.4038448452526835, 2.1950665941551053, 2.1950665941551053, 2.1950665941551053, 1.119400794877744, 1.119400794877744, 2.0008385531203197],
    [2.3356792843222545, 1.3838448452526835, 0.0, 3.3299599465639718, 3.3299599465639718, 3.3299599465639718, 2.09083554180725, 2.09083554180725, 1.0002542040013616],
    [1.099400794877744, 2.1150665941551052, 2.3648500690570953, 0.0, 1.8519655808531272, 1.8519655808531272, 3.0753810004351188, 3.0753810004351188, 4.195905147275425],
    [1.099400794877744, 2.1150665941551052, 2.3648500690570953, 1.7719655808531272, 0.0, 1.8519655808531272, 3.0753810004351188, 3.0753810004351188, 4.195905147275425],
    [1.099400794877744, 2.1150665941551052, 2.3648500690570953, 1.7719655808531272, 1.7719655808531272, 0.0, 3.0753810004351188, 3.0753810004351188, 4.195905147275425],
    [2.1150665941551052, 1.099400794877744, 2.01083554180725, 2.2546511899391777, 2.2546511899391777, 2.2546511899391777, 0.0, 1.8519655808531272, 2.881872377536109],
    [2.1150665941551052, 1.099400794877744, 2.01083554180725, 2.2546511899391777, 2.2546511899391777, 2.2546511899391777, 1.7719655808531272, 0.0, 2.881872377536109],
    [2.2837846345056976, 1.9208385531203198, 0.9802542040013615, 1.68, 1.68, 1.68, 2.097732847201809, 2.097732847201809, 0.0],
]


@pytest.mark.parametrize("smiles", sorted(EMBED_42))
def test_rdkit_embed_matches_rdkit_seed_42(smiles):
    got = chematic.from_smiles(smiles).add_hydrogens().rdkit_embed(random_seed=42)
    repeated = chematic.from_smiles(smiles).add_hydrogens().rdkit_embed(random_seed=42)
    # EMBED_42 retains the historical fixture solely as a compact atom-count
    # oracle. RDKit changed ETKDG bounds after 2026.03, so corpus-scale numeric
    # parity belongs to the versioned A6 artifact lane.
    assert len(got) == len(EMBED_42[smiles])
    assert got == repeated


def test_rdkit_initial_metric_embedding_matches_rdkit_2026_09():
    mol = chematic.from_smiles("CCO").add_hydrogens()
    got = mol.rdkit_embed_initial(random_seed=42)
    expected = [
        [-0.7500043088651186, 0.4084390159458829, 0.06848496479484886],
        [0.44564780848200714, -0.30357044432024505, -0.2925815152031927],
        [1.1754601191726461, -0.458924916623023, 1.1664313728171822],
        [-1.2354609528915097, 1.0831648725953975, -0.5316473904156442],
        [-1.5724182037685812, -1.0383486994621525, 0.45351405445069587],
        [-1.2999968573900929, 0.40241491422299236, 0.40181982062149446],
        [0.5525517516037974, -0.5577005629983048, -0.9232674186573715],
        [0.9120715852990083, -0.80779822184343, -0.5889154542504689],
        [1.7721490583578425, 1.2723240424828839, 0.24616156584245455],
    ]
    for actual_atom, expected_atom in zip(got, expected):
        assert actual_atom == pytest.approx(expected_atom, abs=5e-4)
    assert got == mol.rdkit_embed_initial(random_seed=42)


def test_rdkit_embed_multiple_confs_matches_rdkit_seed_progression():
    mol = chematic.from_smiles("CCO").add_hydrogens()
    conformers = mol.rdkit_embed_multiple_confs(3, random_seed=42)
    assert len(conformers) == 3
    assert conformers[0] == mol.rdkit_embed(random_seed=42)
    expected_first_atoms = [
        [-0.9639577069927618, 0.015102231353043739, 0.0410206427653118],
        [-0.946941193051227, 0.20466995457287168, 0.007322645786689314],
        [-0.8949265887204274, 0.21063088448005707, -0.017107180966725415],
    ]
    for conformer, expected in zip(conformers, expected_first_atoms):
        assert conformer[0] == pytest.approx(expected, abs=1e-4)
    assert mol.rdkit_embed_multiple_confs(0) == []


def test_rdkit_embed_multiple_confs_mmff_matches_individual_optimizer():
    mol = chematic.from_smiles("CCO").add_hydrogens()
    embedded = mol.rdkit_embed_multiple_confs(3, random_seed=42)
    optimized = mol.rdkit_embed_multiple_confs_mmff(3, random_seed=42)
    repeated = mol.rdkit_embed_multiple_confs_mmff(3, random_seed=42)
    assert optimized == repeated
    assert len(optimized) == 3
    status, energy, coords = mol.rdkit_mmff_optimize(embedded[0], max_iters=200)
    assert optimized[0][0] == status
    assert optimized[0][1] == energy
    assert optimized[0][2] == coords


def test_rdkit_embed_multiple_confs_mmff_diagnostics_preserve_results():
    mol = chematic.from_smiles("CCCO").add_hydrogens()
    compact = mol.rdkit_embed_multiple_confs_mmff(4, max_mmff_iters=2)
    detailed = mol.rdkit_embed_multiple_confs_mmff_diagnostics(4, max_mmff_iters=2)
    assert len(compact) == len(detailed)
    for compact_row, detailed_row in zip(compact, detailed):
        status, energy, iterations, gradient_test, coords = detailed_row
        assert compact_row == (status, energy, coords)
        assert iterations == 2
        assert gradient_test is not None and gradient_test >= 0.0


def test_rdkit_embed_multiple_confs_mmff_adaptive_retry_is_deterministic():
    mol = chematic.from_smiles("CCCO").add_hydrogens()
    initial = mol.rdkit_embed_multiple_confs_mmff(4, max_mmff_iters=1)
    adaptive = mol.rdkit_embed_multiple_confs_mmff(
        4,
        max_mmff_iters=1,
        retry_mmff_iters=200,
        retry_candidates=1,
    )
    repeated = mol.rdkit_embed_multiple_confs_mmff(
        4,
        max_mmff_iters=1,
        retry_mmff_iters=200,
        retry_candidates=1,
    )
    assert adaptive == repeated
    assert sum(status == 0 for status, _, _ in initial) == 0
    assert sum(status == 0 for status, _, _ in adaptive) == 1
    assert sum(a[2] != b[2] for a, b in zip(adaptive, initial)) == 1
    assert all(a[1] <= b[1] for a, b in zip(adaptive, initial))


def test_rdkit_embed_multiple_confs_mmff_rescue_skips_converged_ensemble():
    mol = chematic.from_smiles("CCCO").add_hydrogens()
    initial = mol.rdkit_embed_multiple_confs_mmff(4, max_mmff_iters=200)
    rescued = mol.rdkit_embed_multiple_confs_mmff(
        4,
        max_mmff_iters=200,
        retry_mmff_iters=1_000,
        retry_candidates=1,
        retry_only_if_none_converged=True,
    )
    assert any(status == 0 for status, _, _ in initial)
    assert rescued == initial


def test_rdkit_embed_multiple_confs_mmff_gradient_rescue_is_bounded():
    mol = chematic.from_smiles("CCCO").add_hydrogens()
    initial = mol.rdkit_embed_multiple_confs_mmff(4, max_mmff_iters=1)
    rescued = mol.rdkit_embed_multiple_confs_mmff_gradient_rescue(
        4,
        max_mmff_iters=1,
        retry_mmff_iters=200,
        retry_candidates=1,
        minimum_gradient_test=0.0,
    )
    assert sum(status == 0 for status, _, _ in rescued) == 1
    assert sum(a[2] != b[2] for a, b in zip(rescued, initial)) == 1
    assert all(a[1] <= b[1] for a, b in zip(rescued, initial))

    with pytest.raises(ValueError, match="minimum_gradient_test"):
        mol.rdkit_embed_multiple_confs_mmff_gradient_rescue(
            1, minimum_gradient_test=float("nan")
        )


def test_rdkit_bounds_matrix_matches_rdkit():
    got = chematic.from_smiles("CCO").add_hydrogens().rdkit_bounds_matrix()
    assert len(got) == 9
    expected_cells = {
        (0, 1): 1.524,
        (1, 0): 1.504,
        (0, 2): 2.421147938383115,
        (2, 0): 2.3298908238196163,
        (0, 8): 3.3076831787357386,
        (8, 0): 2.1884247952025846,
        (3, 8): 4.197112901631307,
        (8, 3): 1.68,
    }
    for (row, column), expected in expected_cells.items():
        assert got[row][column] == pytest.approx(expected, abs=1e-12)
