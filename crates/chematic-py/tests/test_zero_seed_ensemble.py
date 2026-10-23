"""RDKit #9738: zero-seed ensembles preserve reproducibility and vary shapes."""
import itertools
import math

import pytest
import chematic


@pytest.mark.parametrize("smiles", ["CCCCCCO", "CCCCCC", "NCCCCCC", "CC(C)CCCO"])
@pytest.mark.parametrize("seed", [0, 1, 42])
def test_ensemble_seed_schedule_and_internal_geometry(smiles, seed):
    mol = chematic.from_smiles(smiles).add_hydrogens()
    before = mol.to_moljson()
    ensemble = mol.rdkit_embed_multiple_confs(5, random_seed=seed)
    assert len(ensemble) == 5
    assert ensemble == mol.rdkit_embed_multiple_confs(5, random_seed=seed)
    assert ensemble[:3] == mol.rdkit_embed_multiple_confs(3, random_seed=seed)
    seeds = [0, 2, 3, 4, 5] if seed == 0 else [seed * i for i in range(1, 6)]
    for coordinates, attempt_seed in zip(ensemble, seeds):
        assert coordinates == mol.rdkit_embed(random_seed=attempt_seed)
        assert all(math.isfinite(x) for row in coordinates for x in row)
    heavy = chematic.from_smiles(smiles).atom_count
    pairs = list(itertools.combinations(range(heavy), 2))
    # Internal distances remove rigid translation/rotation as an explanation.
    differences = [
        max(abs(math.dist(a[i], a[j]) - math.dist(b[i], b[j])) for i, j in pairs)
        for a, b in itertools.combinations(ensemble, 2)
    ]
    assert max(differences) > 0.1
    assert mol.to_moljson() == before


def test_zero_seed_mmff_uses_each_embedded_candidate():
    mol = chematic.from_smiles("CCCCCCO").add_hydrogens()
    embedded = mol.rdkit_embed_multiple_confs(3, random_seed=0)
    optimized = mol.rdkit_embed_multiple_confs_mmff(3, random_seed=0, max_mmff_iters=0)
    assert len(optimized) == 3
    for coordinates, result in zip(embedded, optimized):
        assert result == mol.rdkit_mmff_optimize(coordinates, max_iters=0)


def test_zero_seed_empty_and_invalid_ranges():
    mol = chematic.from_smiles("CCO").add_hydrogens()
    assert mol.rdkit_embed_multiple_confs(0, random_seed=0) == []
    with pytest.raises(ValueError):
        mol.rdkit_embed_multiple_confs(2, random_seed=2**31 - 1)
    with pytest.raises(ValueError):
        mol.rdkit_embed_multiple_confs(1, random_seed=-1)
