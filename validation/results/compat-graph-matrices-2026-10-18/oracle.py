"""Live RDKit comparisons for the bounded compatibility NumPy matrices."""

import hashlib
import json
from pathlib import Path
import random

import numpy as np
from rdkit import Chem as R, rdBase
from chematic import rdkit_compat as C

SMILES = [
    "",
    "C",
    "O",
    "N",
    "C=C",
    "C#N",
    "CC(=O)O",
    "c1ccccc1",
    "n1ccccc1",
    "CC.O",
    "C1CCC1",
    "C1CCCCC1",
    "N->[Cu]",
    "[Cu]<-N",
    "[Mo]$[Mo]",
    "[Na+].[Cl-]",
    "[2H]O[2H]",
    "C/C=C/F",
    "F[C@](Cl)(Br)I",
    "O=C=O",
    "CCN(C)C",
    "C1CC2CCC1C2",
    "[NH4+]",
    "C1=CC=CC=C1",
    "COC",
    "[CH3]",
    "[13CH4]",
    "[nH]1cccc1",
    "CC(=O)[O-]",
    "C1CC1.O.N",
    "CS(=O)(=O)C",
    "P(=O)(O)(O)O",
]
ADDED_H = {
    "C",
    "O",
    "N",
    "C=C",
    "C#N",
    "CC(=O)O",
    "c1ccccc1",
    "C1CCC1",
    "C1CCCCC1",
    "COC",
    "[13CH4]",
    "[nH]1cccc1",
}
HERE = Path(__file__).parent
rng = random.Random(9647)
matrices = variants = cells = 0
maximum_error = 0.0
inputs_hash = hashlib.sha256()

for smiles in SMILES:
    initial = C.MolFromSmiles(smiles) if smiles else C.RWMol().GetMol()
    reference = R.MolFromSmiles(smiles)
    assert initial is not None and reference is not None
    for add_h in [False, True] if smiles in ADDED_H else [False]:
        native = C.Mol(initial._mol.rdkit_add_hydrogens()) if add_h else initial
        upstream = R.AddHs(reference) if add_h else R.Mol(reference)
        n = upstream.GetNumAtoms()
        assert native.GetNumAtoms() == n
        order = list(range(n))
        shuffled = order.copy()
        rng.shuffle(shuffled)
        orders = [
            list(o)
            for o in dict.fromkeys(tuple(o) for o in [order, order[::-1], shuffled])
        ]
        for permutation in orders:
            mol = C.RenumberAtoms(native, permutation)
            ref = R.RenumberAtoms(upstream, permutation) if n else R.Mol(upstream)
            variants += 1
            inputs_hash.update(json.dumps([smiles, add_h, permutation]).encode())
            for weighted in [False, True]:
                for force in [False, True]:
                    want = R.GetAdjacencyMatrix(ref, useBO=weighted, force=True)
                    got = C.GetAdjacencyMatrix(mol, useBO=weighted, force=force)
                    assert got.shape == want.shape == (n, n)
                    assert got.dtype == want.dtype
                    assert got.flags.owndata and got.flags.c_contiguous
                    np.testing.assert_array_equal(got, want)
                    matrices += 1
                    cells += n * n
                for atom_weights in [False, True]:
                    for force in [False, True]:
                        want = R.GetDistanceMatrix(
                            ref, useBO=weighted, useAtomWts=atom_weights, force=True
                        )
                        got = C.GetDistanceMatrix(
                            mol, useBO=weighted, useAtomWts=atom_weights, force=force
                        )
                        assert got.shape == want.shape == (n, n)
                        assert got.dtype == want.dtype == np.float64
                        assert got.flags.owndata and got.flags.c_contiguous
                        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)
                        if n:
                            maximum_error = max(
                                maximum_error, float(np.max(np.abs(got - want)))
                            )
                        matrices += 1
                        cells += n * n
summary = dict(
    rdkit=rdBase.rdkitVersion,
    numpy=np.__version__,
    original_structures=len(SMILES),
    added_h_profiles=len(ADDED_H),
    variants=variants,
    matrices=matrices,
    compared_cells=cells,
    max_distance_difference=maximum_error,
    adjacency="exact values and dtype",
    inputs_sha256=inputs_hash.hexdigest(),
    scope="Bounded supported graphs, model indices unchanged; no large unsafe RDKit matrix calls",
)
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
