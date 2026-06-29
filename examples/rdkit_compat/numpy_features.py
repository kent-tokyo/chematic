"""Build a fingerprint feature matrix for scikit-learn / PyTorch.

RDKit equivalent:
    import numpy as np
    from rdkit.Chem import rdMolDescriptors, DataStructs
    arr = np.zeros((nBits,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
"""
import numpy as np

from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import rdMolDescriptors, DataStructs

SMILES = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O", "CCN"]
N_BITS = 1024


def feature_matrix(smis, n_bits=N_BITS):
    rows = []
    for s in smis:
        mol = Chem.MolFromSmiles(s)
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
        rows.append(DataStructs.ConvertToNumpyArray(fp))   # (n_bits,) int8
    return np.vstack(rows)


def main():
    X = feature_matrix(SMILES)
    print(f"feature matrix shape: {X.shape}  dtype: {X.dtype}")
    print(f"bits set per molecule: {X.sum(axis=1).tolist()}")

    assert X.shape == (len(SMILES), N_BITS)
    assert X.dtype == np.int8
    assert X.sum() > 0
    # e.g. ready for: RandomForestClassifier().fit(X, y)
    print("OK: (N, nBits) int8 matrix ready for sklearn / PyTorch")


if __name__ == "__main__":
    main()
