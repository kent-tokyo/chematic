"""Morgan fingerprint + Tanimoto similarity ranking using chematic.rdkit_compat.

RDKit equivalent:
    from rdkit.Chem import rdMolDescriptors, DataStructs
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    DataStructs.BulkTanimotoSimilarity(query_fp, fps)
"""
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import rdMolDescriptors, DataStructs

QUERY = "CC(=O)Oc1ccccc1C(=O)O"          # aspirin
LIBRARY = {
    "salicylic acid": "O=C(O)c1ccccc1O",
    "ibuprofen":      "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "benzene":        "c1ccccc1",
    "paracetamol":    "CC(=O)Nc1ccc(O)cc1",
    "aspirin-copy":   "CC(=O)Oc1ccccc1C(=O)O",
}


def fp(smi):
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles(smi), 2, nBits=2048
    )


def main():
    qfp = fp(QUERY)
    names = list(LIBRARY)
    sims = DataStructs.BulkTanimotoSimilarity(qfp, [fp(LIBRARY[n]) for n in names])

    ranked = sorted(zip(names, sims), key=lambda x: x[1], reverse=True)
    for name, sim in ranked:
        print(f"{sim:.3f}  {name}")

    # The exact-copy must rank first with similarity 1.0.
    assert ranked[0][0] == "aspirin-copy"
    assert ranked[0][1] == 1.0
    print("OK: aspirin-copy ranked first (Tanimoto = 1.0)")


if __name__ == "__main__":
    main()
