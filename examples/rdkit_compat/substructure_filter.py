"""Filter a SMILES list by a SMARTS query using chematic.rdkit_compat.

RDKit equivalent:
    patt = Chem.MolFromSmarts("C(=O)[OH]")
    [s for s in smis if Chem.MolFromSmiles(s).HasSubstructMatch(patt)]
"""
from chematic import rdkit_compat as Chem

LIBRARY = [
    "CCO",                       # ethanol           — no acid
    "CC(=O)O",                   # acetic acid       — carboxylic acid
    "c1ccccc1",                  # benzene           — no acid
    "O=C(O)c1ccccc1",            # benzoic acid      — carboxylic acid
    "CC(=O)Oc1ccccc1C(=O)O",     # aspirin           — carboxylic acid
    "CCN",                       # ethylamine        — no acid
]


def filter_by_smarts(smis, smarts):
    patt = Chem.MolFromSmarts(smarts)
    hits = []
    for s in smis:
        mol = Chem.MolFromSmiles(s)
        if mol is not None and mol.HasSubstructMatch(patt):
            hits.append(s)
    return hits


def main():
    acids = filter_by_smarts(LIBRARY, "C(=O)[OH]")
    for s in acids:
        print("carboxylic acid:", s)

    assert "CC(=O)O" in acids
    assert "O=C(O)c1ccccc1" in acids
    assert "CCO" not in acids
    assert len(acids) == 3
    print(f"OK: {len(acids)}/{len(LIBRARY)} molecules contain a carboxylic acid")


if __name__ == "__main__":
    main()
