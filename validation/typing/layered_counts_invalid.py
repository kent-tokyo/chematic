from chematic import rdkit_compat as Chem


def wrong(mol: Chem.Mol) -> None:
    Chem.LayeredFingerprint(mol, atomCounts=(0, 0, 0))
    Chem.LayeredFingerprint(mol, fpSize=1.5)
    Chem.LayeredFingerprint(mol, setOnlyBits=b"mask")
    Chem.LayeredFingerprint(mol, fromAtoms=iter([0]))
