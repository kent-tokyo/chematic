import chematic as C

def bad_score(mol: C.Mol) -> float:
    return 1.5

C.rdkit_pick_tautomer_scored(["CCO"], lambda m: 1)
C.rdkit_pick_tautomer_scored([C.from_smiles("C")], bad_score)
C.rdkit_pick_tautomer_scored([C.from_smiles("C")], None)
C.rdkit_pick_tautomer_scored([C.from_smiles("C")], lambda m: 1, max_candidates=2.0)
