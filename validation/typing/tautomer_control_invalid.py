import chematic as ch


def wrong_callback(progress: ch.TautomerProgress) -> int:
    return progress.transforms


mol = ch.from_smiles("CC=O")
ch.enumerate_tautomers_checked(mol, callback=wrong_callback)
ch.enumerate_tautomers_checked("CC=O")
ch.enumerate_tautomers_checked(mol, max_tautomers="10")
result = ch.enumerate_tautomers_checked(mol)
wrong_result: str = result.smiles
result.status = "Canceled"
