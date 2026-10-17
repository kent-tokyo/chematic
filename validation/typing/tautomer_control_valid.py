from typing import Tuple
import chematic as ch


def continue_until_three(progress: ch.TautomerProgress) -> bool:
    count: int = progress.transforms
    keys: Tuple[str, ...] = progress.keys
    return count < 100 and len(keys) < 3


result: ch.TautomerEnumeration = ch.enumerate_tautomers_checked(
    ch.from_smiles("CC(=O)CC=O"), callback=continue_until_three
)
smiles: Tuple[str, ...] = result.smiles
