"""Public typed consumer: mypy strict must accept this file."""
from typing import List, Tuple

from chematic import (
    BondFragmentation,
    FragmentAttachment,
    Mol,
    fragment_on_bonds_checked,
    from_smiles,
)

source: Mol = from_smiles("CCO")
counts: List[int] = [0, 0, 0]
result: BondFragmentation = fragment_on_bonds_checked(
    source, (0,), dummy_labels=[(7, 9)], cuts_per_atom=counts
)
atoms: Tuple[int, ...] = result.cuts_per_atom
attachment: FragmentAttachment = result.attachments[0]
dummy_index: int = attachment.dummy
output: Mol = result.molecule
