"""Negative consumer: mypy strict must report exactly these six errors."""
from chematic import fragment_on_bonds_checked, from_smiles

source = from_smiles("CCO")
fragment_on_bonds_checked("CCO", [0])
fragment_on_bonds_checked(source, ["zero"])
fragment_on_bonds_checked(source, [0], cuts_per_atom=(0, 0, 0))
fragment_on_bonds_checked(source, [0], dummy_labels=[("seven", 9)])
result = fragment_on_bonds_checked(source, [0])
result.cuts_per_atom = ()
result.attachments[0].dummy = 99
