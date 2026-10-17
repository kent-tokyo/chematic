"""Exactly six expected strict mypy errors."""
from chematic import from_smiles, replace_core_checked

source, core = from_smiles("CCOC"), from_smiles("CO")
replace_core_checked("CCOC", core, [1, 2])
replace_core_checked(source, None, [1, 2])
replace_core_checked(source, core, ["one", 2])
replace_core_checked(source, core, iter([1, 2]))
result = replace_core_checked(source, core, [1, 2])
result.original_to_output = ()
result.attachments[0].dummy = 99
