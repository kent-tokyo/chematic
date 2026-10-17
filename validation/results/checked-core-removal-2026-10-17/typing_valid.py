"""Strict public consumer for explicit typed core matches and frozen metadata."""
from typing import Optional, Tuple
from chematic import (
    CoreAttachment, CoreReplacement, Mol, from_smiles, replace_core_checked,
)

source: Mol = from_smiles("CCOC")
core: Mol = from_smiles("CO")
result: CoreReplacement = replace_core_checked(source, core, (1, 2))
mapping: Tuple[Optional[int], ...] = result.original_to_output
attachment: CoreAttachment = result.attachments[0]
anchor: int = attachment.anchor
output: Mol = result.molecule
