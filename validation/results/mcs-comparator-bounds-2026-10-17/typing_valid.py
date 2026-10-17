"""Strict positive consumer for enum modes and legacy literal strings."""
from typing import Optional
from chematic import (
    Mol, McsAtomCompare, McsAtomComparison, McsBondCompare, McsBondComparison,
    find_mcs, find_mcs_checked, from_smiles,
)

mols = [from_smiles("CCC"), from_smiles("NNN")]
atom: McsAtomComparison = McsAtomCompare.ANY_HEAVY_ATOM
bond: McsBondComparison = McsBondCompare.ORDER_OR_AROMATIC
result: Optional[Mol] = find_mcs(mols, atom_compare=atom, bond_compare=bond)
checked: tuple[Optional[Mol], bool] = find_mcs_checked(
    mols, atom_compare="any", bond_compare="any"
)
dynamic: McsAtomCompare = McsAtomCompare("any")
