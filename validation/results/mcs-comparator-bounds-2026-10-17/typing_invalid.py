"""Strict negative consumer: exactly six expected argument-type errors."""
from chematic import McsAtomCompare, McsBondCompare, find_mcs, find_mcs_checked

find_mcs([], atom_compare="typo")
find_mcs([], bond_compare="typo")
find_mcs_checked([], atom_compare=McsBondCompare.ANY)
find_mcs_checked([], bond_compare=McsAtomCompare.ANY)
find_mcs([], atom_compare=lambda: True)
find_mcs_checked([], bond_compare=123)
