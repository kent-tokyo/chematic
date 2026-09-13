"""Minimal consumer typing surface for the published Python wheel."""

import chematic


mol = chematic.from_smiles("c1ccccc1")
assert mol.formula == "C6H6"
