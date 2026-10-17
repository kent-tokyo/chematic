"""Typed public XYZ atom data consumer."""

from chematic import XyzAtomDataFrame, parse_xyz_atom_data

frame = parse_xyz_atom_data("1\n\nC 0 0 0 1 2 3\n")
value: float = frame.values[0][0]
text: str = frame.to_extxyz("velocities")
scalar = XyzAtomDataFrame(b"1\n\nH 0 0 0 .25\n", components=1)
