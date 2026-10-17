"""Four intentional type errors; do not execute."""

from chematic import parse_xyz_atom_data

parse_xyz_atom_data(123)
parse_xyz_atom_data("text", components=2)
frame = parse_xyz_atom_data("1\n\nC 0 0 0 1 2 3\n")
frame.to_extxyz(1)
frame.comment = "new"
