"""Exactly six expected strict mypy errors."""
from chematic import Lattice, PeriodicStructure, Site, write_crystal_stl_checked

source = PeriodicStructure(Lattice.cubic(1.), [Site([('C', 1.)], (0., 0., 0.))])
write_crystal_stl_checked(source, [.5])
write_crystal_stl_checked('CIF', [.5], millimeters_per_angstrom=3.)
write_crystal_stl_checked(source, ['large'], millimeters_per_angstrom=3.)
write_crystal_stl_checked(source, [.5], millimeters_per_angstrom='three')
write_crystal_stl_checked(source, [.5], millimeters_per_angstrom=3., subdivisions='two')
result = write_crystal_stl_checked(source, [.5], millimeters_per_angstrom=3.)
result.triangle_count = 1
