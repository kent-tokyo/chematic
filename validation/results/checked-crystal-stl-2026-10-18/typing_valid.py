from typing import Literal
from chematic import Lattice, PeriodicStructure, Site, StlMeshExport, write_crystal_stl_checked

source = PeriodicStructure(Lattice.cubic(1.), [Site([('C', 1.)], (0., 0., 0.))])
result: StlMeshExport = write_crystal_stl_checked(source, [.5], millimeters_per_angstrom=3.)
text: str = result.text
count: int = result.triangle_count
units: Literal['millimeter'] = result.units
scale: float = result.millimeters_per_angstrom
