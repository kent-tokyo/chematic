"""Public typed XSF geometry and existing symmetry status consumer."""

from typing import Optional
from chematic import Lattice, Site, PeriodicStructure, write_xsf_checked

structure = PeriodicStructure(Lattice.cubic(4.0), [Site([("C", 1.0)], (0.0, 0.0, 0.0))])
text: str = write_xsf_checked(structure, allow_metadata_loss=True)
cif_structure = PeriodicStructure.from_cif("CIF input", expand_symmetry=True)
status = cif_structure.symmetry_status
if status is not None:
    complete: bool = status.is_complete_cell
    expanded: bool = status.is_expanded
    count: Optional[int] = status.expanded_site_count
