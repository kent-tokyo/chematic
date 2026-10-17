"""Four intentional public type errors; do not execute."""

from chematic import Lattice, Site, PeriodicStructure, write_xsf_checked

structure = PeriodicStructure(Lattice.cubic(4.0), [Site([("C", 1.0)], (0.0, 0.0, 0.0))])
write_xsf_checked("CIF", allow_metadata_loss=True)
write_xsf_checked(structure, allow_metadata_loss="yes")
PeriodicStructure.from_cif("CIF", expand_symmetry="yes")
status = structure.symmetry_status
if status is not None:
    status.is_complete_cell = True
