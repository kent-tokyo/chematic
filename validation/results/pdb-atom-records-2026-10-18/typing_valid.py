"""Public typed atom attributes, including absent numeric fields."""

from typing import Optional
from chematic import PdbAtomRead, PdbAtomRecord, read_pdb_atom_records

result: PdbAtomRead = read_pdb_atom_records("")
if result.atoms:
    atom: PdbAtomRecord = result.atoms[0]
    occupancy: Optional[float] = atom.occupancy
    position: tuple[float, float, float] = atom.position_angstrom
    original: str = atom.record_text
