"""Four deliberate invalid inputs, expected mypy exit 1."""

from pathlib import Path
from chematic import read_pdb_atom_records

read_pdb_atom_records(b"ATOM")
read_pdb_atom_records(Path("input.pdb"))
read_pdb_atom_records("", max_atoms=1.5)
read_pdb_atom_records("", max_records="2")
