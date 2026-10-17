from typing import Tuple
from chematic import PdbHydrogenName, PdbHydrogenNaming, name_pdb_hydrogens_checked
text: str = ""
report: PdbHydrogenNaming = name_pdb_hydrogens_checked(text, [1, 3])
assignments: Tuple[PdbHydrogenName, ...] = report.assignments
pdb: str = report.pdb
serial: int = assignments[0].serial
