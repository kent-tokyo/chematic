"""Four deliberate invalid types, expected mypy exit 1."""

from chematic import rdkit_compat as C

mol = C.MolFromSmiles("CC")
assert mol is not None
C.GetAdjacencyMatrix("CC")
C.GetAdjacencyMatrix(mol, emptyVal=0.5)
C.GetDistanceMatrix(mol, prefix=b"")
C.GetDistanceMatrix(mol, useAtomWts="true")
