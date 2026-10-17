from chematic.rdkit_compat import DataStructs, ExplicitBitVect

# Expected two argument-type errors: missing source protocol, wrong destination.
DataStructs.ConvertToNumpyArray([])
DataStructs.ConvertToNumpyArray(ExplicitBitVect(8), [0] * 8)
