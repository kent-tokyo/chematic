import numpy as np
from chematic.rdkit_compat import DataStructs, ExplicitBitVect

fp = ExplicitBitVect(8)
array = np.zeros(8, dtype=np.int8)
DataStructs.ConvertToNumpyArray(fp, array)
DataStructs.ConvertToNumpyArray(fp)


class SparseSource:
    def GetNumBits(self) -> int:
        return 8

    def GetOnBits(self) -> list[int]:
        return [2]


DataStructs.ConvertToNumpyArray(SparseSource(), array)
