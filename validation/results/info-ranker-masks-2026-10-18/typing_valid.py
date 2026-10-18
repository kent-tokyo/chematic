from typing import Iterable, SupportsIndex
import numpy as np
from chematic.rdkit_compat import rdInfoTheory as I

class Bits:
    def GetNumBits(self) -> int:
        return 8
    def GetOnBits(self) -> Iterable[SupportsIndex]:
        return iter([1, 3])

r = I.InfoBitRanker(8, 2)
r.SetMaskBits(np.array([1, 3], dtype=np.int64))
r.SetMaskBits(bit for bit in range(4))
r.AccumulateVotes(Bits(), 1)
result: np.ndarray[tuple[int, int], np.dtype[np.float64]] = r.GetTopN(2)
