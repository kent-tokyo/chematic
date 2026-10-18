from typing import Iterable, SupportsIndex
import numpy as np
from chematic.rdkit_compat import rdInfoTheory as I

class Bits:
    def GetNumBits(self) -> int:
        return 8
    def GetOnBits(self) -> Iterable[SupportsIndex]:
        return iter([1, 3])

generator = I.BitCorrMatGenerator()
generator.SetBitList(np.array([1, 3], dtype=np.int64))
generator.CollectVotes(Bits())
counts: np.ndarray[tuple[int], np.dtype[np.float64]] = generator.GetCorrMatrix()
