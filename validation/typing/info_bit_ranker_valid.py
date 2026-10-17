from typing import Iterable, Any
import numpy as np
from chematic.rdkit_compat import rdInfoTheory, ExplicitBitVect


class Sparse:
    def GetNumBits(self) -> int:
        return 4

    def GetOnBits(self) -> Iterable[int]:
        return (0, 2)


ranker = rdInfoTheory.InfoBitRanker(4, 2, rdInfoTheory.InfoType.ENTROPY)
ranker.AccumulateVotes(ExplicitBitVect(4), 0)
ranker.AccumulateVotes(Sparse(), 1)
rows: np.ndarray[Any, np.dtype[np.float64]] = ranker.GetTopN(2)


class NumpySparse:
    def GetNumBits(self) -> np.int64:
        return np.int64(4)

    def GetOnBits(self) -> Iterable[np.int64]:
        return (np.int64(0), np.int64(2))


ranker.AccumulateVotes(NumpySparse(), np.int64(1))
ranker.GetTopN(np.int64(2))
rdInfoTheory.InfoBitRanker(np.int64(4), np.int64(2))
