from typing import Any
import numpy as np
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory
ranker = rdInfoTheory.InfoBitRanker(8, 2, rdInfoTheory.InfoType.CHISQUARE)
vector = ExplicitBitVect(8)
ranker.AccumulateVotes(vector, 0)
rows: np.ndarray[Any, np.dtype[np.float64]] = ranker.GetTopN(8)
