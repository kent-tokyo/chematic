"""Public typed usage, checked independently of legacy unrelated stubs."""

from typing import Iterable, SupportsIndex
import numpy as np
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I

ranker = I.InfoBitRanker(3, 2, I.InfoType.BIASENTROPY)
preferences: Iterable[SupportsIndex] = (np.int64(c) for c in [0])
ranker.SetBiasList(preferences)
ranker.AccumulateVotes(ExplicitBitVect(3), np.int64(0))
ranker.GetTopN(1)
