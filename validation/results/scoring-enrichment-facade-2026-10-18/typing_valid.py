"""Public typed ranked-record usage."""

from typing import Iterable
import numpy as np
from chematic.rdkit_compat import RankedScoreRow, Scoring

rows: Iterable[RankedScoreRow] = (
    (score, label) for score, label in [(1.0, 1), (0.0, 0)]
)
result: list[float] = Scoring.CalcEnrichment(rows, np.int64(1), [0.5, 1.0])
