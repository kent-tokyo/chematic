from typing import Iterable, SupportsIndex
import numpy as np
from chematic.rdkit_compat import MetricBitVector, rdMetricMatrixCalc as M

class Bits:
    def GetNumBits(self) -> int:
        return 8
    def GetOnBits(self) -> Iterable[SupportsIndex]:
        return iter([1, 3])

source: MetricBitVector = Bits()
distances: np.ndarray[tuple[int], np.dtype[np.float64]] = M.GetEuclideanDistMat((iter(row) for row in [[0., 0.], [3., 4.]]))
similarity: np.ndarray[tuple[int], np.dtype[np.float64]] = M.GetTanimotoSimMat([source, Bits()])
