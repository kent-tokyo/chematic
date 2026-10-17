from dataclasses import dataclass
from typing import Any
import numpy as np
from chematic.rdkit_compat import rdAlignment


@dataclass(frozen=True)
class Point:
    x: float
    y: float
    z: float


rows = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
ssd, transform = rdAlignment.GetAlignmentTransform(rows, rows, weights=[1.0, 2.0])
value: float = ssd
matrix: np.ndarray[Any, np.dtype[np.float64]] = transform
rdAlignment.GetAlignmentTransform(np.array(rows), np.array(rows), reflect=True)
rdAlignment.GetAlignmentTransform([Point(0.0, 0.0, 0.0)], [Point(1.0, 2.0, 3.0)])
