from typing import List
import numpy as np
from chematic import align_conformers_checked

rms: List[float] = [2.0]
coordinates = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [[2.0, 3.0, 0.0], [2.0, 4.0, 0.0]]]
result = align_conformers_checked(
    coordinates, atom_ids=[0, 1], conformer_ids=[1, 0], weights=(1.0, 2.0), rms_list=rms
)
result2 = align_conformers_checked(np.zeros((2, 5, 3)), weights=np.ones(5))
