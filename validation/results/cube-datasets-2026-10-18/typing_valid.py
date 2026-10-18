from chematic import CubeDatasetDocument, CubeDatasetSelection, VolumetricGrid, read_cube_datasets
import numpy as np

def select(text: str) -> VolumetricGrid:
    doc: CubeDatasetDocument = read_cube_datasets(text, max_atoms=20)
    selected: CubeDatasetSelection = doc.select(np.int64(0))
    source: str = selected.source_text
    values: tuple[tuple[float, ...], ...] = doc.datasets
    return selected.to_grid()
