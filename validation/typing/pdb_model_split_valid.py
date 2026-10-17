from pathlib import Path
from typing import Optional, Tuple

from chematic import PdbModelFile, PdbModelSplit, PdbSplitLoss, split_pdb_models


def consume(text: str) -> Tuple[str, Optional[int], int]:
    result: PdbModelSplit = split_pdb_models(text, prefix="nmr", max_models=10)
    files: Tuple[PdbModelFile, ...] = result.files
    losses: Tuple[PdbSplitLoss, ...] = result.losses
    source: str = result.source_text
    paths = [Path(row.filename) for row in files]
    assert source and paths and len(losses) >= 0
    return files[0].pdb, files[0].model_num, files[0].atom_count
