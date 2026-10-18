from typing import Iterable, Mapping
from chematic import PdbMergeResult, PdbMergeAtom, merge_pdb_coordinate_blocks

def merge(sources: Iterable[str], maps: Iterable[Mapping[str, str]]) -> str:
    result: PdbMergeResult = merge_pdb_coordinate_blocks(sources, chain_maps=maps)
    atoms: tuple[PdbMergeAtom, ...] = result.atoms
    originals: tuple[str, ...] = result.source_texts
    return result.pdb
