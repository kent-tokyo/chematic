from typing import Tuple
import chematic as ch

source = ch.PeriodicStructure.from_cif("example", expand_symmetry=True)
status = source.symmetry_status
if status is not None:
    complete: bool = status.is_complete_cell
result: ch.ObjMeshExport = ch.write_crystal_obj_checked(source, [.5], subdivisions=2)
text: str = result.text
count: int = result.face_count
omitted: Tuple[str, ...] = result.omitted_features
