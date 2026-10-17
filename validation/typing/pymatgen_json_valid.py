from typing import Any, Dict, List
from chematic import (
    PeriodicStructure,
    PymatgenStructureRecord,
    parse_pymatgen_structures,
)

records: List[PymatgenStructureRecord] = parse_pymatgen_structures(b'{"data": []}')
record = records[0]
source: str = record.source_text
index: int = record.record_index
losses: tuple[str, ...] = record.projection_losses
data: Dict[str, Any] = record.structure_dict()
geometry: PeriodicStructure = record.to_periodic_structure(allow_metadata_loss=True)
