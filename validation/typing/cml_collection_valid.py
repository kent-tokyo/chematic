from typing import List, Optional
from chematic import (
    CmlMoleculeRecord,
    from_smiles,
    read_cml_collection,
    write_cml_collection,
)

output: str = write_cml_collection(
    [from_smiles("CCO").to_cml(), from_smiles("[NH4+]").to_cml().encode()]
)
records: List[CmlMoleculeRecord] = read_cml_collection(output)
index: int = records[0].index
identifier: Optional[str] = records[0].molecule_id
source: str = records[0].source_text
molecule: str = records[0].molecule_xml
