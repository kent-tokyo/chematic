"""Record known existing writer limitations; not a serialization parity gate."""
import json
from pathlib import Path

import chematic as ch
from rdkit import Chem

result = ch.fragment_on_bonds_checked(ch.from_smiles("CCOC"), [0])
smiles = Chem.MolFromSmiles(result.molecule.smiles)
block = result.molecule.to_mol_block()
v2000 = Chem.MolFromMolBlock(block, removeHs=False)
moljson = json.loads(result.molecule.to_moljson())
record = {
    "native_smiles": result.molecule.smiles,
    "expected_canonical": Chem.MolToSmiles(smiles),
    "v2000_canonical": Chem.MolToSmiles(v2000) if v2000 is not None else None,
    "v2000_block": block,
    "moljson_dummy_elements": [a["element"] for a in moljson["atoms"][4:]],
    "scope": "Existing writers; no V2000/MolJSON dummy graph preservation claim",
}
Path(__file__).with_name("export_residuals.json").write_text(
    json.dumps(record, indent=2) + "\n"
)
print(json.dumps({k: v for k, v in record.items() if k != "v2000_block"}, indent=2))
