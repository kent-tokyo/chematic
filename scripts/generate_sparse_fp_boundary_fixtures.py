"""Freeze sparse count fingerprint references using RDKit 2026.03.1."""
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator as generators, rdMolDescriptors
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
root = Path(__file__).resolve().parents[1]
sources = ['', 'C', 'CC', 'CCCC', 'C1CC1', 'C1CCC1', 'C1CCCCC1', 'CC.CN',
           'CC(=O)O', 'C#CC#N', 'C=CC=C', 'NCCO', 'CS(=O)(=O)C', 'CP(=O)(O)O',
           'c1ccccc1', 'c1ccncc1', 'c1cc[nH]c1', 'C1CC2CCC1C2', 'F[C@H](Cl)Br']
ap = generators.GetAtomPairGenerator()
tt = generators.GetTopologicalTorsionGenerator()
rows = []
for source in sources:
    mol = Chem.MolFromSmiles(source)
    assert mol is not None, source
    rows.append(dict(smiles=source,
        atom_pair=sorted(ap.GetSparseCountFingerprint(mol).GetNonzeroElements().items()),
        torsion=sorted(tt.GetSparseCountFingerprint(mol).GetNonzeroElements().items()),
        legacy_torsion=sorted(rdMolDescriptors.GetTopologicalTorsionFingerprint(mol).GetNonzeroElements().items())))
(root/'validation/rdkit-2026.03.1-sparse-fp-boundary.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion, rows=rows), indent=2)+'\n')
print('wrote', len(rows), 'cases')
