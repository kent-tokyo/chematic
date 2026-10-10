"""Freeze V3000 Avalon fingerprints with RDKit 2026.03.1."""
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Avalon import pyAvalonTools
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
root = Path(__file__).resolve().parents[1]
sources = ['CCO', 'c1ccccc1', 'c1cc[nH]c1', '[NH4+]', 'CC(=O)[O-]',
           '[CH3]', '[O][O]', '[13CH3][15NH2]', 'F[C@H](Cl)Br', '*CC*',
           'CC.CN', 'C#N', 'C1CC2CCC1C2', 'CS(=O)(=O)C']
rows = []
for source in sources:
    mol = Chem.MolFromSmiles(source)
    assert mol is not None, source
    block = Chem.MolToMolBlock(mol, forceV3000=True)
    rows.append(dict(smiles=source, block=block,
        on_bits=list(pyAvalonTools.GetAvalonFP(block, isSmiles=False, nBits=512).GetOnBits())))
(root/'validation/rdkit-2026.03.1-avalon-molblock-boundary.json').write_text(
    json.dumps(dict(rdkit_version=rdBase.rdkitVersion, rows=rows), indent=2)+'\n')
print('wrote', len(rows), 'cases')
