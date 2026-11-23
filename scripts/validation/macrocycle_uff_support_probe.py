#!/usr/bin/env python3
"""Compare RDKit and chematic UFF all-atom parameter availability, not construction."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sdf', type=Path)
    parser.add_argument('support', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    support = json.loads(args.support.read_text())
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_hash != support['source_sha256']:
        raise ValueError('source hash mismatch')
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    rows = []
    for item in support['results']:
        row = {'index': item['index'], 'name': item['name']}
        try:
            if item.get('graph_matches') is not True:
                raise ValueError('graph identity was not verified')
            mol = supplier[item['index']]
            if mol is None or mol.GetProp('_Name') != item['name']:
                raise ValueError('reference parse or identity failure')
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
            graph = Chem.MolFromSmiles(smiles)
            rdkit_mol = Chem.AddHs(graph)
            chematic_mol = chematic.from_smiles(smiles).add_hydrogens()
            row['rdkit_uff_all_params'] = bool(AllChem.UFFHasAllMoleculeParams(rdkit_mol))
            row['chematic_uff_all_params'] = bool(chematic_mol.rdkit_uff_has_all_params())
        except Exception as exc:
            row['error'] = repr(exc)
        rows.append(row)
    counts = Counter((r.get('rdkit_uff_all_params'), r.get('chematic_uff_all_params')) for r in rows)
    result = {
        'source_sha256': source_hash,
        'support_sha256': hashlib.sha256(args.support.read_bytes()).hexdigest(),
        'rdkit_version': rdBase.rdkitVersion,
        'partition': support['partition'],
        'counts': [{'rdkit_uff_all_params': a, 'chematic_uff_all_params': b, 'records': n}
                   for (a, b), n in counts.items()],
        'results': rows,
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['counts']))


if __name__ == '__main__':
    main()
