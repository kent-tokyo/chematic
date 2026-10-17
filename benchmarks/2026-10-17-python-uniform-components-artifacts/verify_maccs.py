#!/usr/bin/env python3
"""Dump exact first/repeated MACCS bytes, retaining every parse/evaluation error."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    import chematic
    root = Path('/private/tmp/chematic-python-asa-preparation-20261017')
    paths = [root / 'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi',
             root / 'validation/cosmolkit_comparison/corpora/stress-v1.smi']
    with args.output.open('w') as out:
        for path in paths:
            data = path.read_bytes()
            rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()]
            out.write(json.dumps({'corpus': path.name, 'rows': len(rows),
                                  'sha256': hashlib.sha256(data).hexdigest(),
                                  'version': chematic.__version__}) + '\n')
            for index, smiles in enumerate(rows):
                row = {'corpus': path.name, 'index': index}
                try:
                    mol = chematic.from_smiles(smiles)
                    row['first'] = bytes(mol.maccs()).hex()
                    row['repeated'] = bytes(mol.maccs()).hex()
                except Exception as error:
                    row['error'] = type(error).__name__ + ': ' + str(error)
                out.write(json.dumps(row, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
