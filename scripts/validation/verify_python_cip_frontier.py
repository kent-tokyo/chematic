#!/usr/bin/env python3
"""Retain complete CIP results and normalized labels on both pinned corpora.

The raw source API and the RDKit-model source API are separate lanes. Conversion
into the RDKit model belongs inside the timed operation in any speed benchmark.
RDKit assigns properties in place; label extraction here is outside all timing.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPORA = (ROOT / 'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi',
           ROOT / 'validation/cosmolkit_comparison/corpora/stress-v1.smi')


def source_labels(values):
    return sorted(['b', *sorted(v['bond_atoms']), v['descriptor']]
                   if 'bond_idx' in v else ['a', v['atom_idx'], v['descriptor']]
                  for v in values)


def rdkit_labels(mol):
    return sorted([['a', a.GetIdx(), a.GetProp('_CIPCode')]
                   for a in mol.GetAtoms() if a.HasProp('_CIPCode')]
                  + [['b', *sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())), b.GetProp('_CIPCode')]
                     for b in mol.GetBonds() if b.HasProp('_CIPCode')])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine', choices=('chematic', 'rdkit'), required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.engine == 'chematic':
        import chematic
        parse = chematic.from_smiles
    else:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdCIPLabeler
        RDLogger.DisableLog('rdApp.*')
        parse = Chem.MolFromSmiles
    with args.output.open('w') as out:
        for corpus in CORPORA:
            data = corpus.read_bytes()
            rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()]
            out.write(json.dumps({'corpus': corpus.name, 'rows': len(rows),
                                  'sha256': hashlib.sha256(data).hexdigest()}) + '\n')
            for i, smiles in enumerate(rows):
                row = {'corpus': corpus.name, 'index': i}
                try:
                    m = parse(smiles)
                    if m is None:
                        raise ValueError('parse failed')
                    if args.engine == 'chematic':
                        for lane, molecule in [('raw', m), ('rdkit_model', m.rdkit_parsed())]:
                            values = molecule.cip_stereo(mode='accurate')
                            row[lane] = {'labels': source_labels(values), 'complete': values,
                                         'unresolved': molecule.cip_stereo_unresolved(),
                                         'legacy': molecule.cip_stereo(mode='legacy')}
                    else:
                        rdCIPLabeler.AssignCIPLabels(m)
                        row['labels'] = rdkit_labels(m)
                except Exception as error:
                    row['error'] = type(error).__name__ + ': ' + str(error)
                out.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
            print(corpus.name, len(rows), 'done', flush=True)


if __name__ == '__main__':
    main()
