#!/usr/bin/env python3
"""Select unseen archive slice while excluding every earlier COD scaffold and named corpus."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scaffold(mol):
    heavy = Chem.RemoveHs(mol)
    graph = Chem.MolFromSmiles(Chem.MolToSmiles(heavy, isomericSmiles=True))
    if graph is None:
        raise ValueError('failed to reconstruct molecule')
    return hashlib.sha256(MurckoScaffold.MurckoScaffoldSmiles(mol=graph, includeChirality=False).encode()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sdf', required=True, type=Path)
    p.add_argument('--prior-sdf', action='append', type=Path, required=True)
    p.add_argument('--prior-split', action='append', type=Path, required=True)
    p.add_argument('--corpus-root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    RDLogger.DisableLog('rdApp.*')
    excluded = set()
    names = set()
    invalid = {}
    corpus_paths = set()
    for path in args.prior_split:
        data = json.loads(path.read_text())
        excluded.update(r['scaffold_sha256'] for r in data['rows'] if r.get('scaffold_sha256'))
        for entry in data.get('invalid_exclusions', []):
            invalid[entry['smiles']] = entry
        corpus_paths.update(path for path in data['excluded_manifests_sha256'] if not path.startswith('/'))
    failed_old = 0
    unrepresentable_old = 0
    for path in args.prior_sdf:
        sanitized = Chem.SDMolSupplier(str(path), removeHs=False)
        unsanitized = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
        for mol, raw in zip(sanitized, unsanitized):
            if raw is not None and raw.HasProp('_Name'):
                names.add(raw.GetProp('_Name'))
            if mol is None:
                failed_old += 1
            else:
                try:
                    excluded.add(scaffold(mol))
                    continue
                except Exception:
                    failed_old += 1
            try:
                if raw is None:
                    raise ValueError('unsanitized record unavailable')
                smiles = Chem.MolToSmiles(raw, isomericSmiles=True)
                if Chem.MolFromSmiles(smiles, sanitize=False) is None:
                    raise ValueError('unsanitized graph unavailable')
                invalid[smiles] = {'path': str(path), 'smiles': smiles}
            except Exception:
                unrepresentable_old += 1
    for relative in sorted(corpus_paths):
        path = args.corpus_root / relative
        if path.suffix == '.json':
            items = [r['smiles'] for r in json.loads(path.read_text())['molecules']]
        elif path.suffix == '.jsonl':
            items = [json.loads(line)['smiles'] for line in path.read_text().splitlines() if line.strip()]
        else:
            items = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        for smiles in items:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                invalid[smiles] = {'path': relative, 'smiles': smiles}
            else:
                try:
                    excluded.add(scaffold(mol))
                except Exception:
                    invalid[smiles] = {'path': relative, 'smiles': smiles}
    rows = []
    for index, mol in enumerate(Chem.SDMolSupplier(str(args.sdf), removeHs=False)):
        if mol is None:
            rows.append({'index': index, 'status': 'parse_failure'})
            continue
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
        try:
            key = scaffold(mol)
            split = 'excluded_known_scaffold' if key in excluded or name in names else 'holdout'
            rows.append({'index': index, 'name': name, 'scaffold_sha256': key, 'split': split})
        except Exception as exc:
            rows.append({'index': index, 'name': name, 'status': 'scaffold_failure', 'error': repr(exc)})
    out = {'source_sha256': digest(args.sdf), 'rdkit_version': rdBase.rdkitVersion,
           'prior_sdf_sha256': {str(p): digest(p) for p in args.prior_sdf},
           'prior_split_sha256': {str(p): digest(p) for p in args.prior_split},
           'corpus_sha256': {p: digest(args.corpus_root / p) for p in sorted(corpus_paths)},
           'invalid_exclusions': list(invalid.values()), 'prior_unscorable_records': failed_old,
           'prior_unrepresentable_records': unrepresentable_old,
           'rule': 'exclude prior names and canonical scaffolds across both prior COD slices and named corpora',
           'counts': dict(Counter(r.get('split', r.get('status')) for r in rows)), 'rows': rows}
    args.output.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'counts':out['counts'], 'prior_unscorable_records':failed_old,
                      'prior_unrepresentable_records':unrepresentable_old, 'invalid_exclusions':len(invalid)}))


if __name__ == '__main__':
    main()
