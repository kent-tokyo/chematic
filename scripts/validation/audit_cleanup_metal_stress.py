#!/usr/bin/env python3
"""Audit metal-containing stress rows and explicit-H input variants."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import zipfile

import chematic
from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / 'validation/cosmolkit_comparison/corpora/stress-v1.smi'
METALS = {3, 4, 11, 12, 13, *range(19, 32), *range(37, 51), *range(55, 84), *range(87, 117)}


def hydrogen_count(formula):
    # Formula uses element symbols; exclude He, Hg, Hf, Ho and Hs.
    return sum(int(count or 1) for count in re.findall(r'H(?![a-z])(\d*)', formula))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--candidate-commit', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if chematic.__version__ != '1.2.2' or rdBase.rdkitVersion != '2026.09.1':
        raise ValueError('expected chematic 1.2.2 and RDKit 2026.09.1')
    root = Path(chematic.__file__).parent.parent
    with zipfile.ZipFile(args.wheel) as wheel:
        files = [n for n in wheel.namelist() if n.startswith('chematic/') and not n.endswith('/')]
        if not files or any(not (root / n).is_file() or (root / n).read_bytes() != wheel.read(n) for n in files):
            raise ValueError('installed package differs from wheel')
    RDLogger.DisableLog('rdApp.*')
    counts = Counter()
    rows = []
    for index, line in enumerate(CORPUS.read_text().splitlines()):
        source = line.split()[0]
        mol = Chem.MolFromSmiles(source)
        if mol is None:
            counts['rdkit_input_error'] += 1
            continue
        if not any(atom.GetAtomicNum() in METALS for atom in mol.GetAtoms()):
            continue
        counts['metal_rows'] += 1
        expanded = Chem.MolToSmiles(Chem.AddHs(mol), canonical=False, isomericSmiles=True)
        check = Chem.MolFromSmiles(expanded)
        if check is None or Chem.MolToSmiles(check, isomericSmiles=True) != Chem.MolToSmiles(mol, isomericSmiles=True):
            counts['rdkit_input_variant_mismatch'] += 1
            rows.append({'index': index, 'input': source, 'issue': 'rdkit_input_variant_mismatch'})
            continue
        variants = []
        for text in (source, expanded):
            try:
                mol_in = chematic.from_smiles(text)
                atoms = mol_in.atom_table
                h_in = hydrogen_count(mol_in.formula)
                q_in = sum(a[2] for a in atoms)
                result = mol_in.cleanup()
                atoms = result.atom_table
                variants.append({'outcome': 'ok', 'output': result.rdkit_smiles,
                                 'input_hydrogens': h_in, 'output_hydrogens': hydrogen_count(result.formula),
                                 'input_charge': q_in, 'output_charge': sum(a[2] for a in atoms)})
            except chematic.ChematicInputError as error:
                variants.append({'outcome': 'refused', 'category': error.category, 'code': error.code})
            except Exception as error:
                variants.append({'outcome': 'error', 'type': type(error).__name__, 'message': str(error)})
        counts['compared'] += 1
        comparable = variants[0] == variants[1]
        hydrogen_conserved = all(v['outcome'] != 'ok' or v['input_hydrogens'] == v['output_hydrogens'] for v in variants)
        if not comparable:
            counts['different_variants'] += 1
        if not hydrogen_conserved:
            counts['hydrogen_not_conserved'] += 1
        if all(v['outcome'] == 'refused' for v in variants):
            counts['both_refused'] += 1
        if any(v['outcome'] == 'error' for v in variants):
            counts['unexpected_errors'] += 1
        if not comparable or not hydrogen_conserved or any(v['outcome'] != 'ok' for v in variants):
            rows.append({'index': index, 'input': source, 'expanded': expanded, 'variants': variants,
                         'same_outcome': comparable, 'hydrogen_conserved': hydrogen_conserved})
    report = {'candidate_commit': args.candidate_commit, 'wheel': args.wheel.name,
              'wheel_sha256': hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              'corpus_sha256': hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
              'rdkit_version': rdBase.rdkitVersion, 'installed_package_files_verified': len(files),
              'method': 'metal element set follows cleanup; RDKit AddHs then reparse; compare exact cleanup outcomes and explicit/implicit hydrogen counts; charge changes reported separately',
              'counts': dict(counts), 'exception_rows': rows}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(dict(counts))


if __name__ == '__main__':
    main()
