#!/usr/bin/env python3
"""Exercise the complete strict shard pipeline on two synthetic molecules."""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'scripts' / 'benchmarks'
sys.path.insert(0, str(BENCH))
import shard_platinum_holdout as fixed  # noqa: E402
import shard_platinum_time as timed  # noqa: E402


def run(script, *args):
    command = [sys.executable, str(BENCH / script), *(str(arg) for arg in args)]
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'{script} failed: {result.stderr[-3000:]}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profile', choices=('corrected', 'mixed_checked'), default='corrected')
    args = parser.parse_args()
    wheel_hash = hashlib.sha256(args.wheel.read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        sdf = root / 'synthetic.sdf'
        writer = Chem.SDWriter(str(sdf))
        rows = []
        for index, smiles in enumerate(('F[C@](Cl)(Br)I', 'F/C=C/F')):
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
                raise RuntimeError('synthetic reference embedding failed')
            name = f'synthetic-{index}'
            mol.SetProp('_Name', name)
            writer.write(mol)
            rows.append({'index': index, 'name': name, 'split': 'holdout',
                         'scaffold_sha256': f'synthetic-scaffold-{index}'})
        writer.close()
        source_hash = hashlib.sha256(sdf.read_bytes()).hexdigest()
        split = root / 'split.json'
        split.write_bytes(fixed.encoded({'source_sha256': source_hash,
                                         'invalid_exclusions': [], 'rows': rows}))
        audit = root / 'audit.json'
        audit.write_bytes(fixed.encoded({'source_sha256': source_hash,
                                         'split_sha256': fixed.digest(split.read_bytes()),
                                         'invalid_exclusion_count': 0, 'overlaps': []}))
        fixed_dir = root / 'fixed'
        fixed.prepare(split, audit, fixed_dir, 2, wheel_hash, args.profile)
        fixed_results = []
        for index in range(2):
            suffix = f'{index:03d}'
            result = root / f'fixed-{suffix}.json'
            shard_audit = fixed_dir / f'audit-{suffix}.json'
            fixed_wrapper = ('platinum_mixed_checked_development.py' if args.profile == 'mixed_checked'
                             else 'platinum_development_strict.py')
            run(fixed_wrapper, '--sdf', sdf,
                '--split', fixed_dir / f'split-{suffix}.json',
                '--invalid-exclusion-audit', shard_audit, '--partition', 'holdout',
                '--output', result, '--limit', 1, '--count', 20, '--seed', 42,
                '--profile', 'corrected', '--geometry', '--repair-iters', 400)
            fixed.seal(fixed_dir / 'manifest.json', index, args.wheel, shard_audit, result)
            fixed_results.append(result)
        baseline = root / 'fixed-merged.json'
        fixed_data = fixed.merge(fixed_dir / 'manifest.json', fixed_results, baseline)
        if len(fixed_data['results']) != 2 or any(row.get('status') != 'ok' for row in fixed_data['results']):
            raise RuntimeError('synthetic fixed-count scoring failed')
        run('summarize_platinum_holdout.py', baseline, split, '--output', root / 'fixed-summary.json', '--bootstrap', 100)
        run('summarize_platinum_stereo_strict.py', baseline, '--output', root / 'fixed-stereo.json')
        timed_dir = root / 'timed'
        timed.prepare(fixed_dir / 'manifest.json', baseline, timed_dir)
        timed_results = []
        for index in range(2):
            suffix = f'{index:03d}'
            subset = timed_dir / f'baseline-{suffix}.json'
            result = root / f'timed-{suffix}.json'
            timed_wrapper = ('platinum_time_mixed_strict.py' if args.profile == 'mixed_checked'
                             else 'platinum_time_strict.py')
            run(timed_wrapper, '--sdf', sdf, '--baseline', subset,
                '--output', result, '--limit', 1, '--seconds', 0.2,
                '--max-candidates', 20, '--repair-iters', 400, '--geometry')
            timed.seal(timed_dir / 'manifest.json', index, args.wheel, subset, result)
            timed_results.append(result)
        timed_merged = root / 'timed-merged.json'
        timed_data = timed.merge(timed_dir / 'manifest.json', timed_results, timed_merged)
        if len(timed_data['results']) != 2 or any(row.get('error') or row.get('engine_errors') for row in timed_data['results']):
            raise RuntimeError('synthetic timed scoring failed')
        run('summarize_platinum_time_holdout.py', timed_merged, split, '--output', root / 'timed-summary.json', '--bootstrap', 100)
        run('summarize_platinum_stereo_strict.py', timed_merged, '--output', root / 'timed-stereo.json')
        report = {'rdkit_version': rdBase.rdkitVersion, 'platform': platform.platform(), 'profile': args.profile,
                  'wheel_sha256': wheel_hash, 'synthetic_records': 2,
                  'fixed_stereo': json.loads((root / 'fixed-stereo.json').read_text()),
                  'timed_stereo': json.loads((root / 'timed-stereo.json').read_text()),
                  'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'synthetic_records': report['synthetic_records'],
                      'wheel_sha256': wheel_hash, 'platform': report['platform']}))


if __name__ == '__main__':
    main()
