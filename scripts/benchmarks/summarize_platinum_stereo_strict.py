#!/usr/bin/env python3
"""Separate contradicted and unevaluable declared stereo in strict scorer output."""
import argparse
import json
from pathlib import Path


def summarize(data):
    if data.get('stereo_policy') != 'rdkit-declared-stereo-strict-v1':
        raise ValueError('not a completed strict-stereo scoring run')
    result = {'rows': len(data['results']), 'stereo_policy': data['stereo_policy'], 'engines': {}}
    for engine in ('chematic', 'rdkit'):
        counts = {'conformers': 0, 'contradicted': 0, 'unevaluable': 0,
                  'invalid_either': 0, 'records_with_no_conformers': 0}
        for row in data['results']:
            if f'{engine}_stereo' in row:
                entries = row[f'{engine}_stereo']
            else:
                entries = row.get(engine, {}).get('stereo', [])
            if not entries:
                counts['records_with_no_conformers'] += 1
            for stereo in entries:
                required = ('violated_atoms', 'violated_bonds', 'unevaluable_atoms',
                            'unevaluable_bonds', 'atoms', 'bonds')
                if any(key not in stereo for key in required):
                    raise ValueError('conformer is missing strict-stereo diagnostics')
                violated = bool(stereo['violated_atoms'] or stereo['violated_bonds'])
                unknown = bool(stereo['unevaluable_atoms'] or stereo['unevaluable_bonds'])
                if (set(stereo['atoms']) != set(stereo['violated_atoms']) | set(stereo['unevaluable_atoms'])
                        or set(stereo['bonds']) != set(stereo['violated_bonds']) | set(stereo['unevaluable_bonds'])):
                    raise ValueError('strict-stereo combined status is inconsistent')
                counts['conformers'] += 1
                counts['contradicted'] += violated
                counts['unevaluable'] += unknown
                counts['invalid_either'] += violated or unknown
        result['engines'][engine] = counts
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(json.loads(args.input.read_text()))
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
