#!/usr/bin/env python3
"""Count timed conformers whose declared stereo and recovery were actually scored."""
import argparse
import json
from pathlib import Path


def summarize(data):
    rows = data['results']
    output = {'rows': len(rows), 'engines': {}}
    for engine in ('chematic', 'rdkit'):
        count = {'engine_records': 0, 'missing_engine_records': 0,
                 'row_errors': 0, 'generation_errors': 0, 'scoring_errors': 0,
                 'accepted_conformers': 0, 'stereo_scored_conformers': 0,
                 'unscored_accepted_conformers': 0, 'records_with_unscored_accepted': 0}
        for row in rows:
            if row.get('error'):
                count['row_errors'] += 1
            error = row.get('engine_errors', {}).get(engine)
            if error:
                stage = error.get('stage')
                if stage not in ('generation', 'scoring'):
                    raise ValueError('unknown engine error stage')
                count[stage + '_errors'] += 1
            record = row.get(engine)
            if record is None:
                count['missing_engine_records'] += 1
                continue
            count['engine_records'] += 1
            accepted = record['accepted']
            if not isinstance(accepted, int) or accepted < 0:
                raise ValueError('invalid accepted count')
            stereo = record.get('stereo', [])
            if not isinstance(stereo, list) or len(stereo) > accepted:
                raise ValueError('invalid stereo score count')
            unscored = accepted - len(stereo)
            count['accepted_conformers'] += accepted
            count['stereo_scored_conformers'] += len(stereo)
            count['unscored_accepted_conformers'] += unscored
            count['records_with_unscored_accepted'] += bool(unscored)
        output['engines'][engine] = count
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    result = summarize(json.loads(args.input.read_text()))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
