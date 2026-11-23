#!/usr/bin/env python3
"""Resume and join a prefix checkpoint without rescoring completed records."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def rows(split):
    found = [(r['index'], r['name']) for r in split['rows'] if r.get('split') == 'holdout']
    if len({i for i, _ in found}) != len(found) or len({n for _, n in found}) != len(found):
        raise ValueError('duplicate eligible index or name')
    return sorted(found)


def checkpoint_info(split_path, audit_path, checkpoint_path):
    split, audit, checkpoint = map(read, (split_path, audit_path, checkpoint_path))
    if (audit['split_sha256'] != digest(split_path) or audit['source_sha256'] != split['source_sha256']
        or audit['invalid_exclusion_count'] != len(split['invalid_exclusions']) or audit['overlaps']):
        raise ValueError('invalid original split audit')
    if checkpoint['split_sha256'] != digest(split_path) or checkpoint['source_sha256'] != split['source_sha256']:
        raise ValueError('checkpoint does not belong to the original split')
    expected = rows(split)
    completed = [(r['index'], r['name']) for r in checkpoint['results']]
    if not completed or completed != expected[:len(completed)] or len(completed) >= len(expected):
        raise ValueError('checkpoint is not a nonempty proper prefix of the complete split')
    return split, audit, checkpoint, expected, completed


def prepare(args):
    split, audit, checkpoint, expected, completed = checkpoint_info(args.split, args.audit, args.checkpoint)
    done = {i for i, _ in completed}
    for row in split['rows']:
        if row.get('split') == 'holdout' and row['index'] in done:
            row['split'] = 'excluded_completed_checkpoint'
    split['resume'] = {'original_split_sha256': digest(args.split), 'checkpoint_sha256': digest(args.checkpoint),
                       'completed_records': len(completed), 'remaining_records': len(expected) - len(completed)}
    args.out_split.write_text(json.dumps(split, indent=2) + '\n')
    audit['original_split_sha256'] = audit['split_sha256']
    audit['split_sha256'] = digest(args.out_split)
    audit['resume_checkpoint_sha256'] = digest(args.checkpoint)
    args.out_audit.write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(split['resume']))


def merge(args):
    split, _, checkpoint, expected, completed = checkpoint_info(args.split, args.audit, args.checkpoint)
    resumed = read(args.resumed)
    if resumed['split_sha256'] != digest(args.resume_split):
        raise ValueError('resumed output has an unexpected split')
    resume_split = read(args.resume_split)
    if resume_split['resume']['checkpoint_sha256'] != digest(args.checkpoint):
        raise ValueError('resumed split points to a different checkpoint')
    if (not checkpoint['platform'].startswith('macOS-')
        or not resumed['platform'].startswith('macOS-')
        or 'arm64' not in checkpoint['platform'] or 'arm64' not in resumed['platform']):
        raise ValueError('both scoring portions must use macOS arm64')
    for key in ('rdkit','chematic','source_sha256','script_sha256','profile','geometry_checked',
                'repair_iters','candidate_limit','mmff_iterations','seed','selection'):
        if checkpoint[key] != resumed[key]:
            raise ValueError('incompatible scoring metadata: ' + key)
    if [(r['index'], r['name']) for r in resumed['results']] != expected[len(completed):]:
        raise ValueError('resumed records are missing, duplicated, or misordered')
    checkpoint['results'].extend(resumed['results'])
    checkpoint['resume_provenance'] = {'original_checkpoint_sha256': digest(args.checkpoint),
                                      'resumed_split_sha256': digest(args.resume_split),
                                      'resumed_result_sha256': digest(args.resumed),
                                      'original_platform': checkpoint['platform'],
                                      'resumed_platform': resumed['platform'],
                                      'completed_before_resume': len(completed)}
    args.output.write_text(json.dumps(checkpoint, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'records': len(checkpoint['results']), **checkpoint['resume_provenance']}))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='command', required=True)
    for cmd in ('prepare', 'merge'):
        p = sub.add_parser(cmd)
        for flag in ('split', 'audit', 'checkpoint'):
            p.add_argument('--' + flag, type=Path, required=True)
        if cmd == 'prepare':
            p.add_argument('--out-split', type=Path, required=True)
            p.add_argument('--out-audit', type=Path, required=True)
        else:
            p.add_argument('--resume-split', type=Path, required=True)
            p.add_argument('--resumed', type=Path, required=True)
            p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (prepare if args.command == 'prepare' else merge)(args)


if __name__ == '__main__':
    main()
