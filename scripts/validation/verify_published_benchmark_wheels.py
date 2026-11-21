#!/usr/bin/env python3
"""Compare measured wheel identities to saved official PyPI JSON responses."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MCS = ROOT / 'benchmarks/2026-11-21-v122-published-mcs-uniform.json'
SPEED = ROOT / 'benchmarks/2026-11-21-v122-published-speed-frontier.json'
RELEASE = ROOT / 'validation/results/pypi-v1.2.2-artifact-verification-2026-11-21.json'
SOURCES = {'v120': ('chematic', '1.2.0'), 'v121': ('chematic', '1.2.1'),
           'v122': ('chematic', '1.2.2'), 'rdkit': ('rdkit', '2026.9.2'),
           'cosmolkit': ('cosmolkit', '0.6.2')}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arm in SOURCES:
        parser.add_argument(f'--{arm}-metadata', required=True, type=Path)
    parser.add_argument('--cosmolkit-version', default='0.6.2')
    parser.add_argument('--mcs-record', type=Path, default=MCS)
    parser.add_argument('--speed-record', type=Path, default=SPEED)
    parser.add_argument('--rdkit-current-metadata', type=Path)
    parser.add_argument('--cosmolkit-current-metadata', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    mcs_path = args.mcs_record.resolve()
    speed_path = args.speed_record.resolve()
    mcs = json.loads(mcs_path.read_text())
    speed = json.loads(speed_path.read_text())
    release = json.loads(RELEASE.read_text())
    rows = []
    for arm, (package, version) in SOURCES.items():
        if arm == 'cosmolkit':
            version = args.cosmolkit_version
        path = getattr(args, f'{arm}_metadata')
        metadata = json.loads(path.read_text())
        wheel = mcs['wheel_verification'][arm]
        if metadata['info']['name'].lower() != package or metadata['info']['version'] != version:
            raise ValueError(f'wrong registry metadata: {arm}')
        if arm == 'v122':
            published = {item['filename']: item['digests']['sha256'] for item in metadata['urls']}
            if len(published) != len(metadata['urls']) or published != release['files']:
                raise ValueError('current v1.2.2 PyPI file/hash set differs from release evidence')
        entries = [item for item in metadata['urls'] if item['filename'] == wheel['filename']]
        if len(entries) != 1 or entries[0]['digests']['sha256'] != wheel['sha256']:
            raise ValueError(f'registry wheel missing or hash differs: {arm}')
        if arm.startswith('v'):
            key = {'v120': 'older', 'v121': 'baseline', 'v122': 'candidate'}[arm]
            saved = speed['wheels'][key]
            speed_name = Path(saved['path']).name
            speed_hash = saved['sha256']
        else:
            speed_name = wheel['filename']
            speed_hash = speed['comparator_wheels'].get(speed_name)
        if speed_name != wheel['filename'] or speed_hash != wheel['sha256']:
            raise ValueError(f'speed and MCS wheel differ: {arm}')
        rows.append({'arm': arm, 'package': package, 'version': version,
                     'filename': wheel['filename'], 'sha256': wheel['sha256'],
                     'pypi_url': f'https://pypi.org/pypi/{package}/{version}/json',
                     'pypi_response_sha256': digest(path),
                     'installed_files_verified_in_mcs': wheel['installed_files_verified']})
    latest = {}
    for arm in ('rdkit', 'cosmolkit'):
        path = getattr(args, f'{arm}_current_metadata')
        if path is None:
            continue
        meta = json.loads(path.read_text())
        expected_version = SOURCES[arm][1] if arm == 'rdkit' else args.cosmolkit_version
        if meta['info']['name'].lower() != arm or meta['info']['version'] != expected_version:
            raise ValueError(f'{arm} latest stable version differs from measured wheel')
        latest[arm] = {'version': expected_version, 'pypi_url': f'https://pypi.org/pypi/{arm}/json',
                       'pypi_response_sha256': digest(path)}
    report = {'checked_at': datetime.now(timezone.utc).isoformat(),
              'mcs_record': str(mcs_path.relative_to(ROOT)), 'mcs_record_sha256': digest(mcs_path),
              'speed_record': str(speed_path.relative_to(ROOT)), 'speed_record_sha256': digest(speed_path),
              'release_record': str(RELEASE.relative_to(ROOT)), 'release_record_sha256': digest(RELEASE),
              'v122_complete_file_set_matches': True,
              'method': 'compare exact wheel filename and SHA-256 in PyPI JSON with both measurement records',
              'wheels': rows, 'latest_stable_at_check': latest}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(f"PyPI metadata matches both benchmark records for {len(rows)} wheels")


if __name__ == '__main__':
    main()
