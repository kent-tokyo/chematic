#!/usr/bin/env python3
"""Cross-check retained v1.2.2 release and published-benchmark identities.

This is an offline consistency check. It does not query registries or verify
remote responses; the dated records preserve the earlier live measurements.
"""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.release.check_pypi_dist import check_layout

RESULTS = ROOT / 'validation/results'
REGISTRY = RESULTS / 'published-benchmark-wheel-registry-verification-2026-11-21.json'
CURRENT_REGISTRY = RESULTS / 'published-benchmark-wheel-registry-cos063-2026-11-21.json'
PYPI = RESULTS / 'pypi-v1.2.2-artifact-verification-2026-11-21.json'
CRATES = RESULTS / 'crates-v1.2.2-artifact-verification-2026-11-21.json'
NPM = RESULTS / 'npm-v1.2.2-artifact-verification-2026-11-21.json'
CHANNELS = RESULTS / 'release-channel-verification-v1.2.2.json'
GATES = RESULTS / 'pypi-v1.2.2-published-gates-2026-11-21.json'
MCS = ROOT / 'benchmarks/2026-11-21-v122-published-mcs-uniform.json'
SPEED = ROOT / 'benchmarks/2026-11-21-v122-published-speed-frontier.json'
CURRENT_MCS = ROOT / 'benchmarks/2026-11-21-v122-published-mcs-uniform-cos063.json'
CURRENT_SPEED = ROOT / 'benchmarks/2026-11-21-v122-published-speed-cos063.json'
TAG_COMMIT = 'd7962f9501da9ad2f2cdc5eaddf6d665d79d60bc'
CRATE_NAMES = {'chematic-core', 'chematic-perception', 'chematic-cip', 'chematic-smarts',
               'chematic-smiles', 'chematic-rxn', 'chematic-fp', 'chematic-iupac',
               'chematic-chem', 'chematic-ff', 'chematic-3d', 'chematic-crystal',
               'chematic-depict', 'chematic-inchi', 'chematic-mol', 'chematic',
               'chematic-cli', 'chematic-ewald', 'chematic-mcp', 'chematic-wasm'}
VERSIONS = {'v120': ('chematic', '1.2.0'), 'v121': ('chematic', '1.2.1'),
            'v122': ('chematic', '1.2.2'), 'rdkit': ('rdkit', '2026.9.2'),
            'cosmolkit': ('cosmolkit', '0.6.2')}
EXPECTED_GATES = {
    'Recover v1.2.2 PyPI publication',
    *(f'recover-v122-smoke / Python wheel runtime ({platform})' for platform in
      ('macos-latest / 3.13', 'ubuntu-24.04 / 3.9', 'windows-latest / 3.13')),
    *(f'recover-v122-chemistry / {lane}' for lane in (
        'xsmarts-autoconf linux-x86_64-cp312-sdist', 'xsmarts-autoconf macos-arm64-cp312',
        'xsmarts-autoconf windows-amd64-cp312', 'chematic linux-x86_64-cp39',
        'chematic macos-arm64-cp312', 'chematic macos-x86_64-cp312', 'chematic windows-amd64-cp312',
        'RDKit compare macos-arm64-cp312', 'RDKit compare macos-x86_64-cp312',
        'RDKit compare windows-amd64-cp312', 'RDKit compare linux-x86_64-cp39')),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_measurements(mcs, speed):
    """Verify rotating blocks and the reported MCS comparison gates."""
    require(mcs['complete'] is True and mcs['uniform_smarts'] is True and
            speed['complete'] is True, 'benchmark record incomplete')
    arms = ('v120', 'v121', 'v122', 'rdkit', 'cosmolkit')
    speed_arms = ('older', 'baseline', 'candidate', 'rdkit', 'cosmolkit')
    require(len(mcs['blocks']) == len(speed['blocks']) == 20, 'benchmark block count differs')
    for index, (mcs_block, speed_block) in enumerate(zip(mcs['blocks'], speed['blocks'], strict=True)):
        offset = index % 5
        require(mcs_block['order'] == list(arms[offset:] + arms[:offset]) and
                speed_block['order'] == list(speed_arms[offset:] + speed_arms[:offset]),
                f'benchmark rotation differs in block {index + 1}')
    require(set(speed['agreement']) == {'tpsa', 'num_rings', 'morgan2_2048'}, 'speed metrics differ')
    comparison_names = {'candidate_vs_baseline', 'candidate_vs_older', 'older_vs_rdkit',
                        'baseline_vs_rdkit', 'candidate_vs_rdkit', 'cosmolkit_vs_rdkit'}
    for metric, comparisons in speed['agreement'].items():
        require(set(comparisons) == comparison_names and all(row['passed'] is True and row['rows'] == 1000 and
                row['mismatches'] == 0 for row in comparisons.values()), f'agreement differs: {metric}')
    expected = {'v120': [True, False, True, False, False], 'v121': [True, True, True, True, False],
                'v122': [True] * 5, 'rdkit': [True] * 5, 'cosmolkit': [True] * 5}
    require(len(mcs['cases']) == 5, 'MCS input set differs')
    for arm in arms:
        cases = mcs['case_summary'][arm]
        require(len(cases) == 5, f'MCS case count differs: {arm}')
        for index, summary in enumerate(cases):
            comparable = []
            for block in mcs['blocks']:
                actual, reference = block[arm]['validation'][index], block['rdkit']['validation'][index]
                comparable.append('error' not in actual and 'error' not in reference and
                                  not actual['limited'] and not reference['limited'] and
                                  len(actual['matches']) == len(mcs['cases'][index]) and
                                  len(reference['matches']) == len(mcs['cases'][index]) and
                                  all(actual['matches']) and all(reference['matches']) and
                                  actual['same_query_predicates_as_rdkit'] and
                                  actual['atoms'] == reference['atoms'] and actual['bonds'] == reference['bonds'])
            require(summary['case'] == index and summary['passes_all_blocks'] == all(comparable) == expected[arm][index],
                    f'MCS output contract differs: {arm}/{index}')
            require(summary['first_ns_median'] == statistics.median(block[arm]['first_ns'][index] for block in mcs['blocks']) and
                    summary['repeat_ns_median'] == statistics.median(block[arm]['repeated_ns'][index] for block in mcs['blocks']),
                    f'MCS timing summary differs: {arm}/{index}')
            if arm != 'v122':
                comparison = mcs['paired_repeated_search'][arm][index]
                require(comparison['comparable_by_query_predicates_and_input_matches'] == expected[arm][index] and
                        (comparison['competitor_over_v122'] is not None) == expected[arm][index],
                        f'MCS ratio gate differs: {arm}/{index}')


def check(registry, pypi, crates, npm, channels, gates, mcs, speed, versions=VERSIONS):
    check_measurements(mcs, speed)
    require(channels['release'] == {'version': '1.2.2', 'tag': 'v1.2.2', 'commit': TAG_COMMIT},
            'release tag identity differs')
    require(pypi['version'] == '1.2.2' and pypi['tag_commit'] == TAG_COMMIT and
            pypi['original_ci_artifacts_equal_pypi'] is True and
            pypi['original_workflow_run'] == 44831147798 and
            pypi['wheel_stub_gate'] == '25 wheels passed', 'PyPI release evidence differs')
    check_layout(set(pypi['files']), '1.2.2')
    require(all(re.fullmatch(r'[0-9a-f]{64}', value) for value in pypi['files'].values()),
            'PyPI hash invalid')
    crate_rows = crates['files']
    require(len(crate_rows) == 20 and {row['name'] for row in crate_rows} == CRATE_NAMES,
            'crate set differs')
    for row in crate_rows:
        require(row['version'] == '1.2.2' and row['vcs']['git']['sha1'] == TAG_COMMIT and
                row['vcs']['path_in_vcs'] == 'crates/' + row['name'] and row['bytes'] > 0 and
                re.fullmatch(r'[0-9a-f]{64}', row['sha256']) is not None,
                f"crate provenance differs: {row['name']}")
    require(npm['name'] == '@kent-tokyo/chematic' and npm['version'] == '1.2.2' and
            npm['registry_git_head'] == TAG_COMMIT and npm['computed_integrity_matches'] is True and
            re.fullmatch(r'[0-9a-f]{64}', npm['sha256']) is not None,
            'npm release evidence differs')
    integrity = npm['registry_sha512_integrity']
    require(integrity.startswith('sha512-') and len(base64.b64decode(integrity[7:], validate=True)) == 64,
            'npm integrity invalid')
    require(npm['smoke']['formula'] == 'C6H6' and npm['smoke']['atoms'] == 6,
            'npm smoke differs')
    require(gates['databaseId'] == pypi['recovery_workflow_run'] == 44887442455 and
            gates['event'] == 'workflow_dispatch' and gates['status'] == 'completed' and
            gates['conclusion'] == 'success' and gates['headSha'] == 'e3bbb16fe6916f744df325f494fd7bc64577ffa5',
            'published gate run differs')
    jobs = gates['jobs']
    require(len(jobs) == len(EXPECTED_GATES) and {job['name'] for job in jobs} == EXPECTED_GATES and
            all(job['status'] == 'completed' and job['conclusion'] == 'success' and
                job['url'].startswith(gates['url'] + '/job/') for job in jobs),
            'published runtime or chemistry gate differs')
    rows = registry['wheels']
    by_arm = {row['arm']: row for row in rows}
    require(len(rows) == 5 and set(by_arm) == set(versions), 'benchmark wheel set differs')
    require(registry['v122_complete_file_set_matches'] is True, 'PyPI file-set comparison absent')
    for arm, (package, version) in versions.items():
        row, measured = by_arm[arm], mcs['wheel_verification'][arm]
        require(row['package'] == package and row['version'] == version and
                row['pypi_url'] == f'https://pypi.org/pypi/{package}/{version}/json' and
                re.fullmatch(r'[0-9a-f]{64}', row['pypi_response_sha256']) is not None and
                re.fullmatch(r'[0-9a-f]{64}', row['sha256']) is not None and
                row['installed_files_verified_in_mcs'] == measured['installed_files_verified'] > 0 and
                row['filename'] == measured['filename'] and row['sha256'] == measured['sha256'],
                f'benchmark registry/MCS wheel differs: {arm}')
        if arm.startswith('v'):
            saved = speed['wheels'][{'v120': 'older', 'v121': 'baseline', 'v122': 'candidate'}[arm]]
            require(Path(saved['path']).name == row['filename'] and saved['sha256'] == row['sha256'],
                    f'speed wheel differs: {arm}')
        else:
            require(speed['comparator_wheels'].get(row['filename']) == row['sha256'],
                    f'speed comparator wheel differs: {arm}')
    require(pypi['files'].get(by_arm['v122']['filename']) == by_arm['v122']['sha256'],
            'v1.2.2 measured wheel differs from release file')


def main():
    common = [json.loads(path.read_text()) for path in (PYPI, CRATES, NPM, CHANNELS, GATES)]
    for registry_path, mcs_path, speed_path, cos_version in (
        (REGISTRY, MCS, SPEED, '0.6.2'),
        (CURRENT_REGISTRY, CURRENT_MCS, CURRENT_SPEED, '0.6.3'),
    ):
        registry = json.loads(registry_path.read_text())
        for prefix, path in [('mcs', mcs_path), ('speed', speed_path), ('release', PYPI)]:
            require(registry[f'{prefix}_record'] == str(path.relative_to(ROOT)) and
                    registry[f'{prefix}_record_sha256'] == digest(path), f'{prefix} record hash differs')
        versions = dict(VERSIONS, cosmolkit=('cosmolkit', cos_version))
        check(registry, *common, json.loads(mcs_path.read_text()), json.loads(speed_path.read_text()), versions)
        if cos_version == '0.6.3':
            latest = registry['latest_stable_at_check']
            require({arm: latest[arm]['version'] for arm in ('rdkit', 'cosmolkit')} ==
                    {'rdkit': '2026.9.2', 'cosmolkit': '0.6.3'} and
                    all(latest[arm]['pypi_url'] == f'https://pypi.org/pypi/{arm}/json' and
                        re.fullmatch(r'[0-9a-f]{64}', latest[arm]['pypi_response_sha256']) is not None
                        for arm in ('rdkit', 'cosmolkit')), 'latest comparator record differs')
    print('v1.2.2 release and historical/current published benchmark wheel identities are consistent')


if __name__ == '__main__':
    main()
