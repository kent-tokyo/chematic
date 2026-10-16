from pathlib import Path
import hashlib, json, os, subprocess
ROOT = Path('/private/tmp/chematic-labute-locality-20261016')
OUT = Path(__file__).resolve().parent
PYTHON = '/private/tmp/chematic-ring-frontier-env-20261016/bin/python'
plan = json.loads((OUT / 'measurement-plan.json').read_text())
manifest = json.loads((OUT / 'source-manifest.json').read_text())
comparators = json.loads((OUT / 'comparator-manifest.json').read_text())
proof = json.loads((OUT / 'parity.json').read_text())

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def integrity():
    for path, value in manifest['candidate_sha256'].items(): assert sha(ROOT / path) == value, path
    for value in comparators.values():
        for path, expected in value['native_files'].items(): assert sha(Path(path)) == expected, path
    for arm, value in proof['wheels'].items():
        assert sha(Path(value['wheel'])) == value['sha256']
        for path, expected in value['extensions'].items(): assert sha(OUT / arm / path) == expected

assert proof['outputs']['baseline'] == proof['outputs']['candidate']
assert '0 failed' in (OUT / 'native-tests.log').read_text()
assert 'error:' not in (OUT / 'clippy.log').read_text()
assert 'Refactor budget OK' in (OUT / 'refactor.log').read_text()
assert not (OUT / 'format.log').read_bytes()
assert 'failed' not in (OUT / 'python-tests-final.log').read_text().splitlines()[-1]
for offset in [5000]:
    integrity()
    args = [PYTHON, str(OUT / 'run-python-safe.py'), '--corpus', str(ROOT / plan['corpus']), '--offset', str(offset), '--limit', str(plan['rows_each']), '--blocks', str(plan['blocks']), '--base-revision', plan['base'], '--source-manifest', str(OUT / 'source-manifest.json'), '--output', str(OUT / f'frontier-{offset}.json')]
    for arm in ['baseline', 'candidate']:
        args += ['--' + arm + '-wheel', proof['wheels'][arm]['wheel'], '--' + arm + '-path', str(OUT / arm)]
    with (OUT / f'frontier-{offset}.log').open('w') as log:
        subprocess.run(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    integrity()
    data = json.loads((OUT / f'frontier-{offset}.json').read_text())
    assert data['complete'] and len(data['blocks']) == plan['blocks']
    for engine, value in comparators.items(): assert data['blocks'][0][engine]['version'] == value['runtime_version']
    print(offset, 'complete', flush=True)
