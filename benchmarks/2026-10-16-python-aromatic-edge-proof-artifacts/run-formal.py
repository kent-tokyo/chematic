from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path('/private/tmp/chematic-aromatic-proof-20261016')
OUT = Path(__file__).resolve().parent
PYTHON = '/private/tmp/chematic-ring-frontier-env-20261016/bin/python'
source_bytes = (OUT / 'source-manifest.json').read_bytes()
manifest = json.loads(source_bytes)
comparator_bytes = (OUT / 'comparator-manifest.json').read_bytes()
comparators = json.loads(comparator_bytes)
proof = json.loads((OUT / 'wheel-proof-and-cache-cleanup.json').read_text())['artifacts']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def integrity():
    assert (OUT / 'source-manifest.json').read_bytes() == source_bytes
    assert (OUT / 'comparator-manifest.json').read_bytes() == comparator_bytes
    for name, expected in manifest['candidate_source_sha256'].items():
        assert sha(ROOT / name) == expected, name
    for data in comparators.values():
        for name, expected in data['shared_libraries_sha256'].items():
            assert sha(Path(name)) == expected, name
    for arm, data in proof.items():
        assert sha(next((OUT / (arm + '-wheel')).glob('*.whl'))) == data['wheel_sha256']
        assert sha(next((OUT / arm).rglob('*.so'))) == data['extension_sha256']

assert '274 passed; 0 failed; 1 ignored' in (OUT / 'perception-tests.log').read_text()
assert '2302 passed, 2 skipped' in (OUT / 'python-tests-final.log').read_text()
assert '9 passed' in (OUT / 'runner-tests.log').read_text()
assert 'Finished `release`' in (OUT / 'clippy.log').read_text()
assert not (OUT / 'format.log').read_bytes()
assert 'Refactor budget OK' in (OUT / 'refactor.log').read_text()
parity = json.loads((OUT / 'parity.json').read_text())
assert parity['baseline'] == parity['candidate']
for offset in [0, 5000]:
    integrity()
    args = [PYTHON, str(OUT / 'run-python-safe.py'), '--corpus', str(ROOT / 'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi'), '--offset', str(offset), '--limit', '1000', '--blocks', '20', '--baseline-wheel', str(next((OUT / 'baseline-wheel').glob('*.whl'))), '--candidate-wheel', str(next((OUT / 'candidate-wheel').glob('*.whl'))), '--baseline-path', str(OUT / 'baseline'), '--candidate-path', str(OUT / 'candidate'), '--base-revision', manifest['base_revision'], '--source-manifest', str(OUT / 'source-manifest.json'), '--output', str(OUT / f'frontier-{offset}.json')]
    with (OUT / f'frontier-{offset}.log').open('w') as log:
        subprocess.run(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    integrity()
    data = json.loads((OUT / f'frontier-{offset}.json').read_text())
    assert data['complete'] and len(data['blocks']) == 20
    assert all(data['blocks'][0][engine]['version'] == expected['version'] for engine, expected in comparators.items())
    print(offset, 'complete', flush=True)
