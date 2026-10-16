from pathlib import Path
import gzip
import hashlib
import json
import math
import os
import struct
import subprocess

ROOT = Path('/private/tmp/chematic-aromatic-proof-20261016')
OUT = Path(__file__).resolve().parent
PYTHON = '/private/tmp/chematic-ring-frontier-env-20261016/bin/python'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def run(name, args, target=None):
    env = os.environ.copy()
    if target:
        env['PYTHONPATH'] = str(OUT / target)
    else:
        env.pop('PYTHONPATH', None)
    with (OUT / (name + '.log')).open('w') as stream:
        subprocess.run([PYTHON, *args], env=env, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=True)
    print(name, 'complete', flush=True)

manifest = json.loads((OUT / 'source-manifest.json').read_text())
assert all(sha((ROOT / name).read_bytes()) == expected for name, expected in manifest['candidate_source_sha256'].items())
proof = {}
for target in ['baseline', 'candidate']:
    proof[target] = {}
    for kind in ['descriptor', 'morgan', 'maccs']:
        path = OUT / f'{target}-{kind}.jsonl'
        args = [f'scripts/verify_python_{kind}_frontier.py', '--output', str(path)]
        if kind != 'maccs':
            args.extend(['--engine', 'chematic'])
        run(f'{target}-{kind}', args, target)
        raw = path.read_bytes()
        rows = [json.loads(line) for line in raw.splitlines()]
        outcomes = [row for row in rows if 'index' in row]
        assert len(outcomes) == 14072, (target, kind, len(outcomes))
        def errors(value):
            if isinstance(value, dict):
                return sum(int(k in ('error', 'parse_error')) + errors(v) for k, v in value.items())
            return 0
        failure_count = sum(errors(row) for row in outcomes)
        assert failure_count == 0, (target, kind, failure_count)
        def nonfinite(value):
            if isinstance(value, dict):
                return sum(int(k == 'float_bits' and not math.isfinite(struct.unpack('>d', bytes.fromhex(v))[0])) if k == 'float_bits' else nonfinite(v) for k, v in value.items())
            return 0
        invalid_floats = sum(nonfinite(row) for row in outcomes)
        proof[target][kind] = {'sha256': sha(raw), 'rows': len(outcomes), 'errors': failure_count, 'nonfinite_float_outcomes': invalid_floats}
        path.with_suffix(path.suffix + '.gz').write_bytes(gzip.compress(raw, mtime=0))
        path.unlink()
for kind in ['descriptor', 'morgan', 'maccs']:
    assert proof['baseline'][kind] == proof['candidate'][kind], kind
proof['counts'] = {'descriptor_outcomes': 112576, 'morgan_outputs': 56288, 'maccs_outcomes': 28144}
(OUT / 'parity.json').write_text(json.dumps(proof, indent=2) + '\n')
run('python-tests', ['-m', 'pytest', '-q', 'crates/chematic-py/tests'], 'candidate')
run('runner-tests', ['-m', 'pytest', '-q', 'scripts/tests/test_bench_python_speed_frontier.py', 'scripts/tests/test_cosmolkit_bench_corpus.py'])
assert all(sha((ROOT / name).read_bytes()) == expected for name, expected in manifest['candidate_source_sha256'].items())
print('All validation complete; no timings started', flush=True)
