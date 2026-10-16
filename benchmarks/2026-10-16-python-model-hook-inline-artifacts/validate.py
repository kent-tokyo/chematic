from pathlib import Path
import gzip, hashlib, json, math, os, struct, subprocess, zipfile
ROOT = Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT = Path(__file__).resolve().parent
PYTHON = '/private/tmp/chematic-ring-frontier-env-20261016/bin/python'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def run(name, args, arm):
    env = os.environ.copy()
    env['PYTHONPATH'] = str(OUT / arm)
    env['PYTEST_DISABLE_PLUGIN_AUTOLOAD'] = '1'
    with (OUT / (name + '.log')).open('w') as log:
        subprocess.run([PYTHON, *args], env=env, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(name, 'complete', flush=True)

proof = {'wheels': {}, 'outputs': {}}
for arm in ['baseline', 'candidate']:
    wheel = next((OUT / (arm + '-wheel')).glob('*.whl'))
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(OUT / arm)
        extensions = {str(Path(name)): sha(archive.read(name)) for name in archive.namelist() if name.endswith('.so')}
    assert all(sha((OUT / arm / name).read_bytes()) == expected for name, expected in extensions.items())
    proof['wheels'][arm] = {'wheel': str(wheel), 'sha256': sha(wheel.read_bytes()), 'extensions': extensions}
    proof['outputs'][arm] = {}
    for kind in ['descriptor', 'morgan', 'maccs']:
        path = OUT / f'{arm}-{kind}.jsonl'
        script = OUT / 'verify_maccs.py' if kind == 'maccs' else ROOT / f'scripts/verify_python_{kind}_frontier.py'
        args = [str(script), '--output', str(path)]
        if kind != 'maccs':
            args += ['--engine', 'chematic']
        run(f'{arm}-{kind}', args, arm)
        raw = path.read_bytes()
        outcomes = [json.loads(line) for line in raw.splitlines() if '"index"' in line.decode()]
        assert len(outcomes) == 14072
        def errors(value):
            if isinstance(value, dict):
                return sum(int(k in ('error', 'parse_error')) + errors(v) for k, v in value.items())
            return 0
        assert sum(errors(row) for row in outcomes) == 0
        def nonfinite(value):
            if isinstance(value, dict):
                return sum(int(not math.isfinite(struct.unpack('>d', bytes.fromhex(v))[0])) if k == 'float_bits' else nonfinite(v) for k, v in value.items())
            return 0
        proof['outputs'][arm][kind] = {'sha256': sha(raw), 'rows': len(outcomes), 'errors': 0, 'nonfinite_float_outcomes': sum(nonfinite(row) for row in outcomes)}
        path.with_suffix(path.suffix + '.gz').write_bytes(gzip.compress(raw, mtime=0))
        path.unlink()
assert proof['outputs']['baseline'] == proof['outputs']['candidate']
proof['counts'] = {'descriptor_outcomes': 112576, 'morgan_outputs': 56288, 'maccs_outcomes': 28144}
(OUT / 'parity.json').write_text(json.dumps(proof, indent=2) + '\n')
run('python-tests', ['-m', 'pytest', '-q', '-p', 'no:cacheprovider', 'crates/chematic-py/tests'], 'candidate')
print('All Python validations passed; no formal timings started', flush=True)
