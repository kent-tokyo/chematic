from pathlib import Path
import gzip
import hashlib
import json
import subprocess
import sys

ROOT = Path('/private/tmp/chematic-3d-stereo-guard-20261017')
EVIDENCE = Path(__file__).resolve().parent
manifest_path = EVIDENCE / 'source-manifest.json'
manifest_bytes = manifest_path.read_bytes()
manifest = json.loads(manifest_bytes)
seal = json.loads((EVIDENCE / 'binary-seal.json').read_text())
profiles = [(v['kind'],v['profile'],v['records']) for v in json.loads((EVIDENCE/'measurement-plan.json').read_text())['profiles']]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def integrity():
    assert manifest_path.read_bytes() == manifest_bytes
    for n, expected in manifest['source_files_sha256'].items():
        assert sha((ROOT / n).read_bytes()) == expected, n
    for n, expected in manifest['native_libraries_sha256'].items():
        assert sha(Path(n).read_bytes()) == expected, n
    for n, expected in seal['files_sha256'].items():
        assert sha((EVIDENCE / n).read_bytes()) == expected, n

build_proof=json.loads((EVIDENCE/'binary-source-proof.json').read_text())
assert build_proof['source_contrast_valid'] and build_proof['compiler_rebuilt_production_library']
assert build_proof['candidate_commit']==manifest['candidate_commit']
assert not list(EVIDENCE.glob('paired-*.json')), 'Valid epoch already exists; no timing rerun'
integrity()
validation = json.loads((EVIDENCE / 'validation.json').read_text())
assert validation['passed']
source_checks = {}
for kind, profile, count in profiles:
    suffix = 'sdf' if kind == 'sdf' else 'molblocks'
    corpus = EVIDENCE / f'{kind}-{profile}.{suffix}'
    outputs = {}
    for engine in ['baseline','candidate']:
        for mode in ['snapshot','dump']:
            raw = subprocess.check_output([str(EVIDENCE / f'{engine}-{kind}'), str(corpus), mode, '1', str(count)])
            outputs[engine, mode] = raw
            (EVIDENCE / f'preflight-{kind}-{profile}-{engine}-{mode}.gz').write_bytes(gzip.compress(raw, mtime=0))
    for mode in ['snapshot','dump']:
        assert outputs['baseline',mode] == outputs['candidate',mode], (kind,profile,mode)
    source_checks[f'{kind}-{profile}'] = {'records':count, 'snapshots_equal':True, 'rewrites_equal':True,
        'snapshot_sha256':sha(outputs['baseline','snapshot']), 'rewrite_sha256':sha(outputs['baseline','dump'])}
(EVIDENCE / 'preflight-source-equality.json').write_text(json.dumps(source_checks,indent=2)+'\n')
for kind, profile, count in profiles:
    integrity()
    suffix = 'sdf' if kind == 'sdf' else 'molblocks'
    corpus = EVIDENCE / f'{kind}-{profile}.{suffix}'
    output = EVIDENCE / f'paired-{kind}-{profile}.json'
    command = [sys.executable,str(EVIDENCE/'run-native-safe.py'),('bench_openbabel_sdf_throughput.py' if kind=='sdf' else 'bench_openbabel_v3000_throughput.py'),
        '--baseline',str(EVIDENCE/f'baseline-{kind}'), '--candidate',str(EVIDENCE/f'candidate-{kind}'),
        '--openbabel',str(EVIDENCE/('openbabel-sdf' if kind=='sdf' else 'openbabel-v3000')),'--corpus',str(corpus), '--records',str(count),
        '--blocks','21','--repeats','3','--source-manifest',str(manifest_path),'--output',str(output)]
    with (EVIDENCE/f'paired-{kind}-{profile}.log').open('wb') as log:
        subprocess.run(command,check=True,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    integrity()
    result = json.loads(output.read_text())
    assert result['complete'] and len(result['blocks']) == 21
    assert result['source_parse_snapshots_equal']
    assert result['rewrite_sha256']['baseline'] == result['rewrite_sha256']['candidate']
    print(kind,profile,'complete',flush=True)
