from pathlib import Path
import gzip, hashlib, json, os, subprocess, sys, importlib
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
PYTHON='/private/tmp/chematic-ring-frontier-env-20261016/bin/python'
sys.path.insert(0,str(ROOT/'scripts'))
proof={'scope':'Fresh post-timing output observers only; all observer elapsed fields discarded, no samples added or pooled','profiles':{}}
for lane in ['frontier','bundle']:
 for offset in [0,5000]:
     bench=importlib.import_module('bench_python_scalar_bundle' if lane=='bundle' else 'bench_python_speed_frontier')
     data=json.loads((OUT/f'{lane}-{offset}.json').read_text())
     assert data['complete'] and len(data['blocks'])==20
     observers={}
     for arm in ['baseline','candidate','rdkit','cosmolkit']:
         env=os.environ.copy()
         if arm in ('baseline','candidate'):env['PYTHONPATH']=str(OUT/arm)
         else:env.pop('PYTHONPATH',None)
         args=[PYTHON,str(ROOT/('scripts/bench_python_scalar_bundle.py' if lane=='bundle' else 'scripts/bench_python_speed_frontier.py')),'--worker',arm,'--corpus',str(ROOT/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi'),'--offset',str(offset),'--limit','1000','--ops',','.join(bench.OPS)]
         p=subprocess.run(args,env=env,cwd=ROOT,capture_output=True,text=True,check=True)
         sample=json.loads(p.stdout)
         values={op:record['values'] for op,record in sample['operations'].items()}
         output_hash=bench.sha(json.dumps(values,sort_keys=True).encode())
         assert all(block[arm]['output_sha256']==output_hash and block[arm]['extension_sha256']==sample['extension_sha256'] and block[arm]['version']==sample['version'] for block in data['blocks'])
         observers[arm]=values
         raw=json.dumps({'scope':proof['scope'],'engine':arm,'offset':offset,'version':sample['version'],'output_sha256':output_hash,'extension_sha256':sample['extension_sha256'],'values':values},sort_keys=True).encode()
         path=OUT/f'observer-{lane}-{offset}-{arm}.json.gz';path.write_bytes(gzip.compress(raw,mtime=0))
     checks={}
     differences={}
     for op in bench.OPS:
         row={'candidate_vs_baseline':bench.agreement(observers['baseline'][op],observers['candidate'][op],exact=True)}
         row.update({name+'_vs_rdkit':bench.agreement(observers['rdkit'][op],observers[name][op],exact=False) for name in ['baseline','candidate','cosmolkit']})
         assert row==data['agreement'][op]
         checks[op]=row
         differences[op]=[{'row':i,'rdkit':a,'cosmolkit':b} for i,(a,b) in enumerate(zip(observers['rdkit'][op],observers['cosmolkit'][op])) if not bench.agreement([a],[b],exact=False)['passed']]
     proof['profiles'][f'{lane}-{offset}']={'all_20_output_and_extension_hashes_revalidated':True,'agreement':checks,'all_cosmolkit_differences':differences}
     print(offset,'all observer outputs revalidated',flush=True)
(OUT/'observer-verification.json').write_text(json.dumps(proof,indent=2)+'\n')
