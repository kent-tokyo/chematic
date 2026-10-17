from pathlib import Path
import json,gzip,hashlib,subprocess,sys,os,statistics
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-openbabel-writer-20261017')
sys.path.insert(0,str(R/'scripts'))
from bench_openbabel_sdf_throughput import paired_median_ratio_ci
manifest=json.loads((O/'source-manifest.json').read_text()); seal=json.loads((O/'binary-seal.json').read_text()); plan=json.loads((O/'measurement-plan.json').read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def integrity():
 for name,expected in manifest['source_files_sha256'].items(): assert sha(R/name)==expected,name
 for name,expected in manifest['baseline_source_files_sha256'].items(): assert sha(O/'baseline-source'/name)==expected,name
 for name,expected in manifest['native_libraries_sha256'].items(): assert sha(Path(name))==expected,name
 for name,expected in seal['files_sha256'].items(): assert sha(O/name)==expected,name
assert json.loads((O/'validation.json').read_text())['passed']
assert json.loads((O/'binary-source-proof.json').read_text())['source_contrast_valid']
assert not list(O.glob('paired-*.json')), 'No valid timing reruns'
integrity()
checks={}
for profile in plan['profiles']:
 kind,name,count=profile['kind'],profile['profile'],profile['records']; corpus=O/profile['path'];key=kind+'-'+name
 results={}
 for arm in ['baseline','candidate']:
  for mode in ['snapshot','dump']:
   raw=subprocess.check_output([str(O/(arm+'-'+kind)),str(corpus),mode,'1',str(count)])
   (O/f'preflight-{key}-{arm}-{mode}.gz').write_bytes(gzip.compress(raw,mtime=0));results[arm,mode]=raw
 for mode in ['snapshot','dump']: assert results['baseline',mode]==results['candidate',mode],(key,mode)
 checks[key]={'records':count,'snapshots_equal':True,'rewrites_equal':True,'snapshot_sha256':hashlib.sha256(results['baseline','snapshot']).hexdigest(),'rewrite_sha256':hashlib.sha256(results['baseline','dump']).hexdigest()}
(O/'preflight-source-equality.json').write_text(json.dumps(checks,indent=2)+'\n')
print('All 4384 source snapshots and rewrites exact',flush=True)
for p in plan['profiles']:
 integrity();kind,name,count=p['kind'],p['profile'],p['records'];key=kind+'-'+name
 cmd=[sys.executable,str(R/'scripts'/('bench_openbabel_sdf_throughput.py' if kind=='sdf' else 'bench_openbabel_v3000_throughput.py')),'--baseline',str(O/('baseline-'+kind)),'--candidate',str(O/('candidate-'+kind)),'--openbabel',str(O/('openbabel-sdf' if kind=='sdf' else 'openbabel-v3000')),'--corpus',str(O/p['path']),'--records',str(count),'--blocks','21','--repeats','3','--source-manifest',str(O/'source-manifest.json'),'--output',str(O/('paired-'+key+'.json'))]
 with (O/('paired-'+key+'.log')).open('w') as log: subprocess.run(cmd,cwd=R,check=True,stdout=log,stderr=subprocess.STDOUT)
 integrity();result=json.loads((O/('paired-'+key+'.json')).read_text());assert result['complete'] and len(result['blocks'])==21
 assert result['source_parse_snapshots_equal'] and result['rewrite_sha256']['baseline']==result['rewrite_sha256']['candidate']
 print(key,'complete',flush=True)
source={};negative=[];failures=[];primary={}
for p in plan['profiles']:
 key=p['kind']+'-'+p['profile']; result=json.loads((O/('paired-'+key+'.json')).read_text())
 for mode in ['parse','write','roundtrip']:
  a=[b['baseline'][mode]['elapsed_ns'] for b in result['blocks']];b=[x['candidate'][mode]['elapsed_ns'] for x in result['blocks']]
  ratio=paired_median_ratio_ci(a,b,seed=0);source[key+'/'+mode]=ratio
  if ratio['ci95'][1]<1: negative.append({'lane':key+'/'+mode,'speedup':ratio})
  raw_median=statistics.median(x/y for x,y in zip(a,b,strict=True))
  if raw_median<.99: failures.append({'lane':key+'/'+mode,'speedup':ratio,'unrounded_median':raw_median})
  if key in plan['primary'] and mode=='write': primary[key]=ratio
passed=all(x['ci95'][0]>1 for x in primary.values()) and len(primary)==2 and not failures
decision={'local_adoption_gate_passed':passed,'primary':primary,'source_ratios':source,'negative_intervals':negative,'source_floor_failures':failures,'all_source_snapshots_and_rewrites_exact':True,'source_valid_epochs':10,'universal_superiority':False}
(O/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
print(json.dumps({'adoption':passed,'primary':primary,'floor_failures':failures,'negative_intervals':negative}),flush=True)
