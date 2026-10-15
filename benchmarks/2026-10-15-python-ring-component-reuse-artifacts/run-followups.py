from pathlib import Path
import hashlib,json,subprocess
r=Path('/private/tmp/chematic-ring-component-20261015');s=Path('/private/tmp/chematic-ring-component-source-20261015');py='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';m=json.loads((r/'source-manifest.json').read_text());order=['tpsa','labute_asa','num_rings','morgan2_2048','morgan2_chiral','morgan3_2048','maccs']
def verify():
 assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in m['candidate_sha256'].items())
def run(offset,ops,label,base=None,manifest=None,revision=None):
 verify();base=base or r
 cmd=[py,str(s/'scripts/bench_python_speed_frontier.py'),'--corpus',str(s/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi'),'--offset',str(offset),'--limit','1000','--blocks','20','--ops',','.join(ops),'--baseline-wheel',str(next((base/'baseline-wheel').glob('*.whl'))),'--baseline-path',str(base/'baseline'),'--candidate-wheel',str(next((r/'candidate-wheel').glob('*.whl'))),'--candidate-path',str(r/'candidate'),'--source-manifest',str(manifest or r/'source-manifest.json'),'--base-revision',revision or 'immediate-source-candidate-baseline','--output',str(r/(label+'.json'))]
 with (r/(label+'.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 verify();print(label,'complete',flush=True)
selection={}
for offset in [5000,0]:
 j=json.loads((r/f'frontier-{offset}.json').read_text());negative={k:v['baseline']['speedup'] for k,v in j['summary'].items() if v['baseline'].get('speedup') and v['baseline']['speedup']['ci95'][1]<1}
 ops=[op for op in order if any(k.startswith(op+'/') for k in negative)]
 selection[str(offset)]={'negative_initial_lanes':negative,'operations':ops,'policy':'One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate.'}
(r/'adaptive-selection.json').write_text(json.dumps(selection,indent=2)+'\n')
for offset in [5000,0]:
 ops=selection[str(offset)]['operations']
 if ops:run(offset,ops,f'focused-{offset}')
old=Path('/private/tmp/chematic-kekule-view-20261014');oldmanifest=old/'source-manifest-final.json';olddata=json.loads(oldmanifest.read_text())
assert hashlib.sha256(next((old/'baseline-wheel').glob('*.whl')).read_bytes()).hexdigest()==olddata['wheels']['baseline']
(r/'historical-baseline-source-manifest.json').write_bytes(oldmanifest.read_bytes())
audit={'plan':json.loads((r/'historical-audit-plan.json').read_text()),'historical_manifest_sha256':hashlib.sha256(oldmanifest.read_bytes()).hexdigest(),'candidate_manifest_sha256':hashlib.sha256((r/'source-manifest.json').read_bytes()).hexdigest(),'note':'Comparisons are cumulative source candidates, not isolated attribution or a verified main revision.'}
(r/'historical-audit-manifest.json').write_text(json.dumps(audit,indent=2)+'\n')
run(5000,['num_rings'],'historical-rings-5000',base=old,manifest=r/'historical-audit-manifest.json',revision='cac480fb23cd8adeec6e5e9bed1f18677f87b392-plus-cycle-and-writer-before-kekule-views')
print('All follow-ups complete',flush=True)
