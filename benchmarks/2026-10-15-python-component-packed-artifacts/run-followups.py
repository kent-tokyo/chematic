from pathlib import Path
import hashlib,json,subprocess
r=Path('/private/tmp/chematic-component-packed-20261015');s=Path('/private/tmp/chematic-component-packed-source-20261015');py='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';m=json.loads((r/'source-manifest.json').read_text());order=['tpsa','labute_asa','num_rings','morgan2_2048','morgan2_chiral','morgan3_2048','maccs']
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
print('All adaptive follow-ups complete',flush=True)
