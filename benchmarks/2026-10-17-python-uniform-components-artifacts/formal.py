from pathlib import Path
import hashlib,json,subprocess,sys,os
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-python-asa-preparation-20261017')
plan=json.loads((O/'measurement-plan.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text());comparators=json.loads((O/'comparator-manifest.json').read_text());proof=json.loads((O/'binary-source-proof.json').read_text());seal=json.loads((O/'protocol-seal.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def integrity():
 for p,h in manifest['candidate_sha256'].items():assert sha(R/p)==h,p
 for p,h in manifest['baseline_sha256'].items():assert sha(O/'baseline-source'/p)==h,p
 for item in comparators.values():
  for p,h in item['native_files'].items():assert sha(Path(p))==h,p
 for arm,item in proof['wheels'].items():
  assert sha(Path(item['path']))==item['sha256']
  for p,h in item['extensions'].items():assert sha(O/arm/p)==h,p
 for p,h in seal['files_sha256'].items():assert sha(O/p)==h,p
assert json.loads((O/'native-validation.json').read_text())['passed']
assert json.loads((O/'parity.json').read_text())['passed']
assert json.loads((O/'instrumentation-check.json').read_text())['passed']
assert proof['source_contrast_valid'] and proof['complete_frozen_function_exact']
assert not list(O.glob('frontier-*.json')), 'No valid timing reruns'
for offset in plan['offsets']:
 integrity();cmd=[sys.executable,str(O/'bench-fixed-batches.py'),'--corpus',str(R/plan['corpus']),'--offset',str(offset),'--limit',str(plan['unique_records_each']),'--blocks','20','--base-revision',plan['base_commit'],'--source-manifest',str(O/'source-manifest.json'),'--output',str(O/f'frontier-{offset}.json')]
 for arm in ['baseline','candidate']:cmd+=['--'+arm+'-wheel',proof['wheels'][arm]['path'],'--'+arm+'-path',str(O/arm)]
 env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
 with (O/f'frontier-{offset}.log').open('w') as log:subprocess.run(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 integrity();data=json.loads((O/f'frontier-{offset}.json').read_text());assert data['complete'] and len(data['blocks'])==20
 for engine in ['rdkit','cosmolkit']:assert data['blocks'][0][engine]['version']==comparators[engine]['runtime_version']
 assert data['protocol']['scalar_batch_copies']==20
 print(offset,'all 21 operation/mode lanes complete',flush=True)
sys.path.insert(0,str(R/'validation/cosmolkit_comparison'));from bench_corpus import paired_median_ratio_ci
ratios={};primary={};negative=[];floors=[]
for offset in plan['offsets']:
 data=json.loads((O/f'frontier-{offset}.json').read_text())
 for op in plan['operations']:
  assert data['agreement'][op]['candidate_vs_baseline']['passed'],op
  for mode in plan['modes']:
   base=[b['baseline']['operations'][op][mode] for b in data['blocks']];cand=[b['candidate']['operations'][op][mode] for b in data['blocks']]
   ratio=paired_median_ratio_ci(base,cand,seed=0);key=f'{offset}/{op}/{mode}';ratios[key]=ratio
   assert ratio==data['summary'][op+'/'+mode]['baseline']['speedup']
   if ratio['median']<.99:floors.append({'lane':key,'speedup':ratio})
   if ratio['ci95'][1]<1:negative.append({'lane':key,'speedup':ratio})
   if op=='labute_asa' and mode=='first_use':primary[str(offset)]=ratio
passed=all(r['ci95'][0]>1 for r in primary.values()) and len(primary)==2 and not floors
(O/'decision.json').write_text(json.dumps({'local_adoption_gate_passed':passed,'primary':primary,'source_ratios':ratios,'source_floor_failures':floors,'negative_intervals':negative,'all_14072_source_outcomes_exact':True,'source_valid_epochs_each_slice':1,'universal_superiority':False},indent=2)+'\n')
print(json.dumps({'adoption':passed,'primary':primary,'floor_failures':floors,'negative_intervals':negative}),flush=True)
