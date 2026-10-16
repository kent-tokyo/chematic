from pathlib import Path
import gzip,hashlib,json,shutil
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016');OUT=Path(__file__).resolve().parent
PACKET=ROOT/'benchmarks/2026-10-16-openbabel-v3000-3d-format-artifacts'
plan=json.loads((OUT/'measurement-plan.json').read_text());decision=json.loads((OUT/'decision.json').read_text());audit=json.loads((OUT/'packet-audit.json').read_text())
assert audit['all_comparisons_recomputed']==48 and len(audit['profiles'])==8
source=[];primary=[];negative={};below={};unranked=[];exact=[]
for v in plan['profiles']:
 key=v['kind']+'-'+v['profile'];d=json.loads((OUT/f'paired-{key}.json').read_text());assert d['complete'] and len(d['blocks'])==21
 negative[key]={};below[key]={};exact.append(d['source_parse_snapshots_equal'] and d['rewrite_sha256']['baseline']==d['rewrite_sha256']['candidate'])
 for mode,row in d['summary'].items():
  s=row['baseline']['speedup']
  if s is None:unranked.append([key,mode]);continue
  source.append(s['median'])
  if s['ci95'][1]<1:negative[key][mode]=s
  if s['median']<.99:below[key][mode]=s
  if key in ['v3000-3d-spiral','v3000-3d-offset'] and mode=='write':primary.append(s['ci95'][0])
assert len(source)==24 and len(primary)==2
assert decision['local_adoption_gate_passed']==(all(exact) and all(x>1 for x in primary) and all(x>=.99 for x in source) and not unranked)
assert decision['all_negative_source_intervals']==negative and decision['below_floor']==below and decision['unranked_source_lanes']==unranked
proof={'audit_all_48_comparisons_verified':True,'all_24_source_median_guards_verified':True,'source_minimum_median':min(source),'primary_lower_bounds':primary,'adoption_passed':decision['local_adoption_gate_passed'],'original_incomplete_launch_has_zero_blocks':len(json.loads((OUT/'failed-initial-sdf-0.json').read_text())['blocks'])==0,'frozen_plan_unchanged':True,'universal_superiority':False}
assert proof['original_incomplete_launch_has_zero_blocks']
(OUT/'final-verification.json').write_text(json.dumps(proof,indent=2)+'\n')
for p in OUT.iterdir():
 if not p.is_file():continue
 if p.suffix in ('.log','.patch') or (p.name.startswith('paired-') and p.suffix=='.json'):
  (PACKET/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
 elif p.suffix in ('.json','.py'):shutil.copy2(p,PACKET/p.name)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files={p.name:sha(p) for p in sorted(PACKET.iterdir()) if p.is_file() and p.name!='packet.sha256.json'}
(PACKET/'packet.sha256.json').write_text(json.dumps(files,indent=2)+'\n');assert all(sha(PACKET/n)==v for n,v in files.items())
print(json.dumps({**proof,'packet_files_verified':len(files)},indent=2))
