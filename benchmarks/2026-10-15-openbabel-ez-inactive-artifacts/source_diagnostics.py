from pathlib import Path
import sys,subprocess,json,hashlib,os,errno,time
r=Path('/private/tmp/chematic-openbabel-ez-gate-20261015');s=Path('/private/tmp/chematic-openbabel-ez-gate-source-20261015');old=Path('/private/tmp/chematic-openbabel-v3000-20261014');sys.path.insert(0,str(s/'validation/cosmolkit_comparison'))
from bench_corpus import paired_median_ratio_ci
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((r/'source-manifest.json').read_text());modes=['parse','write','roundtrip'];arms=['baseline','candidate']
reserve=r/'diagnostic-checkpoint-reserve.bin'
with reserve.open('wb') as f:
 for _ in range(16):f.write(b'x'*1024*1024)
def save_report(path, report):
 data=json.dumps(report,indent=2)+'\n';temporary=path.with_suffix('.partial')
 for attempt in range(60):
  try:temporary.write_text(data);temporary.replace(path);return
  except OSError as error:
   if error.errno!=errno.ENOSPC:raise
   if reserve.exists():reserve.unlink();print('Released owned diagnostic reserve; current observations retained in memory',flush=True)
   if attempt==59:raise
   time.sleep(1)

def verify():
 assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items())
 for arm in arms:assert sha(r/(arm+'-v3000'))==m['binaries'][arm]['sha256']

def equality_gate():
 artifacts=r/'continuation.artifacts';artifacts.mkdir(exist_ok=True)
 streams={}
 for profile,p in [('original',old/'corpus-chains.molblocks'),('continued',r/'corpus-continued-chains.molblocks')]:
  for arm in arms:
   for kind,mode in [('molblocks','dump'),('snapshot','snapshot')]:
    data=subprocess.run([str(r/(arm+'-v3000')),str(p),mode,'1','64'],capture_output=True,check=True,timeout=180).stdout
    path=artifacts/(profile+'-'+arm+'.'+kind);path.write_bytes(data);streams[profile,arm,kind]=data
 for kind in ['molblocks','snapshot']:
  ref=streams['original','baseline',kind]
  assert all(v==ref for (profile,arm,k),v in streams.items() if k==kind),kind
 return {'all_original_and_continued_source_snapshots_and_rewrites_byte_identical':True,'sha256':{p.name:sha(p) for p in artifacts.iterdir() if p.is_file()},'records':64}

def run_profile(label,corpus,records,gate):
 verify();expected={};report={'schema':'native-source-only-regression-diagnostic/v1','label':label,'source_manifest_sha256':sha(r/'source-manifest.json'),'binary_sha256':{arm:sha(r/(arm+'-v3000')) for arm in arms},'corpus':{'path':str(corpus),'sha256':sha(corpus),'records':records},'protocol':{'blocks':21,'repeats':5,'modes':modes,'order':'alternating two-source-arm fresh processes, sequential','timer':'unchanged native harness elapsed_ns; input framing, preparation and warm write excluded; whole-process startup excluded; RSS not measured','scope':'source-only diagnostic, never an Open Babel comparison','negative_policy':'one separate all-mode follow-up only when initial CI wholly below 1; no pooling or recursive retries'},'semantic_gate':gate,'blocks':[],'summary':{},'complete':False}
 def save():save_report(r/(label+'.json'),report)
 for block in range(21):
  order=arms[block%2:]+arms[:block%2];row={'order':order,'load_average_start':list(os.getloadavg())}
  for arm in order:
   row[arm]={}
   for mode in modes:
    out=subprocess.run([str(r/(arm+'-v3000')),str(corpus),mode,'5',str(records)],capture_output=True,text=True,check=True,timeout=180)
    v=json.loads(out.stdout);assert v['engine']=='chematic' and v['mode']==mode and v['records']==records*5 and v['repeats']==5 and v['elapsed_ns']>0
    if mode not in expected:expected[mode]=v['output_units']
    assert v['output_units']==expected[mode],(label,arm,mode)
    row[arm][mode]=v
  row['load_average_end']=list(os.getloadavg());report['blocks'].append(row);save()
 for mode in modes:
  stats=paired_median_ratio_ci([b['baseline'][mode]['elapsed_ns'] for b in report['blocks']],[b['candidate'][mode]['elapsed_ns'] for b in report['blocks']],seed=0)
  report['summary'][mode]={'speedup':stats,'strict_speed_win':stats['ci95'][0]>1 and stats['winning_blocks']==21}
 report['complete']=True;save();verify();print(label,report['summary'],flush=True);return report

if __name__=='__main__':
 gate=equality_gate();(r/'continuation-equivalence.json').write_text(json.dumps(gate,indent=2)+'\n')
 initial=run_profile('continued21',r/'corpus-continued-chains.molblocks',64,gate)
 selections={}
 for profile in ['0','5000','chains']:
  data=json.loads((r/f'paired21-{profile}.json').read_text());negative={mode:v['baseline']['speedup'] for mode,v in data['summary'].items() if v['baseline']['speedup'] and v['baseline']['speedup']['ci95'][1]<1}
  if negative:selections[profile]={'negative_initial_lanes':negative,'corpus':str(old/f'corpus-{profile}.molblocks'),'records':64 if profile=='chains' else 1000,'gate':{'source_snapshot_sha256':data['snapshot_sha256'],'byte_identical_rewrites':True,'source_parse_snapshots_equal':data['source_parse_snapshots_equal']}}
 negative={mode:v['speedup'] for mode,v in initial['summary'].items() if v['speedup']['ci95'][1]<1}
 if negative:selections['continued']={'negative_initial_lanes':negative,'corpus':str(r/'corpus-continued-chains.molblocks'),'records':64,'gate':gate}
 (r/'adaptive-selection.json').write_text(json.dumps({'policy':'One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.','profiles':selections},indent=2)+'\n')
 for profile,v in selections.items():run_profile('focused21-'+profile,Path(v['corpus']),v['records'],v['gate'])
 if reserve.exists():reserve.unlink()
 print('All source diagnostics complete',flush=True)
