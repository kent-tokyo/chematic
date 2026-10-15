from pathlib import Path
import json,hashlib,shutil,sys,subprocess,os,time,errno,re
r=Path('/private/tmp/chematic-openbabel-ascii-20261015');s=Path('/private/tmp/chematic-openbabel-ascii-source-20261015');py='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';sys.path.insert(0,str(s/'scripts'));import bench_openbabel_v3000_throughput as bench
output=r/'paired21-chains.json';report=json.loads(output.read_text());original=json.loads(output.read_text());assert not report['complete'] and len(report['blocks'])==14
for name in ['paired21-chains.json','paired21-chains.log','measurement.log','run_formal.log']:
 shutil.copyfile(r/name,r/(Path(name).stem+'-enospc'+Path(name).suffix))
p=output.with_suffix('.partial')
if p.exists():shutil.copyfile(p,r/'failed-checkpoint.partial')
source=report['source_manifest'];binaries={n:Path(v['path']) for n,v in report['binaries'].items()};corpus=Path(report['corpus']['path']);reserve=r/'checkpoint-reserve.bin'
with reserve.open('wb') as f:
 for _ in range(16):f.write(b'x'*1024*1024)
def integrity():
 assert json.loads((r/'source-manifest.json').read_text())==source
 assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in source['frozen_build_source_sha256'].items())
 assert all(hashlib.sha256(Path(n).read_bytes()).hexdigest()==h for n,h in source['native_libraries_sha256'].items())
 assert all(hashlib.sha256(p.read_bytes()).hexdigest()==report['binaries'][n]['sha256'] for n,p in binaries.items())
 assert hashlib.sha256(corpus.read_bytes()).hexdigest()==report['corpus']['sha256']
integrity();checks=bench.inspect(binaries,corpus,64,r/'resume-semantic.artifacts')
for k,v in checks.items():
 if k == 'rewrite_sha256':
  for arm in ['baseline','candidate']:assert report[k][arm]==v[arm]
 else:assert report[k]==v,k
old_ob=(r/'paired21-chains.artifacts/openbabel.molblocks').read_bytes().split(b'$$$$\n')
new_ob=(r/'resume-semantic.artifacts/openbabel.molblocks').read_bytes().split(b'$$$$\n')
assert len(old_ob)==len(new_ob)==65 and old_ob[-1]==new_ob[-1]==b''
for left,right in zip(old_ob[:-1],new_ob[:-1]):
 a=left.splitlines(keepends=True);b=right.splitlines(keepends=True)
 assert len(a)==len(b) and a[:1]==b[:1] and a[2:]==b[2:]
 assert re.fullmatch(rb' OpenBabel[0-9]{10}2D\n',a[1]) and re.fullmatch(rb' OpenBabel[0-9]{10}2D\n',b[1])
(r/'resume-semantic-verification.json').write_text(json.dumps({'semantics_and_source_snapshots_revalidated':True,'source_rewrite_bytes_identical':True,'openbabel_only_program_timestamps_differ':True,'openbabel_records_compared':64,'initial_rewrite_sha256':report['rewrite_sha256'],'resumed_rewrite_sha256':checks['rewrite_sha256']},indent=2)+'\n')
record={'reason':'Transient ENOSPC saving native block checkpoint','saved_blocks':14,'remaining_blocks':7,'unavailable_observation':'Attempted block15 finished all timed operations in memory, but checkpoint open failed. Those unpersisted timings were lost and are explicitly unavailable. Resume recollects block15; all fourteen saved original blocks remain unchanged.','protocol':'Original runner timing_row, observer, eligibility and deterministic bootstrap functions imported unchanged. Source/corpus/library/binary hashes and all semantic/output hashes revalidated before continuing identical rotating order. No performance-selected samples discarded, no pooled candidates.','checkpoint_recovery':'Own 16MiB temporary reserve is released on ENOSPC; the current computed report remains in memory until saved. No native timer or output semantics changed.'};(r/'measurement-recovery.json').write_text(json.dumps(record,indent=2)+'\n')
def save():
 data=json.dumps(report,indent=2)+'\n'
 for attempt in range(60):
  try:p.write_text(data);p.replace(output);return
  except OSError as e:
   if e.errno!=errno.ENOSPC:raise
   if reserve.exists():reserve.unlink();print('Released owned checkpoint reserve; current observations retained in memory',flush=True)
   if attempt==59:raise
   time.sleep(1)
names=list(binaries)
for block in range(14,21):
 order=names[block%3:]+names[:block%3];row={'order':order,'load_average_start':list(os.getloadavg())}
 for name in order:row[name]={mode:bench.timing_row(binaries[name],corpus,mode,5,64) for mode in bench.MODES}
 report['blocks'].append(row);save();print(f'block {block+1}/21',flush=True)
for mode in bench.MODES:
 report['summary'][mode]={}
 for index,name in enumerate(('baseline','openbabel')):
  allowed=bench.eligible(report['semantics'],report['source_parse_snapshots_equal'],mode,name);interval=bench.paired_median_ratio_ci([row[name][mode]['elapsed_ns'] for row in report['blocks']],[row['candidate'][mode]['elapsed_ns'] for row in report['blocks']],seed=index) if allowed else None
  report['summary'][mode][name]={'speedup':interval,'semantic_prerequisite':allowed,'strict_speed_win':bool(interval and interval['ci95'][0]>1 and interval['winning_blocks']==21)}
integrity();assert report['blocks'][:14]==original['blocks'];report['complete']=True;save()
if reserve.exists():reserve.unlink()
print('Native initial matrices completed; all 14 original saved blocks unchanged',flush=True)
with (r/'source_diagnostics.log').open('w') as log:subprocess.run(['python3',str(r/'source_diagnostics.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
print('All continuation and predeclared adaptive diagnostics completed',flush=True)
