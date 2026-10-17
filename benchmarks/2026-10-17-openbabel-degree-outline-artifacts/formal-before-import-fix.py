from pathlib import Path
import json,gzip,hashlib,subprocess,sys,os,statistics,datetime
O=Path(__file__).resolve().parent;C=Path('/private/tmp/chematic-degree-only-steady-evidence-20261017');R=Path('/private/tmp/chematic-degree-only-steady-20261017')
sys.path.insert(0,str(R/'scripts'))
from bench_openbabel_steady_throughput import timing_row
from bench_openbabel_sdf_throughput import paired_median_ratio_ci,observe,semantic_check,eligible,Chem
os.environ.update(CHEMATIC_BENCH_WARMUP_NS='500000000',CHEMATIC_BENCH_TIMED_NS='100000000',CHEMATIC_BENCH_WARM_DUMP='1',CHEMATIC_BENCH_MIN_REPEATS='3')
reference=json.loads((O/'semantic-reference.json').read_text())
plan=json.loads((O/'measurement-plan.json').read_text());manifest=json.loads((O/'source-manifest.json').read_text());seal=json.loads((O/'binary-seal.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(raw):return hashlib.sha256(raw).hexdigest()
def integrity():
 for name,h in manifest['source_files_sha256'].items():assert sha(R/name)==h,name
 for name,h in manifest['baseline_source_files_sha256'].items():assert sha(Path(manifest['baseline_source_root'])/name)==h,name
 for name,h in manifest['native_libraries_sha256'].items():assert sha(Path(name))==h,name
 for name,h in seal['files_sha256'].items():assert sha(O/name)==h,name
 for name,h in reference['prior_packet_files_sha256'].items():assert sha(R/name)==h,name
assert json.loads((O/'validation.json').read_text())['passed']
assert json.loads((O/'binary-source-proof.json').read_text())['source_contrast_valid']
assert not list(O.glob('paired-*.json')),'No timing reruns'
integrity();preflight={}
for profile in plan['profiles']:
 kind,name,count=profile['kind'],profile['profile'],profile['records'];key=kind+'-'+name;corpus=O/profile['path'];old=reference['assessments'][key]
 assert sha(corpus)==profile['sha256']==old['corpus']['sha256']
 streams={};proof={'records':count,'snapshot_sha256':{},'rewrite_sha256':{},'source_parse_snapshots_equal':True,'source_warmed_rewrites_equal':True,'semantics':old['semantics'],'semantic_basis':'Fresh byte-exact match to C source outcomes; C observer assessments retained. SDF assessments recomputed below.'}
 for arm in ['baseline','candidate']:
  for mode,field in [('snapshot','snapshot_sha256'),('dump','rewrite_sha256')]:
   raw=subprocess.check_output([str(O/(arm+'-'+kind)),str(corpus),mode,'1',str(count)])
   streams[arm,mode]=raw;proof[field][arm]=digest(raw)
   assert digest(raw)==old[field]['baseline'],(key,arm,mode,'C outcome changed')
   (O/f'preflight-{key}-{arm}-{mode}.gz').write_bytes(gzip.compress(raw,mtime=0))
 assert streams['baseline','snapshot']==streams['candidate','snapshot'] and streams['baseline','dump']==streams['candidate','dump']
 if kind=='sdf':
  reference=observe(corpus.read_bytes());assert len(reference)==count
  sem={}
  for arm in ['baseline','candidate','openbabel']:
   raw=streams[arm,'dump'] if arm!='openbabel' else subprocess.check_output([str(O/'openbabel-sdf'),str(corpus),'dump','1',str(count)])
   proof['rewrite_sha256'][arm]=digest(raw);assert digest(raw)==old['rewrite_sha256'][arm]
   sem[arm]=semantic_check(reference,observe(raw))
   if arm!='openbabel':
    parsed=json.loads(streams[arm,'snapshot']);canonical=[Chem.MolToSmiles(Chem.MolFromSmiles(row['canonical_smiles'])) for row in parsed]
   else:
    (O/f'preflight-{key}-openbabel-dump.gz').write_bytes(gzip.compress(raw,mtime=0))
    smiles=subprocess.check_output([str(O/'openbabel-sdf'),str(corpus),'smiles','1',str(count)]);(O/f'preflight-{key}-openbabel-smiles.gz').write_bytes(gzip.compress(smiles,mtime=0))
    canonical=[Chem.MolToSmiles(Chem.MolFromSmiles(row.split()[0])) for row in smiles.decode().splitlines()]
   bad=[i for i,(a,b) in enumerate(zip(reference,canonical,strict=True)) if a['graph']!=b]
   sem[arm+'_parse']={'passed':not bad,'rows':len(canonical),'mismatches':len(bad),'examples':bad[:8]}
  assert sem==old['semantics'];proof['semantics']=sem;proof['semantic_basis']='Fresh SDF observer assessment; all rewritten streams also match C hashes'
 preflight[key]=proof
(O/'preflight-source-equality.json').write_text(json.dumps(preflight,indent=2)+'\n');print('All 4384 fresh source snapshots and warmed rewrites exact; SDF semantics recomputed',flush=True)
source={};external={};negative=[];failures=[];primary={}
for profile in plan['profiles']:
 integrity();kind,name,count=profile['kind'],profile['profile'],profile['records'];key=kind+'-'+name;corpus=O/profile['path'];output=O/('paired-'+key+'.json')
 report={'schema':'balanced-steady-throughput/v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'corpus':profile,'source_manifest_sha256':sha(O/'source-manifest.json'),'preflight':preflight[key],'blocks':[],'summary':{},'complete':False}
 def save():output.write_text(json.dumps(report,indent=2)+'\n')
 save()
 for index,order in enumerate(profile['orders']):
  row={'order':order,'load_average_start':list(os.getloadavg())}
  for arm in order:row[arm]={mode:timing_row(O/(arm+'-'+kind),corpus,mode,3,count) for mode in ['parse','write','roundtrip']}
  report['blocks'].append(row);save();print(key,index+1,24,flush=True)
 for mode in ['parse','write','roundtrip']:
  values={arm:[b[arm][mode]['elapsed_ns'] for b in report['blocks']] for arm in profile['orders'][0]}
  ratio=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0);source[key+'/'+mode]=ratio
  median=statistics.median(a/b for a,b in zip(values['baseline'],values['candidate'],strict=True))
  report['summary'][mode]={'source':ratio,'source_unrounded_median':median,'openbabel':None}
  if median<.99:failures.append({'lane':key+'/'+mode,'speedup':ratio,'unrounded_median':median})
  if ratio['ci95'][1]<1:negative.append({'lane':key+'/'+mode,'speedup':ratio})
  if key+'/'+mode in plan['primary']:primary[key+'/'+mode]=ratio
  if 'openbabel' in values:
   allowed=eligible(report['preflight']['semantics'],True,mode,'openbabel');speed=paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1) if allowed else None
   report['summary'][mode]['openbabel']={'semantic_prerequisite':allowed,'speedup':speed,'strict_speed_win':bool(speed and speed['ci95'][0]>1 and speed['winning_blocks']==24)}
   external[key+'/'+mode]=report['summary'][mode]['openbabel']
 integrity();report['complete']=True;save();print(key,'complete',flush=True)
passed=len(primary)==2 and all(x['ci95'][0]>1 for x in primary.values()) and not failures
result={'local_adoption_gate_passed':passed,'primary':primary,'source_ratios':source,'negative_intervals':negative,'source_floor_failures':failures,'fresh_SDF_Open_Babel':external,'all_source_snapshots_and_rewrites_exact':True,'source_valid_epochs':10,'universal_superiority':False}
(O/'decision.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'adoption':passed,'primary':primary,'floor_failures':failures,'negative_intervals':negative}),flush=True)
