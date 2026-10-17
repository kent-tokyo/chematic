"""Fixed balanced native experiment; preflight and formal epochs stay separate."""
from pathlib import Path
import argparse,datetime,gzip,hashlib,itertools,json,os,shutil,statistics,subprocess,sys
E=Path(__file__).resolve().parent
R=Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic')
sys.path.insert(0,str(R/'scripts'))
from bench_openbabel_steady_throughput import timing_row
from bench_openbabel_sdf_throughput import observe,semantic_check,eligible,paired_median_ratio_ci
from rdkit import Chem
os.environ.update(CHEMATIC_BENCH_WARMUP_NS='500000000',CHEMATIC_BENCH_TIMED_NS='100000000',CHEMATIC_BENCH_WARM_DUMP='1',CHEMATIC_BENCH_MIN_REPEATS='3',PYTHONDONTWRITEBYTECODE='1')
MODES=('parse','write','roundtrip')
PLAN=json.loads((E/'measurement-plan.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(raw):return hashlib.sha256(raw).hexdigest()
def atomic(path,value):
    raw=json.dumps(value,indent=2,allow_nan=False).encode()+b'\n';part=path.with_name(path.name+'.partial');part.write_bytes(raw);part.replace(path)
def integrity():
    seal=json.loads((E/'seal.json').read_text())
    for path,h in seal['files_sha256'].items():assert sha(Path(path))==h,path
    for path,h in seal['source_files_sha256'].items():assert sha(R/path)==h,path
    return sha(E/'seal.json')
def native(binary,corpus,mode,count):
    return subprocess.check_output([str(binary),str(corpus),mode,'1',str(count)],timeout=180)
def preflight():
    assert not (E/'preflight.json').exists() and not list(E.glob('paired-*.json'))
    records={}
    for profile in PLAN['profiles']:
        kind,name,count=profile['kind'],profile['profile'],profile['records'];key=kind+'-'+name;corpus=E/profile['path']
        assert sha(corpus)==profile['sha256']
        streams={};record={'records':count,'snapshot_sha256':{},'rewrite_sha256':{},'source_parse_snapshots_equal':True,'source_warmed_rewrites_equal':True,'instrumentation':{}}
        for arm in ('baseline','candidate'):
            for mode,field in [('snapshot','snapshot_sha256'),('dump','rewrite_sha256')]:
                raw=native(E/(arm+'-'+kind),corpus,mode,count);streams[arm,mode]=raw;record[field][arm]=digest(raw)
                (E/f'preflight-{key}-{arm}-{mode}.gz').write_bytes(gzip.compress(raw,mtime=0))
        for mode in ('snapshot','dump'):assert streams['baseline',mode]==streams['candidate',mode],(key,mode,'source difference')
        sem={}
        if kind=='sdf':
            reference=observe(corpus.read_bytes());assert len(reference)==count
            for arm in ('baseline','candidate','openbabel'):
                raw=streams[arm,'dump'] if arm!='openbabel' else native(E/'openbabel-sdf',corpus,'dump',count)
                if arm=='openbabel':
                    record['rewrite_sha256'][arm]=digest(raw);(E/f'preflight-{key}-{arm}-dump.gz').write_bytes(gzip.compress(raw,mtime=0))
                    smiles=native(E/'openbabel-sdf',corpus,'smiles',count);(E/f'preflight-{key}-{arm}-smiles.gz').write_bytes(gzip.compress(smiles,mtime=0))
                    text=[line.split()[0] for line in smiles.decode().splitlines()]
                else:text=[row['canonical_smiles'] for row in json.loads(streams[arm,'snapshot'])]
                canonical=[];refusals=[]
                for index,smiles in enumerate(text):
                    m=Chem.MolFromSmiles(smiles)
                    if m is None:refusals.append(index);canonical.append(None)
                    else:canonical.append(Chem.MolToSmiles(m))
                bad=[i for i,(a,b) in enumerate(zip(reference,canonical,strict=True)) if a['graph']!=b]
                sem[arm+'_parse']={'passed':not bad and not refusals,'rows':len(canonical),'mismatches':len(bad),'all_mismatch_indices':bad,'refusal_indices':refusals}
                sem[arm]=semantic_check(reference,observe(raw))
                atomic(E/f'preflight-{key}-{arm}-observed.json',{'input':reference,'output':observe(raw),'parsed_canonical':canonical})
            record['semantics']=sem
        # All actual modes/arms execute once before any formal epoch. Their samples are excluded.
        for arm in profile['orders'][0]:
            record['instrumentation'][arm]={mode:timing_row(E/(arm+'-'+kind),corpus,mode,3,count) for mode in MODES}
        records[key]=record;print('preflight',key,'passed',flush=True)
    atomic(E/'preflight.json',{'passed':True,'formal_observations':0,'preflight_timings_excluded':True,'profiles':records,'records':sum(x['records'] for x in PLAN['profiles'])})
def formal():
    pre=json.loads((E/'preflight.json').read_text());assert pre['passed'] and pre['formal_observations']==0
    validation=json.loads((E/'validation.json').read_text());assert validation['passed']
    assert json.loads((E/'code-realization.json').read_text())['source_contrast_valid']
    assert not list(E.glob('paired-*.json')) and not (E/'decision.json').exists(),'no timing reruns'
    sealed=integrity();source={};external={};floor=[];primary={};negative=[]
    for profile in PLAN['profiles']:
        assert integrity()==sealed
        kind,name,count=profile['kind'],profile['profile'],profile['records'];key=kind+'-'+name;corpus=E/profile['path'];output=E/f'paired-{key}.json';journal=E/(key+'-journals');journal.mkdir(exist_ok=False)
        report={'schema':'v2000-zero-fields/v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'corpus':profile,'seal_sha256':sealed,'preflight':pre['profiles'][key],'blocks':[],'summary':{},'complete':False}
        atomic(output,report)
        for index,order in enumerate(profile['orders']):
            assert shutil.disk_usage(E).free>=1024**3,'less than1GiB before block; stop without adoption'
            block={'order':order,'load_average_start':list(os.getloadavg())}
            for arm in order:
                block[arm]={}
                for mode in MODES:
                    row=timing_row(E/(arm+'-'+kind),corpus,mode,3,count)
                    atomic(journal/f'{index:02}-{arm}-{mode}.json',{'profile':key,'block':index,'arm':arm,'mode':mode,'timing':row})
                    block[arm][mode]=row
            report['blocks'].append(block);atomic(output,report);print(key,index+1,24,flush=True)
        for mode in MODES:
            values={arm:[block[arm][mode]['elapsed_ns'] for block in report['blocks']] for arm in profile['orders'][0]}
            ratio=paired_median_ratio_ci(values['baseline'],values['candidate'],seed=0);median=statistics.median(a/b for a,b in zip(values['baseline'],values['candidate'],strict=True));source[key+'/'+mode]=ratio
            report['summary'][mode]={'source':ratio,'source_unrounded_median':median,'openbabel':None}
            if median<.99:floor.append({'lane':key+'/'+mode,'unrounded_median':median,'speedup':ratio})
            if ratio['ci95'][1]<1:negative.append({'lane':key+'/'+mode,'speedup':ratio})
            if key+'/'+mode in PLAN['primary']:primary[key+'/'+mode]=ratio
            if 'openbabel' in values:
                allowed=eligible(report['preflight']['semantics'],True,mode,'openbabel');speed=paired_median_ratio_ci(values['openbabel'],values['candidate'],seed=1) if allowed else None
                ranked={'semantic_prerequisite':allowed,'speedup':speed,'strict_speed_win':bool(speed and speed['ci95'][0]>1 and speed['winning_blocks']==24)};report['summary'][mode]['openbabel']=ranked;external[key+'/'+mode]=ranked
        assert integrity()==sealed;report['complete']=True;atomic(output,report);print(key,'complete',flush=True)
    passed=len(primary)==2 and all(r['ci95'][0]>1 for r in primary.values()) and not floor
    atomic(E/'decision.json',{'experiment_complete':True,'local_adoption_gate_passed':passed,'primary':primary,'source_ratios':source,'source_floor_failures':floor,'negative_intervals':negative,'fresh_Open_Babel':external,'source_profiles':10,'rows':pre['records'],'all_source_outputs_exact':True,'universal_superiority':False})
    print('complete; adoption',passed,flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    preflight() if args.preflight_only else formal()
