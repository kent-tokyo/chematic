from pathlib import Path
import subprocess,sys,os,json,hashlib,itertools,random
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017')
sys.path.insert(0,str(R/'scripts'));from bench_openbabel_steady_throughput import timing_row
sys.path.insert(0,str(R/'validation/cosmolkit_comparison'));from bench_corpus import paired_median_ratio_ci
os.environ.update(CHEMATIC_BENCH_WARMUP_NS='500000000',CHEMATIC_BENCH_TIMED_NS='100000000',CHEMATIC_BENCH_WARM_DUMP='1',CHEMATIC_BENCH_MIN_REPEATS='3')
binaries={'baseline':O/'baseline-sdf','candidate':O/'baseline-sdf','openbabel':O/'openbabel-sdf'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sealed={str(p):sha(p) for p in [*binaries.values(),O/'sdf-0.sdf']}
perms=list(itertools.permutations(binaries));random.Random(20261017).shuffle(perms)
for design,count in [('cyclic',21),('balanced',24)]:
 out=O/('same-binary-'+design+'.json');assert not out.exists(),'No diagnostic reruns'
 report={'purpose':'Same-binary control only, excluded from all candidate adoption statistics','design':design,'blocks':[],'binary_sha256':sealed,'same_rust_executable_for_both_labels':True,'complete':False}
 for index in range(count):
  names=list(binaries);order=names[index%3:]+names[:index%3] if design=='cyclic' else list(perms[index%6])
  row={'order':order,'load_average_start':list(os.getloadavg())}
  for arm in order:row[arm]={mode:timing_row(binaries[arm],O/'sdf-0.sdf',mode,3,1000) for mode in ['parse','write','roundtrip']}
  report['blocks'].append(row);out.write_text(json.dumps(report,indent=2)+'\n');print(design,index+1,count,flush=True)
 report['source_label_ratios']={mode:paired_median_ratio_ci([x['baseline'][mode]['elapsed_ns'] for x in report['blocks']],[x['candidate'][mode]['elapsed_ns'] for x in report['blocks']],seed=0) for mode in ['parse','write','roundtrip']}
 for name,digest in sealed.items():assert sha(Path(name))==digest
 report['complete']=True;out.write_text(json.dumps(report,indent=2)+'\n');print(design,report['source_label_ratios'],flush=True)
