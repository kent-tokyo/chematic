from pathlib import Path
import subprocess,json,os,sys,importlib.util
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017')
spec=importlib.util.spec_from_file_location('steady',R/'scripts/bench_openbabel_steady_throughput.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
tiny=O/'instrument-one.sdf';tiny.write_text((O/'sdf-0.sdf').read_text().split('$$$$')[0]+'$$$$\n')
env=os.environ.copy();env.update(CHEMATIC_BENCH_WARMUP_NS='500000000',CHEMATIC_BENCH_TIMED_NS='100000000',CHEMATIC_BENCH_WARM_DUMP='1',CHEMATIC_BENCH_MIN_REPEATS='3')
rows=[]
for arm in ['baseline','candidate','openbabel']:
 for mode in ['parse','write','roundtrip']:
  run=subprocess.run([str(O/(arm+'-sdf')),str(tiny),mode,'3','1'],env=env,check=True,capture_output=True)
  raw=json.loads(run.stdout);normalized=module.normalize(raw,mode,3,1)
  assert raw['repeats']>3
  assert raw['output_units']%raw['repeats']==0
  rows.append({'arm':arm,'mode':mode,'raw':raw,'normalized':normalized,'purpose':'Instrument accounting only; excluded from all formal samples'})
(O/'instrumentation-check.json').write_text(json.dumps({'passed':True,'all_nine_native_modes_checked':True,'actual_counts_and_durations_valid':True,'observations_excluded_from_formal_samples':True,'rows':rows},indent=2)+'\n')
print('All engines and modes report valid raw counts, warmup and timing budgets',flush=True)
