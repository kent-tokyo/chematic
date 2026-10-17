from pathlib import Path
import os,subprocess,hashlib,json,math,struct
R=Path('/private/tmp/chematic-python-asa-preparation-20261017');O=Path(__file__).resolve().parent;P='/private/tmp/chematic-python-asa-env-20261017/bin/python'
s=(R/'benchmarks/2026-10-16-python-scalar-model-reuse-artifacts/verify_maccs.py').read_text().replace('/private/tmp/chematic-sdf-buffer-scan-20261016',str(R));(O/'verify_maccs.py').write_text(s)
proof={'outputs':{},'python_tests':{}}
for arm in ['baseline','candidate']:
 env=os.environ.copy();env.update(PYTHONPATH=str(O/arm),PYTHONDONTWRITEBYTECODE='1')
 for kind,script,args in [('descriptor',R/'scripts/verify_python_descriptor_frontier.py',['--engine','chematic']),('morgan',R/'scripts/verify_python_morgan_frontier.py',['--engine','chematic']),('maccs',O/'verify_maccs.py',[])]:
  output=O/f'{arm}-{kind}.jsonl';cmd=[P,str(script),*args,'--output',str(output)]
  with (O/f'{arm}-{kind}.log').open('w') as log:subprocess.run(cmd,cwd=R,env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
  data=output.read_bytes();rows=[json.loads(x) for x in data.splitlines()];values=[x for x in rows if 'index' in x]
  assert len(values)==14072
  def errors(value):
   if isinstance(value,dict): return sum(1 if k in ['error','parse_error'] else errors(v) for k,v in value.items())
   if isinstance(value,list): return sum(map(errors,value))
   return 0
  failures=sum(map(errors,values));assert failures==0,(arm,kind,failures)
  def nonfinite(value):
   if isinstance(value,dict):
    return sum(int(not math.isfinite(struct.unpack('>d',bytes.fromhex(v))[0])) if k=='float_bits' else nonfinite(v) for k,v in value.items())
   if isinstance(value,list): return sum(map(nonfinite,value))
   return 0
  proof['outputs'].setdefault(arm,{})[kind]={'sha256':hashlib.sha256(data).hexdigest(),'rows':len(values),'errors':failures,'nonfinite_float_outcomes':sum(map(nonfinite,values))}
  print(arm,kind,'all 14072 outputs retained',flush=True)
 assert proof['outputs'][arm]['descriptor']['nonfinite_float_outcomes']==6
 for kind in ['descriptor','morgan','maccs']:
  if arm=='candidate':assert (O/f'baseline-{kind}.jsonl').read_bytes()==(O/f'candidate-{kind}.jsonl').read_bytes(),kind
 cmd=[P,'-m','pytest','-q','-rs','-p','no:cacheprovider','crates/chematic-py/tests']
 (O/f'{arm}-python-test-command.json').write_text(json.dumps({'command':cmd,'env':{'PYTHONPATH':env['PYTHONPATH'],'PYTHONDONTWRITEBYTECODE':'1'}},indent=2)+'\n')
 with (O/f'{arm}-python-tests.log').open('w') as log:subprocess.run(cmd,cwd=R,env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
 proof['python_tests'][arm]=(O/f'{arm}-python-tests.log').read_text().splitlines()[-1]
 print(arm,'full Python suite passed',flush=True)
assert proof['outputs']['baseline']==proof['outputs']['candidate']
proof['passed']=True;proof['counts']={'descriptor_outcomes':112576,'morgan_outputs':56288,'maccs_outcomes':28144}
(O/'parity.json').write_text(json.dumps(proof,indent=2)+'\n')
