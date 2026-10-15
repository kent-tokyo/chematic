from pathlib import Path
import os,subprocess,json,hashlib,shutil,time
r=Path('/private/tmp/chematic-huckel-fifo-20261015');s=Path('/private/tmp/chematic-huckel-fifo-source-20261015');python='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-huckel-fifo-target-20261015';env['CARGO_BUILD_JOBS']='1'
def call(name,cmd,extra=None):
 e=env.copy();e.update(extra or {})
 with (r/(name+'.log')).open('w') as log:subprocess.run(cmd,cwd=s,env=e,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
call('candidate-morgan',[python,'scripts/verify_python_morgan_frontier.py','--engine','chematic','--output',str(r/'candidate-morgan.jsonl')],{'PYTHONPATH':str(r/'candidate')})
for kind in ['descriptors','morgan']:assert (r/('baseline-'+kind+'.jsonl')).read_bytes()==(r/('candidate-'+kind+'.jsonl')).read_bytes(),kind
call('python-tests',[python,'-m','pytest','crates/chematic-py/tests','-q','-p','no:cacheprovider'],{'PYTHONPATH':str(r/'candidate')})
print('ALL VALIDATION PASSED',flush=True)
