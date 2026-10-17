from pathlib import Path
import os,subprocess,json
O=Path(__file__).resolve().parent;P='/private/tmp/chematic-python-asa-env-20261017/bin/python'
env=os.environ.copy();env.update(PYTHONPATH=str(O/'baseline'),PYTHONDONTWRITEBYTECODE='1')
for offset in [0,5000]:
 cmd=[P,str(O/'profile-worker.py'),str(offset)]
 with (O/f'profile-worker-{offset}.json').open('w') as out,(O/f'profile-worker-{offset}.stderr').open('w') as err:
  p=subprocess.Popen(cmd,env=env,stdout=out,stderr=err)
  try:
   sample=subprocess.run(['/usr/bin/sample',str(p.pid),'2','-file',str(O/f'profile-{offset}.txt')],capture_output=True,text=True)
   assert sample.returncode==0,sample.stderr
   assert p.wait(timeout=90)==0
  finally:
   if p.poll() is None:p.terminate();p.wait()
 print(offset,'diagnostic profile completed',flush=True)
