from pathlib import Path
import os,subprocess,time,json
O=Path(__file__).resolve().parent;P='/private/tmp/chematic-python-asa-env-20261017/bin/python';env=os.environ.copy();env.update(PYTHONPATH=str(O/'baseline'),PYTHONDONTWRITEBYTECODE='1')
with (O/'profile-worker-5000-retry1.json').open('w') as out,(O/'profile-worker-5000-retry1.stderr').open('w') as err:
 p=subprocess.Popen([P,str(O/'profile-worker.py'),'5000'],env=env,stdout=out,stderr=err)
 try:
  time.sleep(.3)
  r=subprocess.run(['/usr/bin/sample',str(p.pid),'2','-file',str(O/'profile-5000-retry1.txt')],capture_output=True,text=True)
  (O/'profile-5000-retry1-status.json').write_text(json.dumps({'returncode':r.returncode,'stderr':r.stderr,'scope':'Diagnostic only; initial attachment failure retained. No formal samples.'},indent=2)+'\n')
  assert r.returncode==0,r.stderr
  assert p.wait(timeout=90)==0
 finally:
  if p.poll() is None:p.terminate();p.wait()
