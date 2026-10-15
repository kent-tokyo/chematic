import subprocess,os,time,json
from pathlib import Path
r=Path('/private/tmp/chematic-parity-lazy-20261015');py='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';env=os.environ.copy();env['PYTHONPATH']=str(r/'baseline')
for offset,op in [(5000,'rdkit_tpsa'),(5000,'labute_asa'),(0,'rdkit_tpsa'),(0,'labute_asa')]:
 label=f'sample-{offset}-{op}';log=r/(label+'-worker.jsonl')
 with log.open('w') as f:
  process=subprocess.Popen([py,str(r/'profile-worker.py'),str(offset),op],env=env,stdout=f,stderr=subprocess.STDOUT)
  for _ in range(100):
   if log.stat().st_size:break
   if process.poll() is not None:raise RuntimeError('Worker exited before sampling: '+log.read_text())
   time.sleep(.05)
  ready=json.loads(log.read_text().splitlines()[0]);assert ready['ready']
  with (r/(label+'-sample.log')).open('w') as out:subprocess.run(['/usr/bin/sample',str(process.pid),'5','1','-mayDie','-file',str(r/(label+'.txt'))],stdout=out,stderr=subprocess.STDOUT,check=True)
  assert process.wait(timeout=30)==0
 print(label,'complete',flush=True)
