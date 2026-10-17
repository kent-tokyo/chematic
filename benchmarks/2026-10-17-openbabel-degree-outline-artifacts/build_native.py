from pathlib import Path
import subprocess,sys
O=Path(__file__).resolve().parent
jobs=[]
for step in ['build','native']:
 log=(O/(step+'-driver.log')).open('w');p=subprocess.Popen([sys.executable,str(O/(step+'.py'))],stdout=log,stderr=subprocess.STDOUT);jobs.append((step,p,log));print(step,'starts',flush=True)
failed=[]
for step,p,log in jobs:
 code=p.wait();log.close();print(step,'exit',code,flush=True)
 if code:failed.append((step,code))
assert not failed,failed
