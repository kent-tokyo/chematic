from pathlib import Path
import subprocess,sys
O=Path(__file__).resolve().parent
jobs=[]
for step in ['build','native']:
 log=(O/(step+'-driver.log')).open('w');p=subprocess.Popen([sys.executable,str(O/(step+'.py'))],stdout=log,stderr=subprocess.STDOUT);jobs.append((step,p,log));print(step,'starts',flush=True)
failed=[]
for step,p,log in jobs:
 result=p.wait();log.close();print(step,'exit',result,flush=True)
 if result:failed.append((step,result))
assert not failed,failed
for step in ['profile','seal','formal']:
 print(step,'starts',flush=True)
 with (O/(step+'.log')).open('w') as log:subprocess.run([sys.executable,str(O/(step+'.py'))],stdout=log,stderr=subprocess.STDOUT,check=True)
 print(step,'passed',flush=True)
