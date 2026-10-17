from pathlib import Path
import subprocess,sys,json
O=Path(__file__).resolve().parent
for step in ['build','native','profile','seal','formal']:
 print(step,'starts',flush=True)
 with (O/(step+'-driver.log' if step in ['build','native'] else step+'.log')).open('w') as log:
  subprocess.run([sys.executable,str(O/(step+'.py'))],stdout=log,stderr=subprocess.STDOUT,check=True)
 print(step,'passed',flush=True)
