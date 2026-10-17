from pathlib import Path
import os,time,subprocess,sys,json
O=Path(__file__).resolve().parent
while not (O/'validation.json').exists():
 try:os.kill(18725,0)
 except ProcessLookupError:raise RuntimeError('Validation exited without a passing manifest; formal measurement stopped')
 time.sleep(2)
assert json.loads((O/'validation.json').read_text())['passed']
print('Validation passed; diagnostic profile starts',flush=True)
with (O/'profile.log').open('w') as log:subprocess.run([sys.executable,str(O/'profile.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
print('Diagnostic completed; formal fixed matrices start',flush=True)
with (O/'formal.log').open('w') as log:subprocess.run([sys.executable,str(O/'formal.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
print('Formal matrices completed',flush=True)
