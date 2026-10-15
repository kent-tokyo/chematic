from pathlib import Path
import subprocess,time
r=Path('/private/tmp/chematic-descriptor-profile-20261015')
while True:
 status=(r/'validation.log').read_text()
 if 'Traceback' in status: raise RuntimeError('Validation failed; formal timing must not run')
 if 'ALL VALIDATION PASSED' in status: break
 time.sleep(2)
with (r/'measurement.log').open('w') as log:subprocess.run(['python3',str(r/'measure.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
print('All fixed initial matrices and predeclared adaptive follow-ups complete',flush=True)
