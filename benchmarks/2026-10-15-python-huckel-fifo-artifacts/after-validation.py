from pathlib import Path
import subprocess,time
r=Path('/private/tmp/chematic-huckel-fifo-20261015')
while True:
 status=(r/'validation.log').read_text()
 if 'Traceback' in status:raise RuntimeError('Validation failed; do not start timing')
 if 'ALL VALIDATION PASSED' in status:break
 time.sleep(2)
with (r/'measurement.log').open('w') as f:subprocess.run(['python3',str(r/'measure.py')],stdout=f,stderr=subprocess.STDOUT,check=True)
print('All initial matrices and predeclared adaptive follow-ups complete',flush=True)
