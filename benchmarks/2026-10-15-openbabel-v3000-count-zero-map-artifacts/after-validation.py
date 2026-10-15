from pathlib import Path
import subprocess,time
r=Path('/private/tmp/chematic-openbabel-count-zero-map-20261015')
while True:
 status=(r/'validation-pipeline.log').read_text()
 if 'Traceback' in status:raise RuntimeError('Validation failed; no timing may start')
 if 'ALL VALIDATION PASSED' in status:break
 time.sleep(2)
with (r/'measurement.log').open('w') as f:subprocess.run(['python3',str(r/'measure.py')],stdout=f,stderr=subprocess.STDOUT,check=True)
print('All fixed profiles, continuation diagnostic and selected single follow-ups complete',flush=True)
