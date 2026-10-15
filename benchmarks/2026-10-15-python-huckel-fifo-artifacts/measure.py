from pathlib import Path
import subprocess,json,shutil,hashlib
r=Path('/private/tmp/chematic-huckel-fifo-20261015')
assert 'ALL VALIDATION PASSED' in (r/'validation.log').read_text()
assert (r/'completed-build-cache-cleanup.json').exists()
for stage in ['freeze','run-frontiers','run-followups']:
 with (r/(stage+'.log')).open('w') as f:subprocess.run(['python3',str(r/(stage+'.py'))],stdout=f,stderr=subprocess.STDOUT,check=True)
 print(stage,'passed',flush=True)
