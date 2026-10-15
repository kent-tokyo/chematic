from pathlib import Path
import subprocess
r=Path('/private/tmp/chematic-tpsa-zero-20261015')
assert not list(r.glob('frontier-*.json'))
for stage in ['freeze','run-frontiers','run-followups']:
 with (r/(stage+'.log')).open('w') as f:subprocess.run(['python3',str(r/(stage+'.py'))],stdout=f,stderr=subprocess.STDOUT,check=True)
 print(stage,'passed',flush=True)
