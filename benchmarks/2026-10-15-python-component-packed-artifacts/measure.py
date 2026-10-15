from pathlib import Path
import subprocess,json,shutil,hashlib
r=Path('/private/tmp/chematic-component-packed-20261015')
assert 'ALL VALIDATION PASSED' in (r/'validation.log').read_text()
target=Path('/private/tmp/chematic-component-packed-target-20261015')
size=sum(p.stat().st_size for p in target.rglob('*') if p.is_file())
wheels={arm:hashlib.sha256(next((r/(arm+'-wheel')).glob('*.whl')).read_bytes()).hexdigest() for arm in ['baseline','candidate']}
shutil.rmtree(target)
(r/'completed-build-cache-cleanup.json').write_text(json.dumps({'removed_only_completed_owned_target':str(target),'logical_bytes':size,'retained_measured_wheel_sha256':wheels,'before_timing':True},indent=2)+'\n')
for stage in ['freeze','run-frontiers','run-followups']:
 with (r/(stage+'.log')).open('w') as f:subprocess.run(['python3',str(r/(stage+'.py'))],stdout=f,stderr=subprocess.STDOUT,check=True)
 print(stage,'passed',flush=True)
