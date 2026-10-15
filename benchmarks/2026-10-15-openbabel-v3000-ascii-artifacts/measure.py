from pathlib import Path
import subprocess,shutil,json,hashlib
r=Path('/private/tmp/chematic-openbabel-ascii-20261015');target=Path('/private/tmp/chematic-openbabel-ascii-target-20261015')
assert 'ALL VALIDATION PASSED' in (r/'validation-pipeline.log').read_text()
size=sum(p.stat().st_size for p in target.rglob('*') if p.is_file());retained={a:hashlib.sha256((r/(a+'-v3000')).read_bytes()).hexdigest() for a in ['baseline','candidate']}
shutil.rmtree(target);(r/'cache-cleanup.json').write_text(json.dumps({'removed_completed_owned_target':str(target),'logical_bytes':size,'retained_measured_binaries_sha256':retained,'after_validation_before_timing':True},indent=2)+'\n')
for name in ['freeze','run_formal','source_diagnostics']:
 with (r/(name+'.log')).open('w') as out:subprocess.run(['python3',str(r/(name+'.py'))],stdout=out,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
