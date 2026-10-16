from pathlib import Path
import json, os, subprocess
OUT=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
env=os.environ.copy()
env.update(CARGO_TARGET_DIR='/private/tmp/chematic-model-hook-inline-test-target-20261016',CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_LTO='false',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16',CARGO_PROFILE_RELEASE_STRIP='none')
commands={
 'native-tests':['cargo','test','--offline','--locked','--release','-p','chematic-perception','-p','chematic-chem'],
 'clippy':['cargo','clippy','--offline','--locked','--release','-p','chematic-perception','--all-targets','--','-D','warnings'],
 'format':['cargo','fmt','--all','--','--check'],
 'refactor':['python3','scripts/check_refactor_budget.py'],
}
(OUT/'native-check-commands.json').write_text(json.dumps({'commands':commands,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
for name,cmd in commands.items():
 with (OUT/(name+'.log')).open('w') as log:
  subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
