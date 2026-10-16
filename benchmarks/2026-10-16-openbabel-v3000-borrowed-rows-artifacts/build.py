from pathlib import Path
import json, os, shutil, subprocess, sys
OUT=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
arm=sys.argv[1];assert arm in ('baseline','candidate')
source=OUT/'baseline-source' if arm=='baseline' else ROOT
target=Path('/private/tmp/chematic-v3000-record-stream-target-20261016')
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
cmd=['cargo','build','--offline','--locked','--release','-p','chematic-mol','--example','openbabel_v3000_throughput','--example','openbabel_sdf_throughput']
(OUT/(arm+'-build-command.json')).write_text(json.dumps({'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
with (OUT/(arm+'-build.log')).open('w') as log:subprocess.run(cmd,cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
for kind in ['sdf','v3000']:shutil.copy2(target/'release/examples'/('openbabel_'+kind+'_throughput'),OUT/(arm+'-'+kind))
print(arm,'built',flush=True)
