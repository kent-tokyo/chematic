from pathlib import Path
import json, os, shutil, subprocess, sys
OUT=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/chematic-3d-stereo-guard-20261017')
arm=sys.argv[1];assert arm in ('baseline','candidate')
source=OUT/'baseline-source' if arm=='baseline' else ROOT
target=Path('/private/tmp/chematic-3d-stereo-guard-candidate-verified-target-20261017')
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
cmd=['cargo','build','-vv','--offline','--locked','--release','-p','chematic-mol','--example','openbabel_v3000_throughput','--example','openbabel_sdf_throughput','--example','openbabel_v3000_3d_throughput']
(OUT/(arm+'-build-command.json')).write_text(json.dumps({'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
with (OUT/(arm+'-build.log')).open('w') as log:subprocess.run(cmd,cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
for kind in ['sdf','v3000','v3000-3d']:shutil.copy2(target/'release/examples'/('openbabel_'+kind.replace('-','_')+'_throughput'),OUT/(arm+'-'+kind))
print(arm,'built',flush=True)
