from pathlib import Path
import os,subprocess,json,shutil
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017')
target=O/'candidate-target';assert not target.exists()
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
cmd=['cargo','build','-vv','--offline','--locked','--release','-p','chematic-mol']
for kind in ['sdf','v3000','v3000_3d']:cmd+=['--example','openbabel_'+kind+'_throughput']
(O/'candidate-build-command.json').write_text(json.dumps({'cwd':str(R),'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
with (O/'candidate-build.log').open('w') as log:subprocess.run(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
for kind in ['sdf','v3000','v3000_3d']:shutil.copy2(target/'release/examples'/('openbabel_'+kind+'_throughput'),O/('candidate-'+kind))
print('Fresh outlined candidate built',flush=True)
