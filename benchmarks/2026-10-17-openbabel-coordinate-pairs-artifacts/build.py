from pathlib import Path
import os, subprocess, shutil, json
O=Path(__file__).resolve().parent; R=Path('/private/tmp/chematic-openbabel-writer-20261017')
for arm in ['baseline','candidate']:
 root=O/'baseline-source' if arm=='baseline' else R
 target=O/(arm+'-corrected-target')
 env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
 cmd=['cargo','build','-vv','--offline','--locked','--release','-p','chematic-mol']
 for kind in ['sdf','v3000','v3000_3d']: cmd+=['--example','openbabel_'+kind+'_throughput']
 (O/(arm+'-build-command.json')).write_text(json.dumps({'cwd':str(root),'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
 with (O/(arm+'-build.log')).open('w') as log: subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 for kind in ['sdf','v3000','v3000_3d']: shutil.copy2(target/'release/examples'/('openbabel_'+kind+'_throughput'),O/(arm+'-'+kind))
 print(arm,'built',flush=True)
env['CARGO_TARGET_DIR']=str(O/'baseline-corrected-target')
cmd=['cargo','test','--offline','--locked','--release','-p','chematic-mol','--test','mol2_3d_stereo']
(O/'baseline-prerequisite-command.json').write_text(json.dumps({'cwd':str(O/'baseline-source'),'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
with (O/'baseline-prerequisite-tests.log').open('w') as log:subprocess.run(cmd,cwd=O/'baseline-source',env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
print('baseline shared correctness prerequisite verified',flush=True)
