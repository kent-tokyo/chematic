from pathlib import Path
import os,subprocess,json,shutil
from concurrent.futures import ThreadPoolExecutor
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017')
def build(arm):
 root=O/'baseline-source' if arm=='baseline' else R;target=O/(arm+'-target');assert not target.exists()
 env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
 cmd=['cargo','build','-vv','--offline','--locked','--release','-p','chematic-mol']
 for kind in ['sdf','v3000','v3000_3d']:cmd+=['--example','openbabel_'+kind+'_throughput']
 (O/(arm+'-build-command.json')).write_text(json.dumps({'cwd':str(root),'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
 with (O/(arm+'-build.log')).open('w') as log:subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 for kind in ['sdf','v3000','v3000_3d']:shutil.copy2(target/'release/examples'/('openbabel_'+kind+'_throughput'),O/(arm+'-'+kind))
 print(arm,'built',flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:
 results=[pool.submit(build,arm) for arm in ['baseline','candidate']]
 for result in results:result.result()
for kind in ['sdf','v3000']:
 cmd=['clang++','-O3','-std=c++17','-I/opt/homebrew/include/openbabel3',str(R/'tools'/('openbabel_'+kind+'_throughput.cpp')),'-L/opt/homebrew/lib','-lopenbabel','-o',str(O/('openbabel-'+kind))]
 with (O/('openbabel-'+kind+'-build.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
print('All three engines share the steady-state probe instrumentation',flush=True)
