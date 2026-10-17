from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys
E=Path(__file__).resolve().parent;R=Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic');arm=sys.argv[1];assert arm in ('baseline','candidate');target=E/(arm+'-target');assert not target.exists()
manifest=json.loads((E/(arm+'-source.json')).read_text())
for p,h in manifest['source_files_sha256'].items():assert hashlib.sha256((R/p).read_bytes()).hexdigest()==h,p
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
cmd=['cargo','build','-vv','--offline','--locked','--release','-p','chematic-mol']
for kind in ['sdf','v3000','v3000_3d']:cmd+=['--example','openbabel_'+kind+'_throughput']
(E/(arm+'-build-command.json')).write_text(json.dumps({'command':cmd,'cwd':str(R),'environment':{k:v for k,v in env.items() if k.startswith('CARGO_')},'target_started_empty':True,'rust':subprocess.check_output(['rustc','-Vv'],text=True)},indent=2)+'\n')
with (E/(arm+'-build.log')).open('w') as log:subprocess.run(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
for kind in ['sdf','v3000','v3000_3d']:shutil.copy2(target/'release/examples'/('openbabel_'+kind+'_throughput'),E/(arm+'-'+kind))
for p,h in manifest['source_files_sha256'].items():assert hashlib.sha256((R/p).read_bytes()).hexdigest()==h,p
print(arm,'built; source unchanged',flush=True)
