from pathlib import Path
import subprocess,os,json,shutil,hashlib
O=Path(__file__).resolve().parent;A=Path('/private/tmp/chematic-coordinate-pairs-evidence-20261017');R=Path('/private/tmp/chematic-openbabel-writer-20261017')
previous=json.loads((A/'environment.json').read_text());assert previous['rust']==subprocess.check_output(['rustc','-Vv'],text=True)
old=json.loads((A/'binary-seal.json').read_text())['files_sha256']
for n in ['baseline-sdf','baseline-v3000','baseline-v3000_3d','openbabel-sdf','openbabel-v3000']:
 assert hashlib.sha256((O/n).read_bytes()).hexdigest()==old[n]
(O/'baseline-reuse.json').write_text(json.dumps({'baseline_commit':'34f141266a81264f42990f98acbe25c91d9074af','source_root':str(A/'baseline-source'),'original_build_command':json.loads((A/'baseline-build-command.json').read_text()),'binaries_unchanged_from_sealed_predecessor':{n:old[n] for n in ['baseline-sdf','baseline-v3000','baseline-v3000_3d','openbabel-sdf','openbabel-v3000']},'same_rust_compiler':True,'old_elapsed_samples_reused':False},indent=2)+'\n')
target=O/'candidate-target';assert not target.exists(),'Candidate target must start empty'
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none')
cmd=['cargo','build','-vv','--offline','--locked','--release','-p','chematic-mol']
for kind in ['sdf','v3000','v3000_3d']:cmd+=['--example','openbabel_'+kind+'_throughput']
(O/'candidate-build-command.json').write_text(json.dumps({'cwd':str(R),'command':cmd,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
with (O/'candidate-build.log').open('w') as log:subprocess.run(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
for kind in ['sdf','v3000','v3000_3d']:shutil.copy2(target/'release/examples'/('openbabel_'+kind+'_throughput'),O/('candidate-'+kind))
print('Unchanged sealed baseline reused; new candidate built in empty independent target',flush=True)
