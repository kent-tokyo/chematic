from pathlib import Path
import sys,os,subprocess,json,shutil,zipfile
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-python-asa-preparation-20261017');P='/private/tmp/chematic-python-asa-env-20261017/bin/python'
arm=sys.argv[1];assert arm in ['baseline','candidate'];root=O/'baseline-source' if arm=='baseline' else R
if arm=='baseline':
 root.mkdir()
 for name in ['crates','tools']: shutil.copytree(R/name,root/name)
 for name in ['Cargo.toml','Cargo.lock','README.md','LICENSE','LICENSE-MIT','LICENSE-APACHE']:
  if (R/name).exists(): shutil.copy2(R/name,root/name)
 (O/'baseline-revision.txt').write_text(subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True))
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(O/(arm+'-target')),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none',PYO3_PYTHON=P)
cmd=['/Library/Frameworks/Python.framework/Versions/3.13/bin/maturin','build','--offline','--locked','--release','--manifest-path',str(root/'crates/chematic-py/Cargo.toml'),'--interpreter',P,'--out',str(O/(arm+'-wheel')),'-vv']
(O/(arm+'-build-command.json')).write_text(json.dumps({'command':cmd,'cwd':str(root/'crates/chematic-py'),'env':{k:v for k,v in env.items() if k.startswith(('CARGO_','PYO3_'))}},indent=2)+'\n')
with (O/(arm+'-build.log')).open('w') as log: subprocess.run(cmd,cwd=root/'crates/chematic-py',env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
wheels=list((O/(arm+'-wheel')).glob('*.whl'));assert len(wheels)==1
with zipfile.ZipFile(wheels[0]) as archive: archive.extractall(O/arm)
print(arm,'wheel built',flush=True)
