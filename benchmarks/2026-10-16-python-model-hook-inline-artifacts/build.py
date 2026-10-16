from pathlib import Path
import json, os, subprocess, sys

OUT=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
PYTHON='/private/tmp/chematic-ring-frontier-env-20261016/bin/python'
arm=sys.argv[1]
assert arm in ('baseline','candidate')
source=OUT/'baseline-source' if arm=='baseline' else ROOT
env=os.environ.copy()
env.update(CARGO_TARGET_DIR='/private/tmp/chematic-model-hook-inline-target-20261016',CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_OPT_LEVEL='3',CARGO_PROFILE_RELEASE_LTO='true',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='1',CARGO_PROFILE_RELEASE_STRIP='none',PYO3_PYTHON=PYTHON)
cmd=['maturin','build','--offline','--locked','--release','--manifest-path',str(source/'crates/chematic-py/Cargo.toml'),'--interpreter',PYTHON,'--out',str(OUT/(arm+'-wheel'))]
(OUT/(arm+'-build-command.json')).write_text(json.dumps({'cmd':cmd,'env':{k:v for k,v in env.items() if k.startswith(('CARGO_','PYO3_'))}},indent=2)+'\n')
with (OUT/(arm+'-build.log')).open('w') as log:
    subprocess.run(cmd,env=env,cwd=source,stdout=log,stderr=subprocess.STDOUT,check=True)
print(arm,'built',flush=True)
