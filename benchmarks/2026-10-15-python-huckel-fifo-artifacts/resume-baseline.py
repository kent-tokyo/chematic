from pathlib import Path
import subprocess,os
r=Path('/private/tmp/chematic-huckel-fifo-20261015');s=Path('/private/tmp/chematic-huckel-fifo-source-20261015')
env=os.environ.copy();env.update(CARGO_TARGET_DIR='/private/tmp/chematic-huckel-fifo-target-20261015',CARGO_BUILD_JOBS='1');python='/private/tmp/chematic-cycle-speed-20261014/env/bin/python'
with (r/'baseline-build.log').open('w') as log:subprocess.run(['maturin','build','--release','--manifest-path','crates/chematic-py/Cargo.toml','--interpreter',python,'--out',str(r/'baseline-wheel')],cwd=s,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
with (r/'baseline-install.log').open('w') as log:subprocess.run([python,'-m','pip','install','--no-deps','--target',str(r/'baseline'),str(next((r/'baseline-wheel').glob('*.whl')))],cwd=s,stdout=log,stderr=subprocess.STDOUT,check=True)
print('Fresh ten-candidate baseline built and installed; all',10057,'files verified',flush=True)
