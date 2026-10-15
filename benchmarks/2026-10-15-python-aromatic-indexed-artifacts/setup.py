from pathlib import Path
import json,hashlib,gzip,shutil,ctypes,os,subprocess
r=Path('/private/tmp/chematic-parity-lazy-20261015');s=Path('/private/tmp/chematic-parity-lazy-source-20261015');assert not r.exists() and not s.exists();r.mkdir();s.mkdir();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();active=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic');cleanup=[]
prior=json.loads(Path('/private/tmp/chematic-openbabel-ring-queue-20261015/source-manifest.json').read_text());files=prior['frozen_build_source_sha256'];lib=ctypes.CDLL('/usr/lib/libSystem.B.dylib',use_errno=True);clone=lib.clonefile;clone.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int]
for n,h in files.items():
 origin=active/n;p=s/n;p.parent.mkdir(parents=True,exist_ok=True);assert sha(origin)==h,n
 if clone(os.fsencode(origin),os.fsencode(p),0)!=0:shutil.copyfile(origin,p)
 assert sha(p)==h,n
base={'sha256':files,'scope':'All fourteen accepted source candidates, including inactive E/Z and inline ring queue; all rejected prototypes excluded','last_verified_git_base':'cac480fb23cd8adeec6e5e9bed1f18677f87b392'};(r/'baseline-source.json').write_text(json.dumps(base,indent=2)+'\n')
for n,f in [('CHANGELOG.md','baseline-CHANGELOG.md')]:shutil.copyfile(s/n,r/f)
env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-parity-lazy-baseline-target-20261015';env['CARGO_BUILD_JOBS']='1';py='/private/tmp/chematic-cycle-speed-20261014/env/bin/python'
with (r/'baseline-build.log').open('w') as f:subprocess.run(['maturin','build','--release','--manifest-path','crates/chematic-py/Cargo.toml','--interpreter',py,'--out',str(r/'baseline-wheel')],cwd=s,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
assert all(sha(s/n)==h for n,h in files.items());wheel=next((r/'baseline-wheel').glob('*.whl'));print('Fresh baseline14 wheel built before editing',flush=True)
with (r/'baseline-install.log').open('w') as f:subprocess.run([py,'-m','pip','install','--no-deps','--target',str(r/'baseline'),str(wheel)],stdout=f,stderr=subprocess.STDOUT,check=True)
size=sum(p.stat().st_size for p in Path(env['CARGO_TARGET_DIR']).rglob('*') if p.is_file());shutil.rmtree(env['CARGO_TARGET_DIR']);(r/'baseline-cache-cleanup.json').write_text(json.dumps({'removed_owned_completed_target':env['CARGO_TARGET_DIR'],'logical_bytes':size,'retained_wheel_sha256':sha(wheel),'baseline_files_verified':len(files)},indent=2)+'\n');print('Baseline installed; completed own target removed',flush=True)
