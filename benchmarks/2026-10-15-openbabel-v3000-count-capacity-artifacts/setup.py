from pathlib import Path
import json,hashlib,shutil,os,subprocess
s=Path('/private/tmp/chematic-openbabel-count-capacity-source-20261015');r=Path('/private/tmp/chematic-openbabel-count-capacity-20261015');assert not s.exists() and not r.exists();s.mkdir();r.mkdir();prior=json.loads(Path('/private/tmp/chematic-descriptor-profile-20261015/source-manifest.json').read_text());files=prior['candidate_sha256'];active=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic')
for n,h in files.items():
 data=(active/n).read_bytes();assert hashlib.sha256(data).hexdigest()==h,n;p=s/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
b={'sha256':files,'scope':'All nine previously accepted source candidates including descriptor normalization storage; no rejected experiment','last_verified_git_base':'cac480fb23cd8adeec6e5e9bed1f18677f87b392','live_git_unavailable':True};(r/'baseline-source.json').write_text(json.dumps(b,indent=2)+'\n')
for n in ['CHANGELOG.md','crates/chematic-mol/src/mol3000.rs']:(r/('baseline-'+Path(n).name)).write_bytes((s/n).read_bytes())
env=os.environ.copy();env.update(CARGO_TARGET_DIR='/private/tmp/chematic-openbabel-count-capacity-target-20261015',CARGO_BUILD_JOBS='1')
with (r/'baseline-build.log').open('w') as log:subprocess.run(['cargo','build','--release','--locked','-p','chematic-mol','--example','openbabel_v3000_throughput'],cwd=s,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
p=Path(env['CARGO_TARGET_DIR'])/'release/examples/openbabel_v3000_throughput';shutil.copyfile(p,r/'baseline-v3000');(r/'baseline-v3000').chmod(0o755)
print('Fresh nine-candidate baseline built; all',len(files),'source/input files verified',flush=True)
