from pathlib import Path
import json,hashlib,shutil,ctypes,os,subprocess
active=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic');r=Path('/private/tmp/chematic-huckel-fifo-20261015');s=Path('/private/tmp/chematic-huckel-fifo-source-20261015');assert not (r/"baseline-source.json").exists();r.mkdir(exist_ok=True);s.mkdir(exist_ok=True);assert not any(s.iterdir())
records=[]
for name in ['chematic-descriptor-profile-20261015','chematic-openbabel-count-zero-map-20261015']:
 old=Path('/private/tmp')/name;d=old/'delivery'
 if not d.exists():continue
 index=json.loads((old/'delivery-sha256.json').read_text());assert (old/'delivery-verification.json').exists()
 for n,h in index.items():
  assert hashlib.sha256((d/n).read_bytes()).hexdigest()==h,n
  if hashlib.sha256((active/n).read_bytes()).hexdigest()!=h:
   original=Path('/private/tmp/chematic-descriptor-profile-source-20261015')/n
   assert name=='chematic-descriptor-profile-20261015' and hashlib.sha256(original.read_bytes()).hexdigest()==h,n
 size=sum(p.stat().st_size for p in d.rglob('*') if p.is_file());shutil.rmtree(d);records.append({'removed_verified_duplicate_staging':str(d),'bytes':size,'retained':'complete byte-identical attached packet; all independent source snapshots, measured wheels/binaries, outcomes and logs'})
(r/'duplicate-staging-cleanup.json').write_text(json.dumps(records,indent=2)+'\n')
m=json.loads(Path('/private/tmp/chematic-openbabel-count-zero-map-20261015/source-manifest.json').read_text());files=m['frozen_build_source_sha256'];lib=ctypes.CDLL('/usr/lib/libSystem.B.dylib',use_errno=True);clone=lib.clonefile;clone.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int];cloned=0
for n,h in files.items():
 origin=active/n;assert hashlib.sha256(origin.read_bytes()).hexdigest()==h,n;p=s/n;p.parent.mkdir(parents=True,exist_ok=True)
 if clone(os.fsencode(origin),os.fsencode(p),0)==0:cloned+=1
 else:shutil.copyfile(origin,p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==h,n
(r/'baseline-source.json').write_text(json.dumps({'sha256':files,'scope':'All ten previously accepted source candidates; all rejected experiments excluded from this baseline','last_verified_git_base':'cac480fb23cd8adeec6e5e9bed1f18677f87b392','live_git_unavailable':True,'apfs_cloned_files':cloned},indent=2)+'\n')
for n,name in [('CHANGELOG.md','baseline-CHANGELOG.md'),('crates/chematic-perception/src/rdkit_model_hook.rs','baseline-rdkit_model_hook.rs')]:shutil.copyfile(s/n,r/name)
env=os.environ.copy();env.update(CARGO_TARGET_DIR='/private/tmp/chematic-huckel-fifo-target-20261015',CARGO_BUILD_JOBS='1');python='/private/tmp/chematic-cycle-speed-20261014/env/bin/python'
with (r/'baseline-build.log').open('w') as log:subprocess.run(['maturin','build','--release','--manifest-path','crates/chematic-py/Cargo.toml','--interpreter',python,'--out',str(r/'baseline-wheel')],cwd=s,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
with (r/'baseline-install.log').open('w') as log:subprocess.run([python,'-m','pip','install','--no-deps','--target',str(r/'baseline'),str(next((r/'baseline-wheel').glob('*.whl')))],cwd=s,stdout=log,stderr=subprocess.STDOUT,check=True)
print('Fresh ten-candidate baseline built and installed; all',len(files),'files verified',flush=True)
