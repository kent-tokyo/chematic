from pathlib import Path
import json,hashlib,ctypes,os,shutil,subprocess,time
r=Path('/private/tmp/chematic-openbabel-ring-queue-20261015');s=Path('/private/tmp/chematic-openbabel-ring-queue-source-20261015');assert not r.exists() and not s.exists();r.mkdir();s.mkdir();active=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic');prior=Path('/private/tmp/chematic-openbabel-ez-gate-20261015');m=json.loads((prior/'source-manifest.json').read_text());files=m['frozen_build_source_sha256'];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();lib=ctypes.CDLL('/usr/lib/libSystem.B.dylib',use_errno=True);clone=lib.clonefile;clone.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int]
for n,h in files.items():
 p=s/n;p.parent.mkdir(parents=True,exist_ok=True);assert sha(active/n)==h,n
 if clone(os.fsencode(active/n),os.fsencode(p),0)!=0:shutil.copyfile(active/n,p)
 assert sha(p)==h,n
base={'sha256':files,'scope':'All thirteen accepted source candidates, including inactive E/Z guard; all rejected experiments excluded','last_verified_git_base':m['last_verified_git_base']};(r/'baseline-source.json').write_text(json.dumps(base,indent=2)+'\n')
for n,f in [('CHANGELOG.md','baseline-CHANGELOG.md'),('crates/chematic-perception/src/stereo2d_ez_direction.rs','baseline-ez.rs')]:shutil.copyfile(s/n,r/f)
assert sha(prior/'candidate-v3000')==m['binaries']['candidate']['sha256'];shutil.copyfile(prior/'candidate-v3000',r/'baseline-v3000');(r/'baseline-v3000').chmod(0o755)
for name in ['candidate-build.log','native-library-tests.log','native-integration-tests.log','clippy.log']:shutil.copyfile(prior/name,r/('baseline-'+name))
(r/'baseline-artifact-origin.json').write_text(json.dumps({'source_run':str(prior),'source_manifest_sha256':sha(prior/'source-manifest.json'),'baseline_binary_from':'accepted source13 candidate-v3000','baseline_binary_sha256':sha(r/'baseline-v3000'),'all_current_baseline_files_verified':len(files),'no_timing_reuse':True,'built_before_current_edit':True},indent=2)+'\n')
old=Path('/private/tmp/chematic-openbabel-v3000-20261014')
for offset in [5000,0]:
 name='profile-'+str(offset)
 with (r/(name+'-worker.log')).open('w') as f:
  proc=subprocess.Popen([str(r/'baseline-v3000'),str(old/f'corpus-{offset}.molblocks'),'parse','800','1000'],stdout=f,stderr=subprocess.STDOUT)
  time.sleep(.2)
  with (r/(name+'-sample.log')).open('w') as log:subprocess.run(['/usr/bin/sample',str(proc.pid),'5','1','-mayDie','-file',str(r/(name+'.txt'))],stdout=log,stderr=subprocess.STDOUT,check=True)
  assert proc.wait(timeout=120)==0
 print(name,'complete',flush=True)
print('SETUP AND PROFILING COMPLETE',flush=True)
