from pathlib import Path
import os,subprocess,sys,time,json
E=Path(__file__).resolve().parent;R=Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic')
if '--worker' in sys.argv:
 import chematic
 rows=[s.split()[0] for s in (R/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi').read_text().splitlines() if s.strip()][:1000]
 deadline=time.monotonic()+8
 batches=0
 while time.monotonic()<deadline:
  mols=[chematic.from_smiles(s) for s in rows]
  values=[m.rdkit_ecfp_config(2,2048,include_chirality=True) for m in mols]
  assert all(len(v)==256 for v in values);batches+=1
 print(json.dumps({'batches':batches,'mode':'diagnostic only; parse plus first-use chiral Morgan','excluded_from_formal':True,'module':chematic.__file__}))
else:
 env=os.environ.copy();env['PYTHONPATH']=str(E/'baseline-site');env['PYTHONDONTWRITEBYTECODE']='1'
 with (E/'first-use-profile-worker.log').open('w') as output:
  process=subprocess.Popen([sys.executable,__file__,'--worker'],env=env,stdout=output,stderr=subprocess.STDOUT)
  subprocess.run(['/usr/bin/sample',str(process.pid),'5','1','-file',str(E/'first-use-baseline.sample.txt')],check=True,stdout=subprocess.DEVNULL)
  assert process.wait()==0
