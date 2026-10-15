from pathlib import Path
import json,hashlib,ctypes,os,shutil
r=Path('/private/tmp/chematic-tpsa-zero-20261015');s=Path('/private/tmp/chematic-tpsa-zero-source-20261015');previous=Path('/private/tmp/chematic-parity-lazy-20261015');active=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic');assert not r.exists() and not s.exists();r.mkdir();s.mkdir();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
b=json.loads((previous/'baseline-source.json').read_text());lib=ctypes.CDLL('/usr/lib/libSystem.B.dylib',use_errno=True);clone=lib.clonefile;clone.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int]
for n,h in b['sha256'].items():
 p=s/n;p.parent.mkdir(parents=True,exist_ok=True);assert sha(active/n)==h,n
 if clone(os.fsencode(active/n),os.fsencode(p),0)!=0:shutil.copyfile(active/n,p)
 assert sha(p)==h,n
(r/'baseline-source.json').write_text(json.dumps(b,indent=2)+'\n');shutil.copytree(previous/'baseline',r/'baseline');shutil.copytree(previous/'baseline-wheel',r/'baseline-wheel')
for n in ['baseline-build.log','baseline-install.log','comparator-environment.txt']:shutil.copyfile(previous/n,r/n)
for n,f in [('CHANGELOG.md','baseline-CHANGELOG.md'),('crates/chematic-chem/src/descriptors/tpsa.rs','baseline-tpsa.rs')]:shutil.copyfile(s/n,r/f)
wheel=next((r/'baseline-wheel').glob('*.whl'));assert sha(wheel)==sha(next((previous/'baseline-wheel').glob('*.whl')))
(r/'baseline-reuse.json').write_text(json.dumps({'reused_exact_fresh_baseline14_wheel_from':str(previous),'source_files_verified':len(b['sha256']),'wheel_sha256':sha(wheel),'first_experiment_rejected':True,'no_candidate_wheel_or_timing_reuse':True,'baseline_built_before_either_experiment_edit':True},indent=2)+'\n')
for n in ['run-frontiers.py','run-followups.py','run-python-safe.py','summarize.py','measure.py']:
 (r/n).write_text((previous/n).read_text().replace('parity-lazy','tpsa-zero'))
plan=json.loads((previous/'measurement-plan.json').read_text());plan.update(intended_improvement='RDKit-default cold TPSA: skip aromatic perception when the already-selected model view has no nitrogen or oxygen; preserve the existing zero fold result and all N/O paths',baseline_reuse='Exact freshly built baseline14 wheel from the separately archived rejected indexed-normalization experiment; all source hashes unchanged; new candidate wheel and all new timings',scope=plan['scope']+' Previous experiment and its negatives remain separate; this new chemical predicate does not resolve them.')
(r/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Baseline14 source restored in a separate private clone; exact existing baseline wheel reused with provenance')
