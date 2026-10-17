from pathlib import Path
import shutil,subprocess,json,sys,hashlib
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017');P=Path('/private/tmp/chematic-openbabel-writer-20261017')
base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip();assert base=='35b932f911851f51c5cdb89fa371f1c2c4a1b11a'
B=O/'baseline-source';B.mkdir()
for d in ['crates','tools']:shutil.copytree(R/d,B/d)
for f in ['Cargo.toml','Cargo.lock','README.md','LICENSE','LICENSE-MIT','LICENSE-APACHE']:
 if (R/f).exists():shutil.copy2(R/f,B/f)
p=R/'crates/chematic-perception/src/stereo_validation.rs';s=p.read_text();a=s.index('pub fn stereo_centers(');b=s.index('\n/// Summarise',a);frozen=s[a:b].replace('pub fn stereo_centers(','fn frozen_stereo_centers(',1)
needle='''        // Skip aromatics and obvious non-centers.
        if atom.aromatic {
            continue;
        }
'''
assert s.count(needle)==1;s=s.replace(needle,'''        // Every emitted shape needs at least three heavy neighbors,
        // including a constrained ring nitrogen with a lone pair.
        if atom.aromatic || mol.degree(idx) < 3 {
            continue;
        }
''')
s+='\n#[cfg(test)]\n#[path = "stereo_degree_tests.rs"]\nmod degree_tests;\n';p.write_text(s)
t=(P/'crates/chematic-perception/src/stereo_degree_tests.rs').read_text();a=t.index('fn frozen_stereo_centers(');b=t.index('\nfn compare(',a);t=t[:a]+frozen+t[b:]
(R/'crates/chematic-perception/src/stereo_degree_tests.rs').write_text(t)
sys.path.insert(0,str(O));from patch_harness import apply
apply(B);apply(R)
shutil.copy2(O/'bench_openbabel_steady_throughput.py',R/'scripts/bench_openbabel_steady_throughput.py')
shutil.copy2(O/'test_bench_openbabel_steady_throughput.py',R/'scripts/tests/test_bench_openbabel_steady_throughput.py')
common=[str(p.relative_to(B)) for p in B.rglob('*') if p.is_file() and ('steady_timing' in p.name or p.name in ['openbabel_sdf_throughput.rs','openbabel_v3000_throughput.rs','openbabel_v3000_3d_throughput.rs','openbabel_sdf_throughput.cpp','openbabel_v3000_throughput.cpp'])]
for name in common:assert (B/name).read_bytes()==(R/name).read_bytes(),name
(O/'source-proof.json').write_text(json.dumps({'base_commit':base,'production_change':'Necessary degree>=3 precondition only; numerical formatter unchanged','complete_frozen_classifier_sha256':hashlib.sha256(frozen.replace('fn frozen_stereo_centers(','pub fn stereo_centers(',1).encode()).hexdigest(),'common_instrumentation_files':common,'both_arms_have_byte_identical_common_instrumentation':True,'old_elapsed_samples_reused':False},indent=2)+'\n')
plan=json.loads((O/'prospective-plan.json').read_text());shutil.copy2(O/'prospective-plan.json',O/'prospective-plan-original.json');plan['baseline_commit']=base;plan['baseline_update']='Use newly merged #981/#986 main before any C samples; earlier sources and results remain unchanged.'
A=Path('/private/tmp/chematic-coordinate-pairs-evidence-20261017');old=json.loads((A/'measurement-plan.json').read_text());plan['profiles']=old['profiles']
for profile in plan['profiles']:shutil.copy2(A/profile['path'],O/profile['path'])
(O/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Latest-main baseline, degree-only production change, symmetric timing patch prepared')
