from pathlib import Path
import subprocess,json,shutil,hashlib
R=Path('/private/tmp/chematic-openbabel-writer-20261017');O=Path(__file__).resolve().parent;A=Path('/private/tmp/chematic-coordinate-pairs-evidence-20261017')
base=json.loads((A/'source-proof.json').read_text())['base_commit']
p=R/'crates/chematic-perception/src/stereo_validation.rs';s=p.read_text();start=s.index('pub fn stereo_centers(');end=s.index('\n/// Summarise',start)
frozen=s[start:end].replace('pub fn stereo_centers(', 'fn frozen_stereo_centers(',1)
old=subprocess.check_output(['git','show',base+':crates/chematic-perception/src/stereo_validation.rs'],cwd=R,text=True);bs=old.index('pub fn stereo_centers(');be=old.index('\n/// Summarise',bs)
assert s[start:end]==old[bs:be]
needle='''        // Skip aromatics and obvious non-centers.
        if atom.aromatic {
            continue;
        }
'''
replacement='''        // Every emitted shape needs at least three heavy neighbors,
        // including a constrained ring nitrogen with a lone pair.
        if atom.aromatic || mol.degree(idx) < 3 {
            continue;
        }
'''
assert s.count(needle)==1;s=s.replace(needle,replacement)
s+='\n#[cfg(test)]\n#[path = "stereo_degree_tests.rs"]\nmod degree_tests;\n';p.write_text(s)
t=Path('/private/tmp/chematic-stereo-shape-filter-20261017/crates/chematic-perception/src/stereo_shape_tests.rs').read_text()
ts=t.index('fn frozen_stereo_centers(');te=t.index('\nfn compare(',ts)
t=t[:ts]+frozen+t[te:]
(R/'crates/chematic-perception/src/stereo_degree_tests.rs').write_text(t)
(O/'source-proof.json').write_text(json.dumps({'base_commit':base,'predecessor_candidate_commit':'774ddc180d7389b95ed30ab7fbae6bcab65f2d64','predecessor_evidence_commit':'4c374951','predecessor_gate_passed':False,'frozen_classifier_sha256':hashlib.sha256(old[bs:be].encode()).hexdigest(),'transformation':'Retain decimal-pair formatting and reject degree<3 before local Vec/H-count/signature work. Every returned shape has at least three heavy neighbors; find_sssr still runs first. No changed valid-center branch or public API.','baseline_source_reused':str(A/'baseline-source'),'prior_timing_samples_reused':False},indent=2)+'\n')
(O/'baseline-source').symlink_to(A/'baseline-source',target_is_directory=True)
plan=json.loads((A/'measurement-plan.json').read_text());plan['candidate']='Decimal-pair coordinate writes plus a necessary degree>=3 classifier guard before branch allocation/signatures'
plan['primary']=['sdf-0/write','sdf-5000/write','v3000_3d-branched-spiral/parse','v3000_3d-branched-offset/parse']
plan['adoption']='Exact source snapshots/rewrites; frozen complete classifier matches 14072 graphs and 6144 five-atom/H overrides, numeric formatter oracles retained; all four primary paired CI lower bounds>1; every unrounded source median>=0.99; disclose negative intervals'
plan['build']='Reuse the exact sealed independent baseline binaries/source/build logs from formatter-only experiment; no prior elapsed samples. Candidate gets a new empty target, Rust opt3/full LTO/1 CGU/1 job, incremental false, strip none. Verify actual production code for both changes.'
plan['predecessor']='Formatter-only candidate rejected at branched-spiral parse median0.983815; add necessary degree guard to eliminate already-profiled impossible-center branch work. New source and independent timing epoch, no pooling or threshold relaxation.'
for p in plan['profiles']:shutil.copy2(A/p['path'],O/p['path'])
for n in ['baseline-build.log','baseline-build-command.json','baseline-prerequisite-command.json','baseline-prerequisite-tests.log','prerequisite-targeted-tests.log','prerequisite-commit.log','source-copy-correction.json']:
 shutil.copy2(A/n,O/n)
for n in ['baseline-sdf','baseline-v3000','baseline-v3000_3d','openbabel-sdf','openbabel-v3000']:shutil.copy2(A/n,O/n)
(O/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
for n in ['formal.py','native.py','profile.py']:
 text=(A/n).read_text()
 if n=='formal.py':
  text=text.replace("if key in plan['primary'] and mode=='write': primary[key]=ratio","if key+'/'+mode in plan['primary']: primary[key+'/'+mode]=ratio").replace('len(primary)==2','len(primary)==4')
 if n=='native.py':
  text=text.replace("packages=['chematic-mol']","packages=['chematic-perception','chematic-mol','chematic-chem','chematic-fp']")
  text=text.replace("CARGO_TARGET_DIR=str(O/'native-target')","CARGO_TARGET_DIR=str(A/'native-target')")
  text=text.replace("O=Path(__file__).resolve().parent; R=", "A=Path('/private/tmp/chematic-coordinate-pairs-evidence-20261017'); O=Path(__file__).resolve().parent; R=")
 if n=='profile.py':
  text=text.replace("str(O/'baseline-sdf'),str(O/'sdf-5000.sdf'),'write','1000','1000'", "str(O/'baseline-v3000_3d'),str(O/'v3000_3d-branched-spiral.molblocks'),'parse','100','64'")
  text=text.replace('prepared read ancestry from writer attribution and all diagnostic elapsed values from formal matrices','prepared read ancestry from parser attribution and all diagnostic elapsed values from formal matrices')
 (O/n).write_text(text)
print('New candidate plan fixed with four primary lanes; old timing samples excluded')
