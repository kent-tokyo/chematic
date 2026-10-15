from pathlib import Path
import hashlib,json,difflib,subprocess,tempfile,os
from datetime import datetime,timezone
s=Path('/private/tmp/chematic-openbabel-integer-source-20261015');r=Path('/private/tmp/chematic-openbabel-integer-20261015');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();base=json.loads((r/'baseline-source.json').read_text());prior=json.loads(Path('/private/tmp/chematic-openbabel-count-zero-map-20261015/source-manifest.json').read_text())
original_baseline=json.loads(Path('/private/tmp/chematic-openbabel-ascii-20261015/source-manifest.json').read_text())
assert base['sha256']==original_baseline['baseline_source_files_sha256']
assert sha(r/'baseline-v3000')==original_baseline['binaries']['baseline']['sha256']
changed=sorted(n for n,h in base['sha256'].items() if sha(s/n)!=h);assert changed==['CHANGELOG.md','crates/chematic-mol/src/mol2000.rs','crates/chematic-mol/src/mol2000/writer.rs','crates/chematic-mol/src/mol3000/writer.rs'],changed
new=['crates/chematic-mol/src/mol2000/v3000_integer_tests.rs'];names=changed+new;oldfiles={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-mol/src/mol2000.rs':'baseline-mol2000.rs','crates/chematic-mol/src/mol2000/writer.rs':'baseline-v2000-writer.rs','crates/chematic-mol/src/mol3000/writer.rs':'baseline-v3000-writer.rs'};patch=''
for n in names:
 old=(r/oldfiles[n]).read_text().splitlines(True) if n in oldfiles else [];now=(s/n).read_text().splitlines(True);patch+=f'diff --git a/{n} b/{n}\n'
 if not old:patch+='new file mode 100644\n'
 patch+=''.join(difflib.unified_diff(old,now,fromfile='a/'+n if old else '/dev/null',tofile='b/'+n))
p=r/'candidate.patch';p.write_text(patch)
with tempfile.TemporaryDirectory(prefix='chematic-ascii-patch-') as name:
 temp=Path(name)
 for n,old in oldfiles.items():q=temp/n;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes((r/old).read_bytes())
 subprocess.run(['git','init','-q',str(temp)],check=True,timeout=10);subprocess.run(['git','apply',str(p)],cwd=temp,check=True,timeout=10);assert all((temp/n).read_bytes()==(s/n).read_bytes() for n in names)
files={n:sha(s/n) for n in sorted(list(base['sha256'])+new)};libs=prior['native_libraries_sha256'];old=Path('/private/tmp/chematic-openbabel-v3000-20261014');ob=old/'openbabel-v3000'
for p,h in libs.items():assert sha(Path(p))==h,p
assert sha(ob)==prior['binaries']['openbabel']['sha256']
for profile,info in prior['corpora'].items():assert sha(old/f'corpus-{profile}.molblocks')==info['output_sha256']
m={'schema':'openbabel-v3000-integer-source/v1','created_utc':datetime.now(timezone.utc).isoformat(),'source_version':'1.1.4','last_verified_git_base':base['last_verified_git_base'],'live_git_unavailable':True,'baseline_scope':base['scope'],'baseline_source_files_sha256':base['sha256'],'frozen_build_source_sha256':files,'source_files_sha256':{n:files[n] for n in ['Cargo.toml','Cargo.lock',*names,'crates/chematic-mol/examples/openbabel_v3000_throughput.rs','scripts/bench_openbabel_v3000_throughput.py','scripts/bench_openbabel_sdf_throughput.py','validation/cosmolkit_comparison/bench_corpus.py','tools/openbabel_v3000_throughput.cpp']},'only_changed_existing_source':changed,'added_source_files':[],'added_test_files':new,'candidate_patch_sha256':sha(r/'candidate.patch'),'patch_reconstruction_exact':True,'compiler':subprocess.check_output(['rustc','-Vv'],text=True),'cargo':subprocess.check_output(['cargo','--version'],text=True).strip(),'build':'Candidate built in a fresh dedicated target, release opt3/LTO/codegen-units1/jobs1. Exact baseline11 native binary was built freshly before editing in the preceding rejected ASCII trial and is reused after verifying all 10058 source hashes and binary identity. No samples are pooled.', 'baseline_artifact_origin':json.loads((r/'baseline-artifact-origin.json').read_text()),'rustflags':os.environ.get('RUSTFLAGS'),'native_libraries_sha256':libs,'corpora':prior['corpora'],'continuation_diagnostic':json.loads((r/'continuation-plan.json').read_text()),'corpus_provenance':'validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json; original source corpora/licenses and exact generated MOL blocks unchanged','binaries':{name:{'path':str(p),'sha256':sha(p)} for name,p in {'baseline':r/'baseline-v3000','candidate':r/'candidate-v3000','openbabel':ob}.items()},'prior_experiment_evidence':'The rejected ASCII lexer and all rejected direct-line, flat-geometry, two-spelling zero-z, exact/rounded-only capacity, FIFO and SmallVec trials remain archived and excluded. Latest packed component facts are included in both arms.'}
(r/'source-manifest.json').write_text(json.dumps(m,indent=2)+'\n');print('Frozen',len(files),'source/input files; exact patch reconstruction and pinned libraries/corpora verified')
