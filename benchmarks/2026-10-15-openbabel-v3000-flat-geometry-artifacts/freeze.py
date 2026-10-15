from pathlib import Path
import hashlib,json,difflib,subprocess,tempfile,os
from datetime import datetime,timezone
s=Path('/private/tmp/chematic-openbabel-flat-geometry-source-20261015');r=Path('/private/tmp/chematic-openbabel-flat-geometry-20261015');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();base=json.loads((r/'baseline-source.json').read_text());prior=json.loads(Path('/private/tmp/chematic-openbabel-token-20261015/source-manifest.json').read_text())
changed=sorted(n for n,h in base['sha256'].items() if sha(s/n)!=h);assert changed==['CHANGELOG.md','crates/chematic-mol/src/mol2000.rs','crates/chematic-mol/src/mol3000.rs'],changed
new=['crates/chematic-mol/src/mol3000/flat_geometry_tests.rs','validation/benchmark_corpora/flat-geometry-speed-v1.recipe.json','validation/benchmark_corpora/flat-geometry-speed-v1.manifest.json'];names=changed+new;oldfiles={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-mol/src/mol2000.rs':'baseline-mol2000.rs','crates/chematic-mol/src/mol3000.rs':'baseline-mol3000.rs'};patch=''
for n in names:
 old=(r/oldfiles[n]).read_text().splitlines(True) if n in oldfiles else [];now=(s/n).read_text().splitlines(True);patch+=f'diff --git a/{n} b/{n}\n'
 if not old:patch+='new file mode 100644\n'
 patch+=''.join(difflib.unified_diff(old,now,fromfile='a/'+n if old else '/dev/null',tofile='b/'+n))
p=r/'candidate.patch';p.write_text(patch)
with tempfile.TemporaryDirectory(prefix='chematic-lines-patch-') as name:
 temp=Path(name)
 for n,old in oldfiles.items():q=temp/n;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes((r/old).read_bytes())
 subprocess.run(['git','init','-q',str(temp)],check=True,timeout=10);subprocess.run(['git','apply',str(p)],cwd=temp,check=True,timeout=10)
 assert all((temp/n).read_bytes()==(s/n).read_bytes() for n in names)
files={n:sha(s/n) for n in sorted(list(base['sha256'])+new)};libs=prior['native_libraries_sha256'];old=Path('/private/tmp/chematic-openbabel-v3000-20261014');ob=old/'openbabel-v3000'
for p,h in libs.items():assert sha(Path(p))==h,p
assert sha(ob)==prior['binaries']['openbabel']['sha256']
for profile,info in prior['corpora'].items():assert sha(old/f'corpus-{profile}.molblocks')==info['output_sha256']
m={'schema':'openbabel-v3000-flat-geometry-source/v1','created_utc':datetime.now(timezone.utc).isoformat(),'source_version':'1.1.4','last_verified_git_base':base['last_verified_git_base'],'live_git_unavailable':True,'baseline_scope':base['scope'],'baseline_source_files_sha256':base['sha256'],'frozen_build_source_sha256':files,'source_files_sha256':{n:files[n] for n in ['Cargo.toml','Cargo.lock',*names,'crates/chematic-mol/examples/openbabel_v3000_throughput.rs','scripts/bench_openbabel_v3000_throughput.py','scripts/bench_openbabel_sdf_throughput.py','validation/cosmolkit_comparison/bench_corpus.py','tools/openbabel_v3000_throughput.cpp']},'only_changed_existing_source':changed,'added_source_files':new[1:],'added_test_files':new[:1],'candidate_patch_sha256':sha(r/'candidate.patch'),'patch_reconstruction_exact':True,'compiler':subprocess.check_output(['rustc','-Vv'],text=True),'cargo':subprocess.check_output(['cargo','--version'],text=True).strip(),'build':'cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput; fresh dedicated target; baseline source verified before editing; reuse unchanged dependencies between arms; CARGO_BUILD_JOBS=1','rustflags':os.environ.get('RUSTFLAGS'),'native_libraries_sha256':libs,'corpora':prior['corpora'],'continuation_diagnostic':json.loads((r/'continuation-plan.json').read_text()),'corpus_provenance':'validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json; original source corpora/licenses and exact generated MOL blocks unchanged','binaries':{name:{'path':str(p),'sha256':sha(p)} for name,p in {'baseline':r/'baseline-v3000','candidate':r/'candidate-v3000','openbabel':ob}.items()}}
m['rejected_direct_line_experiment']={'manifest_sha256':sha(Path('/private/tmp/chematic-openbabel-lines-20261015/source-manifest.json')),'decision':'rejected_for_adoption; no original-profile parse improvement and reproduced writer regressions; all evidence retained separately'}
fixture=json.loads((s/'validation/benchmark_corpora/flat-geometry-speed-v1.manifest.json').read_text())
assert sha(s/fixture['path'])==fixture['sha256']
assert fixture['oracle']['mol2000_sha256']==base['sha256']['crates/chematic-mol/src/mol2000.rs']
assert fixture['oracle']['mol3000_sha256']==base['sha256']['crates/chematic-mol/src/mol3000.rs']
(r/'source-manifest.json').write_text(json.dumps(m,indent=2)+'\n');print('Frozen',len(files),'files; exact patch reconstruction; libraries, harness and all original corpora unchanged')
