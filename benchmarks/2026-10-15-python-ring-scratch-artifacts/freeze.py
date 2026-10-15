from pathlib import Path
import hashlib,json,difflib,subprocess,platform
s=Path('/private/tmp/chematic-ring-scratch-source-20261015');r=Path('/private/tmp/chematic-ring-scratch-20261015');b=json.loads((r/'baseline-source.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
changed=sorted(n for n,h in b['sha256'].items() if sha(s/n)!=h)
assert changed==['CHANGELOG.md','crates/chematic-perception/src/rdkit_parity/view.rs','crates/chematic-perception/src/sssr/candidates.rs'],changed
new=['crates/chematic-perception/src/rdkit_parity/component_scratch_tests.rs','crates/chematic-perception/src/sssr/ring_scratch_tests.rs'];names=changed+new
oldnames={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-perception/src/rdkit_parity/view.rs':'baseline-view.rs','crates/chematic-perception/src/sssr/candidates.rs':'baseline-candidates.rs'}
patch=''
for n in names:
 old=(r/oldnames[n]).read_text().splitlines(keepends=True) if n in oldnames else []
 patch+='diff --git a/'+n+' b/'+n+'\n'
 if not old:patch+='new file mode 100644\n'
 patch+=''.join(difflib.unified_diff(old,(s/n).read_text().splitlines(keepends=True),fromfile='a/'+n if old else '/dev/null',tofile='b/'+n))
(r/'ring-scratch.patch').write_text(patch)
check=r/'patch-reconstruction';check.mkdir()
for n,f in oldnames.items():
 p=check/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((r/f).read_bytes())
subprocess.run(['git','init','-q',str(check)],check=True,timeout=10)
subprocess.run(['git','apply',str(r/'ring-scratch.patch')],cwd=check,check=True,timeout=10)
assert all(sha(check/n)==sha(s/n) for n in names)
files={n:sha(s/n) for n in sorted(list(b['sha256'])+new)}
m={'schema':'ring-scratch-source-candidate/v1','source_version':'1.1.4','baseline_source':b,'candidate_sha256':files,'changed_production_files':changed[1:],'changed_documentation_files':['CHANGELOG.md'],'added_test_files':new,'added_public_api':[],'patch_sha256':sha(r/'ring-scratch.patch'),'patch_reconstruction_exact':True,'no_live_git_or_main_verification':True,'versions':{k:subprocess.check_output(v,text=True).strip() for k,v in {'rustc':['rustc','--version'],'cargo':['cargo','--version'],'maturin':['maturin','--version']}.items()},'host':platform.platform(),'build':{'profile':'release','opt_level':3,'lto':True,'codegen_units':1,'features':['extension-module','native-inchi'],'jobs':1,'note':'Exact original ten-candidate baseline wheel and all source/input hashes reused from the preceding rejected FIFO trial, built before both trials. Original baseline build logs retained. Candidate built freshly after inline-scratch change; timing observations are separate, never pooled.'},'wheels':{arm:sha(next((r/(arm+'-wheel')).glob('*.whl'))) for arm in ['baseline','candidate']},'output_parity':{kind:{'byte_identical':(r/('baseline-'+kind+'.jsonl')).read_bytes()==(r/('candidate-'+kind+'.jsonl')).read_bytes(),'sha256':sha(r/('candidate-'+kind+'.jsonl'))} for kind in ['descriptors','morgan']}}
assert all(x['byte_identical'] for x in m['output_parity'].values())
(r/'source-manifest.json').write_text(json.dumps(m,indent=2)+'\n');print('Frozen',len(files),'source/input files; exact outputs and patch reconstruction passed')
