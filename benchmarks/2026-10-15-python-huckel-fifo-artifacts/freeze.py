from pathlib import Path
import hashlib,json,difflib,subprocess,platform
s=Path('/private/tmp/chematic-huckel-fifo-source-20261015');r=Path('/private/tmp/chematic-huckel-fifo-20261015');b=json.loads((r/'baseline-source.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
changed=sorted(n for n,h in b['sha256'].items() if sha(s/n)!=h)
assert changed==['CHANGELOG.md','crates/chematic-perception/src/rdkit_model_hook.rs'],changed
new=['crates/chematic-perception/src/rdkit_model_hook_storage_tests.rs'];names=changed+new
oldnames={'CHANGELOG.md':'baseline-CHANGELOG.md','crates/chematic-perception/src/rdkit_model_hook.rs':'baseline-rdkit_model_hook.rs'}
patch=''
for n in names:
 old=(r/oldnames[n]).read_text().splitlines(keepends=True) if n in oldnames else []
 patch+='diff --git a/'+n+' b/'+n+'\n'
 if not old:patch+='new file mode 100644\n'
 patch+=''.join(difflib.unified_diff(old,(s/n).read_text().splitlines(keepends=True),fromfile='a/'+n if old else '/dev/null',tofile='b/'+n))
(r/'huckel-fifo.patch').write_text(patch)
check=r/'patch-reconstruction';check.mkdir()
for n,f in oldnames.items():
 p=check/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((r/f).read_bytes())
subprocess.run(['git','init','-q',str(check)],check=True,timeout=10)
subprocess.run(['git','apply',str(r/'huckel-fifo.patch')],cwd=check,check=True,timeout=10)
assert all(sha(check/n)==sha(s/n) for n in names)
files={n:sha(s/n) for n in sorted(list(b['sha256'])+new)}
m={'schema':'huckel-fifo-source-candidate/v1','source_version':'1.1.4','baseline_source':b,'candidate_sha256':files,'changed_production_files':changed[1:],'changed_documentation_files':['CHANGELOG.md'],'added_test_files':new,'added_public_api':[],'patch_sha256':sha(r/'huckel-fifo.patch'),'patch_reconstruction_exact':True,'no_live_git_or_main_verification':True,'versions':{k:subprocess.check_output(v,text=True).strip() for k,v in {'rustc':['rustc','--version'],'cargo':['cargo','--version'],'maturin':['maturin','--version']}.items()},'host':platform.platform(),'build':{'profile':'release','opt_level':3,'lto':True,'codegen_units':1,'features':['extension-module','native-inchi'],'jobs':1,'note':'Fresh ten-candidate baseline built before editing. Initial baseline build hit ENOSPC; retry log retained. Verified duplicate staging and uncompressed prior outcomes were reclaimed after exact compressed-packet comparison; every source, wheel and outcome retained. Only owned completed build targets are reclaimed between validation/build/timing phases.'},'wheels':{arm:sha(next((r/(arm+'-wheel')).glob('*.whl'))) for arm in ['baseline','candidate']},'output_parity':{kind:{'byte_identical':(r/('baseline-'+kind+'.jsonl')).read_bytes()==(r/('candidate-'+kind+'.jsonl')).read_bytes(),'sha256':sha(r/('candidate-'+kind+'.jsonl'))} for kind in ['descriptors','morgan']}}
assert all(x['byte_identical'] for x in m['output_parity'].values())
(r/'source-manifest.json').write_text(json.dumps(m,indent=2)+'\n');print('Frozen',len(files),'source/input files; exact outputs and patch reconstruction passed')
