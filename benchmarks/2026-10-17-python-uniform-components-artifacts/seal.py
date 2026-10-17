from pathlib import Path
import hashlib,json,subprocess,tarfile,io,platform,sys,importlib
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-python-asa-preparation-20261017')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base=(O/'baseline-revision.txt').read_text().strip();candidate=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
tracked=subprocess.check_output(['git','ls-files'],cwd=R,text=True).splitlines()
files={p:sha(R/p) for p in tracked if (R/p).is_file() and (p.startswith(('crates/','tools/','scripts/')) or p in ['Cargo.toml','Cargo.lock'] or p.startswith('validation/benchmark_corpora/') or p.startswith('validation/cosmolkit_comparison/'))}
baseline={str(p.relative_to(O/'baseline-source')):sha(p) for p in (O/'baseline-source').rglob('*') if p.is_file()}
paths=['crates','tools']+[p for p in ['Cargo.toml','Cargo.lock','README.md','LICENSE','LICENSE-MIT','LICENSE-APACHE'] if (O/'baseline-source'/p).is_file()]
archive=subprocess.check_output(['git','archive',base,*paths],cwd=R)
with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
 gitfiles={m.name:hashlib.sha256(tar.extractfile(m).read()).hexdigest() for m in tar.getmembers() if m.isfile()}
assert baseline==gitfiles, set(baseline)^set(gitfiles)
manifest={'base_commit':base,'candidate_commit':candidate,'candidate_sha256':files,'baseline_sha256':baseline}
(O/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
comparators={}
for name in ['rdkit','cosmolkit','numpy']:
 module=importlib.import_module(name);root=Path(module.__file__).parent
 native={str(p):sha(p) for p in root.rglob('*') if p.is_file() and p.suffix in ['.so','.dylib']}
 assert native,name
 comparators[name]={'runtime_version':module.__version__,'module':module.__file__,'native_files':native}
(O/'comparator-manifest.json').write_text(json.dumps(comparators,indent=2)+'\n')
wheels={};machine={}
for arm in ['baseline','candidate']:
 wheel=list((O/(arm+'-wheel')).glob('*.whl'));assert len(wheel)==1
 so=list((O/arm).rglob('*.so'));assert len(so)==1
 wheels[arm]={'path':str(wheel[0]),'sha256':sha(wheel[0]),'extensions':{str(p.relative_to(O/arm)):sha(p) for p in so}}
 log=(O/(arm+'-build.log')).read_text();assert 'rustc --crate-name chematic_perception ' in log
 symbols=subprocess.check_output(['nm',str(so[0])],text=True)
 targets=[line.split()[-1] for line in symbols.splitlines() if '19chematic_perception12rdkit_parity15parity_shortcut' in line and (line.split()[-1].startswith('__RNvNt') or line.split()[-1].startswith('__RNCNvNt'))]
 assert len(targets)==2,targets
 instructions=[]
 for idx,symbol in enumerate(targets):
  data=subprocess.check_output(['/Library/Developer/CommandLineTools/usr/bin/llvm-objdump','--disassemble-symbols='+symbol,str(so[0])])
  assert len(data)>1000
  (O/f'{arm}-shortcut-{idx}-disassembly.txt').write_bytes(data)
  import re
  body='\n'.join(re.sub(r'^\s*[0-9a-f]+:\s*','',line) for line in data.decode().splitlines() if re.match(r'^\s*[0-9a-f]+:',line))
  instructions.append({'symbol':symbol,'sha256':hashlib.sha256(body.encode()).hexdigest(),'instructions':len(body.splitlines())})
 machine[arm]=instructions
assert [m['sha256'] for m in machine['baseline']] != [m['sha256'] for m in machine['candidate']]
# Verify the actual frozen test body independently of its runtime outcomes.
source=subprocess.check_output(['git','show',base+':crates/chematic-perception/src/rdkit_parity/view.rs'],cwd=R,text=True)
old=source[source.index('fn parity_shortcut('):source.index('\n/// Whether the RDKit-parity aromatic view',source.index('fn parity_shortcut('))].strip()
test=(R/'crates/chematic-perception/src/rdkit_parity/uniform_component_tests.rs').read_text();frozen=test[test.index('fn frozen_parity_shortcut('):test.index('\nfn compare(')].replace('fn frozen_parity_shortcut(','fn parity_shortcut(',1).strip();assert old==frozen
(O/'binary-source-proof.json').write_text(json.dumps({'source_contrast_valid':True,'compiler_rebuilt_production_library':True,'independent_targets':True,'all_baseline_source_files_match_git':True,'complete_frozen_function_exact':True,'base_commit':base,'candidate_commit':candidate,'shortcut_machine_code':machine,'wheels':wheels},indent=2)+'\n')
inputs={p.name:sha(p) for p in O.iterdir() if p.is_file() and p.name in ['measurement-plan.json','source-manifest.json','comparator-manifest.json','binary-source-proof.json','instrumentation-check.json','bench-fixed-batches.py','formal.py','build.py','implement.py','prepare-protocol.py','check-instrumentation.py']}
(O/'protocol-seal.json').write_text(json.dumps({'files_sha256':inputs},indent=2)+'\n')
(O/'environment.json').write_text(json.dumps({'host':platform.platform(),'python':sys.version,'rust':subprocess.check_output(['rustc','-Vv'],text=True),'comparators':{k:v['runtime_version'] for k,v in comparators.items()}},indent=2)+'\n')
print('Source, wheels, actual shortcut machine code and common instrumentation sealed',flush=True)
