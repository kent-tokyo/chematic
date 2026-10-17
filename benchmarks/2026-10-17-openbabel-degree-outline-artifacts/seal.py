from pathlib import Path
import hashlib,json,subprocess,re
O=Path(__file__).resolve().parent;C=Path('/private/tmp/chematic-degree-only-steady-evidence-20261017');R=Path('/private/tmp/chematic-degree-only-steady-20261017')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert json.loads((O/'validation.json').read_text())['passed']
assert json.loads((O/'instrumentation-check.json').read_text())['passed']
old=json.loads((O/'reused-baseline-source-manifest.json').read_text());oldseal=json.loads((O/'reused-baseline-binary-seal.json').read_text());plan=json.loads((O/'measurement-plan.json').read_text())
assert old['base_commit']==plan['base_commit']=='35b932f911851f51c5cdb89fa371f1c2c4a1b11a'
source={name:sha(R/name) for name in old['source_files_sha256']}
assert [name for name,value in source.items() if value!=old['source_files_sha256'][name]]==['crates/chematic-perception/src/stereo_validation.rs']
for name,value in old['baseline_source_files_sha256'].items():assert sha(C/'baseline-source'/name)==value
for name,value in old['native_libraries_sha256'].items():assert sha(Path(name))==value
for name,value in old['common_probe_files_sha256'].items():assert sha(R/name)==sha(C/'baseline-source'/name)==value
for name in ['baseline-sdf','baseline-v3000','baseline-v3000_3d','openbabel-sdf','openbabel-v3000','baseline-build.log']:
 assert sha(O/name)==oldseal['files_sha256'][name]
subprocess.run(['git','diff','--exit-code'],cwd=R,check=True,capture_output=True)
candidate=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
log=(O/'candidate-build.log').read_text()
for crate in ['chematic_mol','chematic_perception']:assert 'rustc --crate-name '+crate+' ' in log
binary=O/'candidate-v3000_3d';symbols=subprocess.check_output(['nm',str(binary)],text=True)
def owner(suffix):
 found=[line.split()[-1] for line in symbols.splitlines() if line.split() and line.split()[-1].endswith(suffix)];assert len(found)==1,(suffix,found)
 raw=subprocess.check_output(['/Library/Developer/CommandLineTools/usr/bin/llvm-objdump','--disassemble-symbols='+found[0],str(binary)])
 (O/(suffix+'-disassembly.txt')).write_bytes(raw)
 instructions='\n'.join(re.sub(r'^\s*[0-9a-f]+:\s*','',line) for line in raw.decode().splitlines() if re.match(r'^\s*[0-9a-f]+:',line))
 return {'symbol':found[0],'instruction_sha256':hashlib.sha256(instructions.encode()).hexdigest(),'instructions':len(instructions.splitlines())},raw
reader,raw=owner('perceive_tetrahedral_from_3d_with_mode');classifier,_=owner('stereo_centers')
assert classifier['symbol'].lstrip('_') in raw.decode(), '3D reader must call the outlined classifier'
assert reader['instructions']<json.loads((O/'reused-baseline-binary-source-proof.json').read_text())['inlined_classifier_owner']['baseline']['instructions']
proof={'base_commit':old['base_commit'],'candidate_commit':candidate,'source_contrast_valid':True,'baseline_binary_reused_from_C_without_elapsed_observations':True,'baseline_source_and_binaries_still_match_C_seals':True,'candidate_target_started_empty':True,'candidate_production_libraries_recompiled':True,'common_probe_patch_identical':True,'only_change_from_C_source':'inline(never) on stereo_centers','classifier_is_separate_called_function':True,'candidate_reader':reader,'candidate_classifier':classifier}
(O/'binary-source-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
manifest={'base_commit':old['base_commit'],'candidate_commit':candidate,'source_files_sha256':source,'baseline_source_root':str(C/'baseline-source'),'baseline_source_files_sha256':old['baseline_source_files_sha256'],'native_libraries_sha256':old['native_libraries_sha256'],'common_probe_files_sha256':old['common_probe_files_sha256']}
(O/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
import platform,sys,rdkit,numpy
oldenv=json.loads((O/'environment.json').read_text())
fresh={'host':platform.platform(),'python':sys.version,'rdkit':rdkit.__version__,'numpy':numpy.__version__,'rust':subprocess.check_output(['rustc','-Vv'],text=True),'clang':subprocess.check_output(['clang++','--version'],text=True),'openbabel':subprocess.check_output(['obabel','-V'],text=True),'native_prefix':oldenv['native_prefix']}
for key in ['rust','clang','openbabel','rdkit','numpy']:assert fresh[key]==oldenv[key],key
(O/'environment.json').write_text(json.dumps(fresh,indent=2)+'\n')
files={p.name:sha(p) for p in O.iterdir() if p.is_file() and (p.name.startswith(('baseline-','candidate-','openbabel-','reused-baseline-')) or p.suffix in ['.sdf','.molblocks'] or p.name in ['measurement-plan.json','source-manifest.json','binary-source-proof.json','formal.py','seal.py','build.py','native.py','instrument.py','validation.json','instrumentation-check.json','environment.json','semantic-reference.json'])}
(O/'binary-seal.json').write_text(json.dumps({'files_sha256':files},indent=2)+'\n')
print(json.dumps(proof),flush=True)
