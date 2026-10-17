from pathlib import Path
import hashlib,json,subprocess,platform,sys,re
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-degree-only-steady-20261017')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
proof=json.loads((O/'source-proof.json').read_text());base=proof['base_commit'];common=proof['common_instrumentation_sha256']
assert json.loads((O/'validation.json').read_text())['passed'];assert json.loads((O/'instrumentation-check.json').read_text())['passed']
candidate=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
tracked=subprocess.check_output(['git','ls-files'],cwd=R,text=True).splitlines()
source={p:sha(R/p) for p in tracked if (R/p).is_file() and (p.startswith(('crates/','tools/')) or p in ['Cargo.lock','Cargo.toml'] or p.startswith('scripts/bench_openbabel_') or p=='validation/cosmolkit_comparison/bench_corpus.py')}
baseline={str(p.relative_to(O/'baseline-source')):sha(p) for p in (O/'baseline-source').rglob('*') if p.is_file()}
for name,h in baseline.items():
 if name in common:assert h==common[name]==sha(R/name),name
 else:assert hashlib.sha256(subprocess.check_output(['git','show',base+':'+name],cwd=R)).hexdigest()==h,name
assert (R/'crates/chematic-mol/src/mol2000/writer.rs').read_bytes()==subprocess.check_output(['git','show',base+':crates/chematic-mol/src/mol2000/writer.rs'],cwd=R)
prefix=Path('/opt/homebrew/lib/libopenbabel.dylib').resolve().parents[1];native={str(p):sha(p) for p in prefix.rglob('*') if p.is_file()}
manifest={'base_commit':base,'candidate_commit':candidate,'source_files_sha256':source,'baseline_source_files_sha256':baseline,'native_libraries_sha256':native,'common_probe_files_sha256':common,'baseline_production_matches_git_except_common_probes':True}
(O/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
owners={}
for arm in ['baseline','candidate']:
 log=(O/(arm+'-build.log')).read_text()
 for crate in ['chematic_mol','chematic_perception']:assert 'rustc --crate-name '+crate+' ' in log
 binary=O/(arm+'-v3000_3d');symbols=subprocess.check_output(['nm',str(binary)],text=True)
 found=[line.split()[-1] for line in symbols.splitlines() if line.split() and line.split()[-1].endswith('perceive_tetrahedral_from_3d_with_mode')];assert len(found)==1,found
 raw=subprocess.check_output(['/Library/Developer/CommandLineTools/usr/bin/llvm-objdump','--disassemble-symbols='+found[0],str(binary)]);assert len(raw)>1000
 (O/(arm+'-stereo-owner-disassembly.txt')).write_bytes(raw)
 instructions='\n'.join(re.sub(r'^\s*[0-9a-f]+:\s*','',line) for line in raw.decode().splitlines() if re.match(r'^\s*[0-9a-f]+:',line))
 owners[arm]={'symbol':found[0],'instruction_sha256':hashlib.sha256(instructions.encode()).hexdigest(),'instructions':len(instructions.splitlines())}
assert owners['baseline']['instruction_sha256']!=owners['candidate']['instruction_sha256']
(O/'binary-source-proof.json').write_text(json.dumps({'source_contrast_valid':True,'both_targets_started_empty':True,'both_production_libraries_recompiled_in_both_arms':True,'base_commit':base,'candidate_commit':candidate,'baseline_git_matches_except_seven_symmetric_probe_files':True,'common_probe_patch_matches_both_arms':True,'numeric_formatter_unchanged_from_base':True,'inlined_classifier_owner':owners},indent=2)+'\n')
import rdkit,numpy
(O/'environment.json').write_text(json.dumps({'host':platform.platform(),'python':sys.version,'rdkit':rdkit.__version__,'numpy':numpy.__version__,'rust':subprocess.check_output(['rustc','-Vv'],text=True),'clang':subprocess.check_output(['clang++','--version'],text=True),'openbabel':subprocess.check_output(['obabel','-V'],text=True),'native_prefix':str(prefix)},indent=2)+'\n')
sealed={p.name:sha(p) for p in O.iterdir() if p.is_file() and (p.name.startswith(('baseline-','candidate-','openbabel-')) or p.suffix in ['.molblocks','.sdf'] or p.name in ['source-proof.json','measurement-plan.json','source-manifest.json','binary-source-proof.json','formal.py','seal.py','setup.py','patch_harness.py','build.py','run.py','native.py','profile.py','instrument.py','validation.json','native-check-commands.json','instrumentation-check.json','environment.json','diagnostic-profile-scope.json'])}
(O/'binary-seal.json').write_text(json.dumps({'files_sha256':sealed},indent=2)+'\n')
print('Both actual production compilations and symmetric steady probes verified',flush=True)
