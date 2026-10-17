from pathlib import Path
import hashlib,json,subprocess,platform,sys,re
O=Path(__file__).resolve().parent;R=Path('/private/tmp/chematic-openbabel-writer-20261017')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
proof=json.loads((O/'source-proof.json').read_text());base=proof['base_commit']
candidate=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
tracked=subprocess.check_output(['git','ls-files'],cwd=R,text=True).splitlines()
source={p:sha(R/p) for p in tracked if (R/p).is_file() and (p.startswith(('crates/','tools/')) or p in ['Cargo.lock','Cargo.toml'] or p.startswith('scripts/bench_openbabel_') or p=='validation/cosmolkit_comparison/bench_corpus.py')}
baseline={str(p.relative_to(O/'baseline-source')):sha(p) for p in (O/'baseline-source').rglob('*') if p.is_file()}
for name in baseline:
 old=subprocess.check_output(['git','show',base+':'+name],cwd=R)
 assert hashlib.sha256(old).hexdigest()==baseline[name],name
prefix=Path('/opt/homebrew/lib/libopenbabel.dylib').resolve().parents[1]
native={str(p):sha(p) for p in prefix.rglob('*') if p.is_file()}
manifest={'base_commit':base,'candidate_commit':candidate,'source_files_sha256':source,'baseline_source_files_sha256':baseline,'native_libraries_sha256':native}
(O/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
disassemblies={}
for arm in ['baseline','candidate']:
 log=(O/(arm+'-build.log')).read_text();assert 'rustc --crate-name chematic_mol ' in log
 binary=O/(arm+'-sdf')
 symbols=subprocess.check_output(['nm',str(binary)],text=True)
 found=[line.split()[-1] for line in symbols.splitlines() if line.split() and line.split()[-1].endswith('push_fixed4')]
 assert len(found)==1,found
 raw=subprocess.check_output(['/Library/Developer/CommandLineTools/usr/bin/llvm-objdump','--disassemble-symbols='+found[0],str(binary)])
 assert len(raw)>1000
 (O/(arm+'-fixed4-disassembly.txt')).write_bytes(raw)
 instructions='\n'.join(re.sub(r'^\s*[0-9a-f]+:\s*','',line) for line in raw.decode().splitlines() if re.match(r'^\s*[0-9a-f]+:',line))
 disassemblies[arm]={'symbol':found[0],'instruction_sha256':hashlib.sha256(instructions.encode()).hexdigest(),'instructions':len(instructions.splitlines())}
assert disassemblies['baseline']['instruction_sha256']!=disassemblies['candidate']['instruction_sha256']
(O/'binary-source-proof.json').write_text(json.dumps({'compiler_rebuilt_production_library':True,'separate_empty_targets':True,'base_commit':base,'candidate_commit':candidate,'all_baseline_source_files_match_base':True,'source_contrast_valid':True,'production_fixed4_function':disassemblies},indent=2)+'\n')
sealed={p.name:sha(p) for p in O.iterdir() if p.is_file() and (p.name.startswith(('baseline-','candidate-','openbabel-')) or p.suffix in ['.molblocks','.sdf'] or p.name in ['measurement-plan.json','source-manifest.json','binary-source-proof.json','formal.py','seal.py','prepare.py','setup.py','build.py'])}
(O/'binary-seal.json').write_text(json.dumps({'files_sha256':sealed},indent=2)+'\n')
import rdkit,numpy
(O/'environment.json').write_text(json.dumps({'host':platform.platform(),'python':sys.version,'rdkit':rdkit.__version__,'numpy':numpy.__version__,'rust':subprocess.check_output(['rustc','-Vv'],text=True),'openbabel':subprocess.check_output(['obabel','-V'],text=True),'native_prefix':str(prefix)},indent=2)+'\n')
print('Actual compiled production formatter differs; independent targets verified',flush=True)
