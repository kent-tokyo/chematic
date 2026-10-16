from pathlib import Path
import hashlib, importlib.metadata, json, os, platform, subprocess, sys
OUT=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/chematic-python-ring-first-use-20261017')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base=json.loads((OUT/'baseline-source-manifest.json').read_text())
for path,expected in base['sha256'].items(): assert sha(OUT/'baseline-source'/path)==expected,path
paths=set(base['sha256'])|{'crates/chematic-perception/src/relevant_cycles_limits_tests.rs'}
paths.update(['scripts/bench_python_speed_frontier.py','scripts/verify_python_descriptor_frontier.py','scripts/verify_python_morgan_frontier.py','validation/cosmolkit_comparison/bench_corpus.py','validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi','validation/cosmolkit_comparison/corpora/stress-v1.smi'])
candidate={p:sha(ROOT/p) for p in sorted(paths)}
changed={p for p in base['sha256'] if base['sha256'][p]!=candidate[p]}
assert changed=={'crates/chematic-perception/src/relevant_cycles.rs','crates/chematic-perception/src/lib.rs','crates/chematic-chem/src/descriptors/rings_and_filters.rs'},changed
(OUT/'source-manifest.json').write_text(json.dumps({'base':base['base'],'candidate_sha256':candidate,'changed_build_files':sorted(changed)},indent=2)+'\n')
comparators={}
for name in ['rdkit','cosmolkit']:
 module=__import__(name)
 directory=Path(module.__file__).parent
 natives={str(p):sha(p) for p in sorted(directory.rglob('*')) if p.is_file() and p.suffix in ('.so','.dylib','.dll','.pyd')}
 assert natives
 comparators[name]={'module':str(directory),'runtime_version':module.__version__,'distribution_version':importlib.metadata.version(name),'native_files':natives}
(OUT/'comparator-manifest.json').write_text(json.dumps(comparators,indent=2)+'\n')
environment={'platform':platform.platform(),'machine':platform.machine(),'python':sys.version,'load_average':os.getloadavg(),'base':base['base']}
for cmd in ['rustc -Vv','cargo -V','maturin -V']:
 environment[cmd]=subprocess.check_output(cmd.split(),text=True)
(OUT/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
print('Sealed',len(candidate),'source/input files;',sum(len(x['native_files']) for x in comparators.values()),'package-local native-file entries',flush=True)
