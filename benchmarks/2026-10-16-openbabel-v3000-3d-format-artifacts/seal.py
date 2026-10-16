from pathlib import Path
import hashlib,json,os,platform,subprocess,sys
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016');OUT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
baseline=json.loads((OUT/'baseline-source-manifest.json').read_text())
assert all(sha(OUT/'baseline-source'/p)==v for p,v in baseline['source_files_sha256'].items())
paths=set(baseline['source_files_sha256'])
paths.update(p for p in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines() if p.startswith(('scripts/','validation/benchmark_corpora/','validation/cosmolkit_comparison/bench_corpus.py')) and (ROOT/p).is_file())
paths.update(json.loads((OUT/'existing-fixtures-restored.json').read_text()))
paths.add('crates/chematic-mol/src/mol3000/coordinate_format_tests.rs')
source={p:sha(ROOT/p) for p in sorted(paths)}
old=json.loads(Path('/private/tmp/chematic-v3000-record-stream-evidence-20261016/source-manifest.json').read_text())
native=old['native_libraries_sha256'];assert all(sha(Path(p))==v for p,v in native.items())
manifest={'base_revision':baseline['base_revision'],'source_files_sha256':source,'native_libraries_sha256':native,'changed':['crates/chematic-mol/src/mol3000/writer.rs'],'added':['crates/chematic-mol/examples/openbabel_v3000_3d_throughput.rs','crates/chematic-mol/src/mol3000/coordinate_format_tests.rs'],'baseline_build_snapshot':'Only Git-tracked inputs plus common 3D probe; no ignored Python bytecode','instrumentation_equivalence':'Frozen baseline probe and candidate probe equal after rustfmt normalization; raw hashes in frozen-writer-proof.json'}
(OUT/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
files=[f'{a}-{k}' for a in ['baseline','candidate'] for k in ['sdf','v3000','v3000-3d']]+['openbabel-sdf','openbabel-v3000']
files += [p.name for p in OUT.iterdir() if p.suffix in ('.sdf','.molblocks')]
(OUT/'binary-seal.json').write_text(json.dumps({'files_sha256':{n:sha(OUT/n) for n in sorted(files)}},indent=2)+'\n')
env={'platform':platform.platform(),'python':sys.version,'load_average':os.getloadavg(),'base_revision':baseline['base_revision']}
for cmd in [['rustc','-Vv'],['cargo','-V'],['obabel','-V']]:env[' '.join(cmd)]=subprocess.check_output(cmd,text=True)
from rdkit import Chem
env['rdkit_observer']=Chem.rdBase.rdkitVersion
(OUT/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
print('Sealed',len(source),'source/input files;',len(native),'native library entries;',len(files),'binaries/corpora')
