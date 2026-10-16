from pathlib import Path
import hashlib, json, os, platform, subprocess, sys
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
OLD=Path('/private/tmp/chematic-v3000-loss-work-evidence-20261016')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
base=json.loads((OUT/'baseline-source-manifest.json').read_text())
for p,v in base['source_files_sha256'].items():assert sha(OUT/'baseline-source'/p)==v,p
prior=json.loads((OLD/'source-manifest.json').read_text())
paths=set(base['source_files_sha256'])|set(prior['source_files_sha256'])|{'crates/chematic-mol/src/mol3000/record_stream_tests.rs'}
paths.update(json.loads((OUT/'existing-fixtures-restored.json').read_text()))
sources={p:sha(ROOT/p) for p in sorted(paths)}
changed=[p for p,v in base['source_files_sha256'].items() if sources[p]!=v]
assert changed==['crates/chematic-mol/src/mol3000.rs'],changed
native=prior['native_libraries_sha256']
assert all(sha(Path(p))==v for p,v in native.items())
manifest={'base_revision':base['base_revision'],'source_files_sha256':sources,'native_libraries_sha256':native,'changed':changed,'added':['crates/chematic-mol/src/mol3000/record_stream_tests.rs']}
(OUT/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
binaries=[f'{arm}-{kind}' for arm in ['baseline','candidate','openbabel'] for kind in ['sdf','v3000']]
corpora=[p.name for p in OUT.iterdir() if p.suffix in ('.sdf','.molblocks')]
(OUT/'binary-seal.json').write_text(json.dumps({'files_sha256':{p:sha(OUT/p) for p in sorted(binaries+corpora)}},indent=2)+'\n')
environment={'platform':platform.platform(),'python':sys.version,'load_average':os.getloadavg(),'base_revision':base['base_revision']}
for command in [['rustc','-Vv'],['cargo','-V'],['obabel','-V']]:environment[' '.join(command)]=subprocess.check_output(command,text=True)
from rdkit import Chem
environment['rdkit_observer']=Chem.rdBase.rdkitVersion
(OUT/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
sys.path.insert(0,str(ROOT/'scripts'))
import bench_openbabel_sdf_throughput as sdf
original=sdf.observe((OUT/'v3000-chains.molblocks').read_bytes())
continued=sdf.observe((OUT/'v3000-continued.molblocks').read_bytes())
equality=sdf.semantic_check(original,continued)
assert equality['passed'],equality
(OUT/'continuation-observer-equality.json').write_text(json.dumps(equality,indent=2)+'\n')
print('Sealed',len(sources),'source/input files;',len(native),'native library entries; continuation fields unchanged')
