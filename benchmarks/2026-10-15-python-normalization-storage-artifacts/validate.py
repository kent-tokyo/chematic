from pathlib import Path
import os,subprocess,json,hashlib,shutil
r=Path('/private/tmp/chematic-descriptor-profile-20261015');s=Path('/private/tmp/chematic-descriptor-profile-source-20261015');python='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-descriptor-profile-target-20261015';env['CARGO_BUILD_JOBS']='1'
def call(name,cmd,extra=None):
 e=env.copy();e.update(extra or {})
 with (r/(name+'.log')).open('w') as log: subprocess.run(cmd,cwd=s,env=e,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
b=json.loads((r/'baseline-source.json').read_text());assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in b['sha256'].items())
for kind in ['descriptors','morgan']:
 call('baseline-'+kind,[python,'scripts/verify_python_'+('descriptor' if kind=='descriptors' else kind)+'_frontier.py','--engine','chematic','--output',str(r/('baseline-'+kind+'.jsonl'))],{'PYTHONPATH':str(r/'baseline')})
shutil.rmtree(env['CARGO_TARGET_DIR']);(r/'baseline-cache-cleanup.json').write_text(json.dumps({'deleted':'owned completed baseline Cargo target only','retained':'baseline wheel, source, profile and all logs'}))
call('optimize',['python3',str(r/'optimize.py')])
call('format-fragments',['rustfmt','--edition','2024','--config','skip_children=true','crates/chematic-perception/src/aromaticity/model.rs','crates/chematic-perception/src/aromaticity/normalization_storage_tests.rs'])
call('format',['cargo','fmt','--all','--check'])
call('refactor',['python3','scripts/check_refactor_budget.py'])
call('native-tests',['cargo','test','--release','-p','chematic-chem','-p','chematic-perception','--lib'])
call('clippy',['cargo','clippy','--release','-p','chematic-chem','-p','chematic-perception','--all-targets','--','-D','warnings'])
shutil.rmtree(env['CARGO_TARGET_DIR']);(r/'native-cache-cleanup.json').write_text(json.dumps({'deleted':'owned completed native validation Cargo target only','retained':'all native validation logs and source'}))
call('candidate-build',['maturin','build','--release','--manifest-path','crates/chematic-py/Cargo.toml','--interpreter',python,'--out',str(r/'candidate-wheel')])
wheel=next((r/'candidate-wheel').glob('*.whl'))
call('candidate-install',[python,'-m','pip','install','--no-deps','--target',str(r/'candidate'),str(wheel)])
for kind in ['descriptors','morgan']:
 call('candidate-'+kind,[python,'scripts/verify_python_'+('descriptor' if kind=='descriptors' else kind)+'_frontier.py','--engine','chematic','--output',str(r/('candidate-'+kind+'.jsonl'))],{'PYTHONPATH':str(r/'candidate')})
 assert (r/('baseline-'+kind+'.jsonl')).read_bytes()==(r/('candidate-'+kind+'.jsonl')).read_bytes(),kind
call('python-tests',[python,'-m','pytest','crates/chematic-py/tests','-q','-p','no:cacheprovider'],{'PYTHONPATH':str(r/'candidate')})
print('ALL VALIDATION PASSED',flush=True)
