from pathlib import Path
import os,subprocess,json,hashlib
r=Path('/private/tmp/chematic-aromatic-components-20261015');s=Path('/private/tmp/chematic-aromatic-components-source-20261015');python='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-aromatic-components-target-20261015';env['CARGO_BUILD_JOBS']='1'
def call(name,cmd,extra=None):
 e=env.copy();e.update(extra or {})
 with (r/(name+'.log')).open('w') as log: subprocess.run(cmd,cwd=s,env=e,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
b=json.loads((r/'baseline-source.json').read_text());assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in b['sha256'].items())
wheel=next((r/'baseline-wheel').glob('*.whl'))
call('baseline-install',[python,'-m','pip','install','--no-deps','--target',str(r/'baseline'),str(wheel)])
for kind in ['descriptors','morgan']:
 call('baseline-'+kind,[python,'scripts/verify_python_'+('descriptor' if kind=='descriptors' else kind)+'_frontier.py','--engine','chematic','--output',str(r/('baseline-'+kind+'.jsonl'))],{'PYTHONPATH':str(r/'baseline')})
call('optimize',['python3',str(r/'optimize.py')])
call('format-new-test',['rustfmt','--edition','2024','crates/chematic-perception/src/rdkit_parity/component_typing_tests.rs'])
call('format',['cargo','fmt','--all','--check'])
call('refactor',['python3','scripts/check_refactor_budget.py'])
# Reuse release dependencies where possible to keep disk use bounded.
native={'CARGO_ENCODED_RUSTFLAGS':'-C\x1flink-arg=-undefined\x1f-C\x1flink-arg=dynamic_lookup'}
call('native-tests',['cargo','test','--release','-p','chematic-chem','-p','chematic-perception','--lib'],native)
call('clippy',['cargo','clippy','--release','-p','chematic-chem','-p','chematic-perception','--all-targets','--','-D','warnings'],native)
call('candidate-build',['maturin','build','--release','--manifest-path','crates/chematic-py/Cargo.toml','--interpreter',python,'--out',str(r/'candidate-wheel')])
wheel=next((r/'candidate-wheel').glob('*.whl'))
call('candidate-install',[python,'-m','pip','install','--no-deps','--target',str(r/'candidate'),str(wheel)])
for kind in ['descriptors','morgan']:
 call('candidate-'+kind,[python,'scripts/verify_python_'+('descriptor' if kind=='descriptors' else kind)+'_frontier.py','--engine','chematic','--output',str(r/('candidate-'+kind+'.jsonl'))],{'PYTHONPATH':str(r/'candidate')})
 assert (r/('baseline-'+kind+'.jsonl')).read_bytes()==(r/('candidate-'+kind+'.jsonl')).read_bytes(),kind
call('python-tests',[python,'-m','pytest','crates/chematic-py/tests','-q','-p','no:cacheprovider'],{'PYTHONPATH':str(r/'candidate')})
print('ALL VALIDATION PASSED',flush=True)
