from pathlib import Path
import os,subprocess,json,hashlib,shutil,time
r=Path('/private/tmp/chematic-tpsa-zero-20261015');s=Path('/private/tmp/chematic-tpsa-zero-source-20261015');python='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-tpsa-zero-target-20261015';env['CARGO_BUILD_JOBS']='1'
def call(name,cmd,extra=None):
 e=env.copy();e.update(extra or {})
 with (r/(name+'.log')).open('w') as log:subprocess.run(cmd,cwd=s,env=e,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
call('format-test',['rustfmt','--edition','2024','crates/chematic-chem/src/descriptors/tpsa.rs','crates/chematic-chem/src/descriptors/tpsa_zero_tests.rs'])
call('format',['cargo','fmt','--all','--check'])
call('refactor',['python3','scripts/check_refactor_budget.py'])
call('native-tests',['cargo','test','--release','-p','chematic-chem','--lib'])
call('clippy',['cargo','clippy','--release','-p','chematic-chem','--all-targets','--','-D','warnings'])
(r/'native-cache-reuse.json').write_text(json.dumps({'retained_owned_validation_target_for_cargo_fingerprint_reuse':env['CARGO_TARGET_DIR'],'candidate_wheel_built_after_source_edit':True,'no_wheel_reuse':True}))
call('candidate-build',['maturin','build','--release','--manifest-path','crates/chematic-py/Cargo.toml','--interpreter',python,'--out',str(r/'candidate-wheel')])
call('candidate-install',[python,'-m','pip','install','--no-deps','--target',str(r/'candidate'),str(next((r/'candidate-wheel').glob('*.whl')))])
for kind in ['descriptors','morgan']:
 call('baseline-'+kind,[python,'scripts/verify_python_'+('descriptor' if kind=='descriptors' else kind)+'_frontier.py','--engine','chematic','--output',str(r/('baseline-'+kind+'.jsonl'))],{'PYTHONPATH':str(r/'baseline')})
 call('candidate-'+kind,[python,'scripts/verify_python_'+('descriptor' if kind=='descriptors' else kind)+'_frontier.py','--engine','chematic','--output',str(r/('candidate-'+kind+'.jsonl'))],{'PYTHONPATH':str(r/'candidate')})
 assert (r/('baseline-'+kind+'.jsonl')).read_bytes()==(r/('candidate-'+kind+'.jsonl')).read_bytes(),kind
call('python-tests',[python,'-m','pytest','crates/chematic-py/tests','-q','-p','no:cacheprovider'],{'PYTHONPATH':str(r/'candidate')});print('ALL VALIDATION PASSED',flush=True)
