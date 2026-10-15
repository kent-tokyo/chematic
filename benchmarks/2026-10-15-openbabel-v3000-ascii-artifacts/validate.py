from pathlib import Path
import os,subprocess,shutil,hashlib,json
r=Path('/private/tmp/chematic-openbabel-ascii-20261015');s=Path('/private/tmp/chematic-openbabel-ascii-source-20261015');env=os.environ.copy();env['CARGO_BUILD_JOBS']='1';target=Path('/private/tmp/chematic-openbabel-ascii-baseline-target-20261015');env['CARGO_TARGET_DIR']=str(target)
def call(name,cmd):
 with (r/(name+'.log')).open('w') as out:subprocess.run(cmd,cwd=s,env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
base=json.loads((r/'baseline-source.json').read_text());assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in base['sha256'].items())
call('baseline-build',['cargo','build','--release','--locked','-p','chematic-mol','--example','openbabel_v3000_throughput'])
exe=target/'release/examples/openbabel_v3000_throughput';shutil.copyfile(exe,r/'baseline-v3000');(r/'baseline-v3000').chmod(0o755);assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in base['sha256'].items());shutil.rmtree(target);print('Baseline binary retained; completed owned target removed before source edit',flush=True)
env['CARGO_TARGET_DIR']='/private/tmp/chematic-openbabel-ascii-target-20261015'
call('optimize',['python3',str(r/'optimize.py')])
call('format-new',['rustfmt','--edition','2024','crates/chematic-mol/src/mol3000/ascii_whitespace_tests.rs'])
call('format',['cargo','fmt','--all','--check'])
call('refactor',['python3','scripts/check_refactor_budget.py'])
call('native-library-tests',['cargo','test','--release','--locked','-p','chematic-mol','--lib'])
call('native-integration-tests',['cargo','test','--release','--locked','-p','chematic-mol','--test','stereo2d_ez_reader_integration','--test','stereo_reader_integration','--test','conformer_3d_io','--test','square_planar_mol_io','--test','v3000_semantic_interop','--test','v3000_lookup_contract'])
call('clippy',['cargo','clippy','--release','--locked','-p','chematic-mol','--all-targets','--','-D','warnings'])
call('runner-tests',['/private/tmp/chematic-cycle-speed-20261014/env/bin/python','-m','pytest','scripts/tests/test_bench_openbabel_v3000_throughput.py','scripts/tests/test_bench_openbabel_sdf_throughput.py','-q','-p','no:cacheprovider'])
call('candidate-build',['cargo','build','--release','--locked','-p','chematic-mol','--example','openbabel_v3000_throughput'])
exe=Path(env['CARGO_TARGET_DIR'])/'release/examples/openbabel_v3000_throughput';shutil.copyfile(exe,r/'candidate-v3000');(r/'candidate-v3000').chmod(0o755);print('ALL VALIDATION PASSED',flush=True)
