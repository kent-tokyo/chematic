from pathlib import Path
import subprocess,os,shutil,json,hashlib
s=Path('/private/tmp/chematic-openbabel-lines-source-20261015');r=Path('/private/tmp/chematic-openbabel-lines-20261015');env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-openbabel-lines-target-20261015';env['CARGO_BUILD_JOBS']='1';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();
bin=Path(env['CARGO_TARGET_DIR'])/'release/examples/openbabel_v3000_throughput'
commands=[('native-tests-retry.log',['cargo','test','--release','--locked','-p','chematic-mol']),('clippy.log',['cargo','clippy','--release','--locked','-p','chematic-mol','--all-targets','--','-D','warnings']),('runner-tests.log',['/private/tmp/chematic-cycle-speed-20261014/env/bin/python','-m','pytest','scripts/tests/test_bench_openbabel_v3000_throughput.py','scripts/tests/test_bench_openbabel_sdf_throughput.py','-q','-p','no:cacheprovider']),('candidate-build.log',['cargo','build','--release','--locked','-p','chematic-mol','--example','openbabel_v3000_throughput'])]
for name,cmd in commands:
 with (r/name).open('w') as f:subprocess.run(cmd,cwd=s,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
shutil.copyfile(bin,r/'candidate-v3000');(r/'candidate-v3000').chmod(0o755);print('ALL VALIDATION PASSED',flush=True)
