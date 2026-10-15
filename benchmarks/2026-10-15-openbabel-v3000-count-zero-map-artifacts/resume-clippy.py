from pathlib import Path
import subprocess,os,shutil,json,hashlib,time
s=Path('/private/tmp/chematic-openbabel-count-zero-map-source-20261015');r=Path('/private/tmp/chematic-openbabel-count-zero-map-20261015');env=os.environ.copy();env['CARGO_TARGET_DIR']='/private/tmp/chematic-openbabel-count-zero-map-target-20261015';env['CARGO_BUILD_JOBS']='1'
commands=[('clippy.log',['cargo','clippy','--release','--locked','-p','chematic-mol','--all-targets','--','-D','warnings']),('runner-tests.log',['/private/tmp/chematic-cycle-speed-20261014/env/bin/python','-m','pytest','scripts/tests/test_bench_openbabel_v3000_throughput.py','scripts/tests/test_bench_openbabel_sdf_throughput.py','-q','-p','no:cacheprovider']),('candidate-build.log',['cargo','build','--release','--locked','-p','chematic-mol','--example','openbabel_v3000_throughput'])]
for name,cmd in commands:
 with (r/name).open('w') as out:subprocess.run(cmd,cwd=s,env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
exe=Path(env['CARGO_TARGET_DIR'])/'release/examples/openbabel_v3000_throughput';shutil.copyfile(exe,r/'candidate-v3000');(r/'candidate-v3000').chmod(0o755);print('ALL VALIDATION PASSED',flush=True)
