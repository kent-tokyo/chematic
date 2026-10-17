from pathlib import Path
import json,os,re,subprocess
A=Path('/private/tmp/chematic-coordinate-pairs-evidence-20261017'); O=Path(__file__).resolve().parent; R=Path('/private/tmp/chematic-openbabel-writer-20261017')
env=os.environ.copy();env.update(CARGO_TARGET_DIR=str(A/'native-target'),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_LTO='false',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16')
packages=['chematic-perception','chematic-mol','chematic-chem','chematic-fp']
test=['cargo','test','--offline','--locked','--release','--lib','--tests','--no-fail-fast'];clippy=['cargo','clippy','--offline','--locked','--release','--all-targets']
for package in packages: test+=['-p',package];clippy+=['-p',package]
clippy+=['--','-D','warnings']
commands={'native-tests':test,'clippy':clippy,'format':['cargo','fmt','--all','--','--check'],'refactor':['python3','scripts/check_refactor_budget.py'],'runner-tests':['/Library/Frameworks/Python.framework/Versions/3.13/bin/python3','-m','pytest','-q','-p','no:cacheprovider','scripts/tests/test_bench_openbabel_sdf_throughput.py','scripts/tests/test_bench_openbabel_v3000_throughput.py']}
(O/'native-check-commands.json').write_text(json.dumps({'commands':commands,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
for name,cmd in commands.items():
 with (O/(name+'.log')).open('w') as log: subprocess.run(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
rows=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',(O/'native-tests.log').read_text())
assert rows
(O/'validation.json').write_text(json.dumps({'passed':True,'native':{'passed':sum(int(x[0]) for x in rows),'failed':sum(int(x[1]) for x in rows),'ignored':sum(int(x[2]) for x in rows),'suites':len(rows)},'runner_log':(O/'runner-tests.log').read_text()},indent=2)+'\n')
