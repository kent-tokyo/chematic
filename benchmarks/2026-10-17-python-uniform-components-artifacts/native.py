from pathlib import Path
import json,os,re,subprocess
R=Path('/private/tmp/chematic-python-asa-preparation-20261017');O=Path(__file__).resolve().parent
E=os.environ.copy();E.update(CARGO_TARGET_DIR=str(O/'native-target'),CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_LTO='false',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16',PYO3_PYTHON='/private/tmp/chematic-python-asa-env-20261017/bin/python')
packages=['chematic-perception','chematic-smiles','chematic-chem','chematic-fp','chematic-mol'];test=['cargo','test','--offline','--locked','--release','--lib','--tests'];clippy=['cargo','clippy','--offline','--locked','--release','--all-targets']
for p in packages:test+=['-p',p];clippy+=['-p',p]
clippy+=['--','-D','warnings']
cmds={'native-tests':test,'clippy':clippy,'format':['cargo','fmt','--all','--','--check'],'refactor':['python3','scripts/check_refactor_budget.py'],'runner-tests':['/private/tmp/chematic-python-asa-env-20261017/bin/python','-m','pytest','-q','-p','no:cacheprovider','scripts/tests/test_bench_python_speed_frontier.py']}
(O/'native-check-commands.json').write_text(json.dumps({'commands':cmds,'env':{k:v for k,v in E.items() if k.startswith(('CARGO_','PYO3_'))}},indent=2)+'\n')
for name,cmd in cmds.items():
 with (O/(name+'.log')).open('w') as log: subprocess.run(cmd,cwd=R,env=E,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
rows=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',(O/'native-tests.log').read_text());assert rows
(O/'native-validation.json').write_text(json.dumps({'passed':True,'native':{'passed':sum(int(x[0]) for x in rows),'failed':sum(int(x[1]) for x in rows),'ignored':sum(int(x[2]) for x in rows),'suites':len(rows)}},indent=2)+'\n')
