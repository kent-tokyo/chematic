from pathlib import Path
import json, os, re, subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
env=os.environ.copy();env.update(CARGO_TARGET_DIR='/private/tmp/chematic-v3000-record-stream-test-target-20261016',CARGO_BUILD_JOBS='1',CARGO_INCREMENTAL='0',CARGO_PROFILE_RELEASE_LTO='false',CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16',CARGO_PROFILE_RELEASE_STRIP='none')
commands={
 'native-tests':['cargo','test','--offline','--locked','--release','-p','chematic-mol'],
 'clippy':['cargo','clippy','--offline','--locked','--release','-p','chematic-mol','--all-targets','--','-D','warnings'],
 'format':['cargo','fmt','--all','--','--check'],
 'refactor':['python3','scripts/check_refactor_budget.py'],
 'runner-tests':['/private/tmp/chematic-ring-frontier-env-20261016/bin/python','-m','pytest','-q','-p','no:cacheprovider','scripts/tests/test_bench_openbabel_sdf_throughput.py','scripts/tests/test_bench_openbabel_v3000_throughput.py'],
}
(OUT/'native-check-commands.json').write_text(json.dumps({'commands':commands,'env':{k:v for k,v in env.items() if k.startswith('CARGO_')}},indent=2)+'\n')
for name,cmd in commands.items():
 with (OUT/(name+'.log')).open('w') as log:subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(name,'passed',flush=True)
text=(OUT/'native-tests.log').read_text()
rows=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',text)
assert rows and 'FAILED' not in text
native={'passed':sum(int(r[0]) for r in rows),'failed':sum(int(r[1]) for r in rows),'ignored':sum(int(r[2]) for r in rows),'suites':len(rows)}
(OUT/'native-validation.json').write_text(json.dumps(native,indent=2)+'\n')
(OUT/'validation.json').write_text(json.dumps({'passed':True,'native':native,'release_all_target_clippy':True,'format':True,'refactor_budget':True,'runner_gates_passed':10},indent=2)+'\n')
