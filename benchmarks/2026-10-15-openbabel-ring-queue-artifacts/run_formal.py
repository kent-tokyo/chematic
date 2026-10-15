from pathlib import Path
import subprocess,json,hashlib
s=Path('/private/tmp/chematic-openbabel-ring-queue-source-20261015');r=Path('/private/tmp/chematic-openbabel-ring-queue-20261015');old=Path('/private/tmp/chematic-openbabel-v3000-20261014');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();source=r/'source-manifest.json';m=json.loads(source.read_text())
def check_source():
 for n,h in m['frozen_build_source_sha256'].items():assert sha(s/n)==h,n
check_source()
for profile in ['0','5000','chains']:
 output=r/f'paired21-{profile}.json'
 cmd=['/private/tmp/chematic-cycle-speed-20261014/env/bin/python',str(r/'run-native-safe.py'),'--baseline',str(r/'baseline-v3000'),'--candidate',str(r/'candidate-v3000'),'--openbabel',str(old/'openbabel-v3000'),'--corpus',str(old/f'corpus-{profile}.molblocks'),'--records','64' if profile=='chains' else '1000','--blocks','21','--repeats','5','--source-manifest',str(source),'--output',str(output)]
 with (r/f'paired21-{profile}.log').open('w') as f:subprocess.run(cmd,cwd=s,stdout=f,stderr=subprocess.STDOUT,check=True)
 check_source();v=json.loads(output.read_text());assert v['complete'] and len(v['blocks'])==21 and v['source_parse_snapshots_equal'];a=output.with_suffix('.artifacts');assert (a/'baseline.molblocks').read_bytes()==(a/'candidate.molblocks').read_bytes()
 print(profile,'complete', {lane:cmps['baseline']['speedup'] for lane,cmps in v['summary'].items()},flush=True)
