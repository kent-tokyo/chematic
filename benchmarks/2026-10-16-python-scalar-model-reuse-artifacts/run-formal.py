from pathlib import Path
import hashlib, json, subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
PYTHON='/private/tmp/chematic-ring-frontier-env-20261016/bin/python'
plan=json.loads((OUT/'measurement-plan.json').read_text())
manifest=json.loads((OUT/'source-manifest.json').read_text())
comparators=json.loads((OUT/'comparator-manifest.json').read_text())
proof=json.loads((OUT/'parity.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def integrity():
 for p,expected in manifest['candidate_sha256'].items():assert sha(ROOT/p)==expected,p
 for item in comparators.values():
  for p,expected in item['native_files'].items():assert sha(Path(p))==expected,p
 for arm,item in proof['wheels'].items():
  assert sha(Path(item['wheel']))==item['sha256']
  for p,expected in item['extensions'].items():assert sha(OUT/arm/p)==expected,p
assert proof['outputs']['baseline']==proof['outputs']['candidate']
assert '0 failed' in (OUT/'native-tests.log').read_text()
assert 'error:' not in (OUT/'clippy.log').read_text()
assert 'Refactor budget OK' in (OUT/'refactor.log').read_text()
assert not (OUT/'format.log').read_bytes()
assert 'failed' not in (OUT/'python-tests.log').read_text().splitlines()[-1]
for offset in plan['offsets']:
 for lane in ['frontier','bundle']:
  integrity()
  cmd=[PYTHON,str(OUT/'run-python-safe.py')]
  if lane=='bundle':cmd+=['--bundle']
  cmd+=['--corpus',str(ROOT/plan['corpus']),'--offset',str(offset),'--limit',str(plan['rows_each']),'--blocks',str(plan['blocks']),'--base-revision',plan['base'],'--source-manifest',str(OUT/'source-manifest.json'),'--output',str(OUT/f'{lane}-{offset}.json')]
  for arm in ('baseline','candidate'):
   cmd+=['--'+arm+'-wheel',proof['wheels'][arm]['wheel'],'--'+arm+'-path',str(OUT/arm)]
  with (OUT/f'{lane}-{offset}.log').open('w') as log:
   subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
  integrity()
  data=json.loads((OUT/f'{lane}-{offset}.json').read_text())
  assert data['complete'] and len(data['blocks'])==plan['blocks']
  for engine,value in comparators.items():assert data['blocks'][0][engine]['version']==value['runtime_version']
  print(lane,offset,'complete',flush=True)
