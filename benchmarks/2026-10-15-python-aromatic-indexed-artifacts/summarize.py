from pathlib import Path
import json
r=Path('/private/tmp/chematic-parity-lazy-20261015')
for p in sorted(list(r.glob('frontier-*.json'))+list(r.glob('focused-*.json'))):
 j=json.loads(p.read_text());print(p.name,'blocks',len(j['blocks']),'complete',j['complete'])
 if not j['complete']:continue
 for k,v in j['summary'].items():
  x=v['baseline']['speedup'];print(k, 'base', None if not x else (round(x['median'],4),[round(a,4) for a in x['ci95']],x['winning_blocks']), 'comparators', {n:None if not v[n]['speedup'] else round(v[n]['speedup']['median'],3) for n in ['rdkit','cosmolkit']})
