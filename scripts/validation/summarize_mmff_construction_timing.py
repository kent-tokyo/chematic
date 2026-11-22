#!/usr/bin/env python3
"""Compare A-B-B-A constructor runs without dropping failures."""
import json
from pathlib import Path
import statistics
import sys

files = [json.loads(Path(p).read_text()) for p in sys.argv[1:5]]
assert len(files) == 4
rows = [{r['name']: r for r in f['results']} for f in files]
assert all(len(x) == len(rows[0]) and x.keys() == rows[0].keys() for x in rows)
errors = {name: [r[name].get('error') for r in rows] for name in rows[0] if any('error' in r[name] for r in rows)}
results = []
for name in rows[0]:
    if name in errors:
        continue
    a = statistics.median(rows[0][name]['samples_ns'] + rows[3][name]['samples_ns'])
    b = statistics.median(rows[1][name]['samples_ns'] + rows[2][name]['samples_ns'])
    results.append({'name': name, 'v8_ns': a, 'v9_ns': b, 'ratio': b/a})
ratios = sorted(x['ratio'] for x in results)
summary = {'records': len(rows[0]), 'errors': errors, 'paired_records': len(results),
           'ratio_p50': statistics.median(ratios), 'ratio_p95': ratios[min(len(ratios)-1, int(len(ratios)*.95))],
           'sum_v8_medians_ms': sum(x['v8_ns'] for x in results)/1e6,
           'sum_v9_medians_ms': sum(x['v9_ns'] for x in results)/1e6}
out = {'classification': 'development-only synthetic-coordinate constructor microbenchmark',
       'run_order': ['v8','v9','v9','v8'], 'summary': summary, 'results': results}
Path(sys.argv[5]).write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(summary))
assert not errors
