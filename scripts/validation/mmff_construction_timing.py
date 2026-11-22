#!/usr/bin/env python3
"""Measure repeated MMFF field construction on the exposed COD development set."""
import gzip
import json
from pathlib import Path
import platform
import statistics
import time

import chematic

with gzip.open('validation/results/cod-linux-development-v8-2026-11-22.json.gz', 'rt') as f:
    source = json.load(f)
results = []
for row in source['results']:
    if row.get('status') != 'ok':
        continue
    item = {'name': row['name']}
    try:
        molecule = chematic.from_smiles(row['smiles']).add_hydrogens()
        n = len(molecule.mmff94_numeric_atom_types())
        coords = [[((i * 37) % 101)/17, ((i * 53) % 103)/19, ((i * 67) % 107)/23] for i in range(n)]
        molecule.mmff_force_field(coords)
        samples = []
        for _ in range(7):
            start = time.perf_counter_ns()
            molecule.mmff_force_field(coords)
            samples.append(time.perf_counter_ns() - start)
        item['median_ns'] = statistics.median(samples)
        item['samples_ns'] = samples
    except Exception as e:
        item['error'] = repr(e)
    results.append(item)
output = {'platform': platform.platform(), 'sample': 'exposed COD development',
          'repeats': 7, 'results': results}
Path(__import__('sys').argv[1]).write_text(json.dumps(output, indent=2)+'\n')
print(json.dumps({'rows': len(results), 'errors': sum('error' in x for x in results),
                  'sum_medians_ms': sum(x.get('median_ns',0) for x in results)/1e6}))
