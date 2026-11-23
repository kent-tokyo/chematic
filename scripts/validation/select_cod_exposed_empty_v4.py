#!/usr/bin/env python3
"""Select the already-exposed Linux v4 empty-ensemble cases for development."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

prior = Path('validation/results/cod-second-holdout-linux-fixed-v4-2026-11-22.json.gz')
split_path = Path('validation/results/cod-round2-scaffold-split-2026-11-22.json.gz')
assert hashlib.sha256(prior.read_bytes()).hexdigest() == 'f4b752f0557fe5d3da99f5dd0da4aa52cb0ab0b035b1e708502103d787b6d5cf'
assert hashlib.sha256(split_path.read_bytes()).hexdigest() == 'eaa5dd432a3c22320bbb27d0e7d62600b9848d52d76fa59f6e4e87a7b73b91fd'
prior_data = json.loads(gzip.decompress(prior.read_bytes()))
split = json.loads(gzip.decompress(split_path.read_bytes()))
assert prior_data['source_sha256'] == split['source_sha256']
selected = {r['index']: r['name'] for r in prior_data['results']
            if r.get('status') == 'ok' and r.get('chematic_count') == 0 and r.get('rdkit_count', 0) > 0}
assert len(selected) == 52
rows = []
for row in split['rows']:
    if row['index'] in selected:
        assert row['name'] == selected[row['index']]
        assert row['split'] == 'holdout'
        rows.append({**row, 'split': 'development'})
assert len(rows) == 52 and len({r['index'] for r in rows}) == 52
out = {'source_sha256': split['source_sha256'], 'rows': rows,
       'selection': 'exposed prior Linux v4 holdout rows with zero chematic and nonzero RDKit candidates; relabeled development',
       'prior_result_gzip_sha256': hashlib.sha256(prior.read_bytes()).hexdigest(),
       'original_split_gzip_sha256': hashlib.sha256(split_path.read_bytes()).hexdigest()}
Path(sys.argv[1]).write_text(json.dumps(out, indent=2) + '\n')
print('selected', len(rows), 'sha256', hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
